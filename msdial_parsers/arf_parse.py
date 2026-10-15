"""MS-DIAL ``.arf``（アライン後のサンプル別ピーク）の単体パーサ。

``AlignmentResult_<timestamp>_PeakProperties.arf`` は、アラインメントの各スポット
（= 1 特徴量）について、全サンプル分の ``AlignmentChromPeakFeature`` 行を持つ。
サンプル間比較（正規化・PCA・差次的解析）の主データ源はこれ。同じフォルダの
``*_DriftSopts.arf`` はイオンモビリティ用で、通常は空（数十バイト）。

ファイル構造:
  トップレベルの MessagePack ストリームに 1 個以上の ``[header, payload]`` が並ぶ。
  payload は「展開後サイズ（msgpack int）＋ LZ4 ブロック」。展開すると、
  スポットのグループ（= サンプル行のリスト）を内包するコンテナ配列が現れる。
  1 ブロックが展開後 1 GB 近くになることがあるので、グループ単位で逐次に解く。

行（AlignmentChromPeakFeature）の Key 番号:
  0 FileID / 1 FileName / 2 MasterPeakID（負 = gap-fill）/ 3 PeakID / 9 MS2RawSpectrumID /
  10 MS2RawSpectrumID2CE（非空 dict = MS/MS 取得済み）/ 15 ChromXs（RT 等）/
  18 PeakHeightTop / 20 PeakAreaAboveZero / 21 PeakAreaAboveBaseline / 22 Mass /
  23 IonMode / 24 Name / 37 PeakShape = [EstimatedNoise, SignalToNoise, ...]

使い方（Python）::

    from msdial_parsers.arf_parse import read_arf, peak_table, intensity_matrix
    spots = read_arf("AlignmentResult_..._PeakProperties.arf")
    rows = peak_table(spots)                          # スポット × サンプルの縦持ち
    matrix, samples, spot_ids = intensity_matrix(spots, prop="height")

使い方（CLI）::

    python -m msdial_parsers arf <path.arf> [--csv peaks.csv] [--wide-csv matrix.csv --prop height]

注意:
  - gap-fill（``is_gap_filled=True``）は未検出サンプルに MS-DIAL が補間した値で、
    実検出ではない。検出率を数えるときは除く。
  - スポットの ``Name`` は代表サンプルの行から取ったもの。同定の正準は ``.arf2``。

依存: msgpack, lz4（行列化のみ numpy）。切り出し元: metabolomix/arf/reader.py。
"""

from __future__ import annotations

import io
import sys
from enum import Enum

import lz4.block
import msgpack


class IonMode(Enum):
    Positive = 0
    Negative = 1
    Both = 2
    Unknown = -1

    @classmethod
    def from_int(cls, val):
        try:
            return cls(val).name
        except Exception:
            return cls.Unknown.name


class ArfMemoryError(MemoryError):
    """読み切る前に RAM が尽きたことと、そこまでに読めた規模を運ぶ。"""

    def __init__(self, *, spots: int, rows: int):
        super().__init__(f"{spots} スポット / {rows} 行を読んだところで RAM が尽きた")
        self.spots = spots
        self.rows = rows


# --------------------------------------------------------------------------
# 低水準: MessagePack + LZ4 の展開
# --------------------------------------------------------------------------

def _decode(value):
    return value.decode("utf-8", errors="ignore") if isinstance(value, bytes) else value


def _unpacker(data: bytes) -> msgpack.Unpacker:
    return msgpack.Unpacker(
        io.BytesIO(data),
        raw=False,
        # POS モードの .arf はヘッダ相当フィールドに不正な UTF-8 を含むことがある。
        # strict だと 1 バイトで全体が止まるので、不正バイトは無視して読み進める。
        unicode_errors="ignore",
        strict_map_key=False,
        max_buffer_size=1024 * 1024 * 1024,
    )


def _decompress_lz4_payload(compressed: bytes) -> bytes:
    unpacker = _unpacker(compressed)
    size = next(unpacker)
    return lz4.block.decompress(compressed[unpacker.tell():], uncompressed_size=size)


def _payload_from_top_level_object(obj):
    if isinstance(obj, msgpack.ExtType):
        return obj.data
    if isinstance(obj, (list, tuple)) and len(obj) == 2:
        _header, payload = obj
        if isinstance(payload, msgpack.ExtType):
            return payload.data
        if isinstance(payload, (bytes, bytearray)):
            return bytes(payload)
    return None


def _iter_decompressed_blocks(data: bytes):
    """トップレベルの各ブロックを展開済みバイト列のまま返す（msgpack は解かない）。"""
    for block_index, obj in enumerate(_unpacker(data)):
        payload = _payload_from_top_level_object(obj)
        obj = None  # 圧縮ペイロードの写しを展開後まで抱えない
        if payload is None:
            print(f"[WARNING] top-level object {block_index} is not an LZ4 MessagePack block; skipped",
                  file=sys.stderr)
            continue
        decompressed = _decompress_lz4_payload(payload)
        payload = None
        yield block_index, decompressed


def _is_peak_group(item) -> bool:
    """item が 1 スポット分のグループ（サンプル行のリスト）か。"""
    return isinstance(item, list) and len(item) > 0 and all(isinstance(r, list) for r in item)


def _starts_array(first_byte: int) -> bool:
    """msgpack の型バイトが配列（fixarray / array16 / array32）か。"""
    return 0x90 <= first_byte <= 0x9F or first_byte in (0xDC, 0xDD)


def _iter_peak_groups_in_stream(decompressed: bytes):
    """展開済みストリームからグループを 1 個ずつ返す。

    - グループ要素を 1 つ以上含む項目はコンテナ。グループ要素だけを返す。
    - グループ要素を含まず、要素が全部行なら、項目自体が単独グループ。
    - それ以外（メタデータ・配列でない項目）は読み飛ばす。
    """
    unpacker = _unpacker(decompressed)
    end = len(decompressed)
    try:
        while unpacker.tell() < end:
            if not _starts_array(decompressed[unpacker.tell()]):
                unpacker.skip()
                continue
            length = unpacker.read_array_header()
            found_group = False
            others = []
            for _ in range(length):
                element = unpacker.unpack()
                if _is_peak_group(element):
                    found_group = True
                    yield element
                elif not found_group:
                    others.append(element)
            if not found_group and _is_peak_group(others):
                yield others
    except msgpack.OutOfData:
        return  # 末尾が切れている。そこまでのグループは返し済み。


def iter_peak_groups(data: bytes):
    """``(block_index, local_index, group)`` を順に返す。group は生の行リスト。"""
    for block_index, decompressed in _iter_decompressed_blocks(data):
        local_index = 0
        for group in _iter_peak_groups_in_stream(decompressed):
            yield block_index, local_index, group
            local_index += 1


# --------------------------------------------------------------------------
# 行（AlignmentChromPeakFeature）の解釈
# --------------------------------------------------------------------------

def _convert_to_times(data: list) -> dict:
    result = {}
    for t in data[:4]:
        if not (isinstance(t, (list, tuple)) and len(t) >= 2):
            continue
        if isinstance(t[1], (list, tuple)) and len(t[1]) > 0:
            if t[0] == 1:
                result["rt"] = t[1][0]
            elif t[0] == 2:
                result["ri"] = t[1][0]
            elif t[0] == 3:
                result["m/z"] = t[1][0]
            elif t[0] == 4:
                result["dt"] = t[1][0]
    return result


def _rt_of(row: list):
    if len(row) > 15 and isinstance(row[15], list):
        return _convert_to_times(row[15]).get("rt")
    if len(row) > 16 and isinstance(row[16], list):
        return _convert_to_times(row[16]).get("rt")
    return None


def _is_alignment_peak_row(row) -> bool:
    return (isinstance(row, list) and len(row) > 25
            and isinstance(row[18], (int, float)))


def _file_name_from_row(row: list) -> str | None:
    for item in row[:10]:
        if isinstance(item, (bytes, str)):
            decoded = _decode(item)
            if decoded:
                return decoded
    return None


def peak_row(row: list) -> dict:
    """``AlignedPeakProperties`` の 1 行（1 サンプル分）を辞書にする。読めない行は空 dict。"""
    if not _is_alignment_peak_row(row):
        return {}
    peak_shape = row[37] if len(row) > 37 and isinstance(row[37], list) else []
    master_peak_id = row[2]
    ms2_ce = row[10]
    return {
        "file_id": row[0],
        "file_name": _file_name_from_row(row),
        "peak_id": row[3],
        "master_peak_id": master_peak_id,
        "height": row[18],
        "area": row[20],
        "area_above_baseline": row[21],
        "m_z": row[22],
        "rt": _rt_of(row),
        "signal_to_noise": peak_shape[1] if len(peak_shape) > 1 else None,
        "estimated_noise": peak_shape[0] if len(peak_shape) > 0 else None,
        "ms2_raw_id": row[9],
        "is_msms": bool(ms2_ce) if isinstance(ms2_ce, dict) else False,
        # MasterPeakID が負（-2）のサンプルは gap-fill（未検出 → 補間値）
        "is_gap_filled": isinstance(master_peak_id, int) and master_peak_id < 0,
    }


#: 行を細くするときに残す Key。15（rt）と 37（PeakShape）は読む部分だけに縮めて残す。
#: **peak_row で新しい Key を読むときは、ここにも足すこと**（足さないと None になる）。
_RETAINED_ROW_KEYS = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 18, 20, 21, 22)


def _slim_row(row, shared_names: dict):
    """行を peak_row が同じ答えを返す範囲で細くする（生の行は 1 行 6.5 KB ある）。"""
    if not isinstance(row, list):
        return row
    slim = [None] * len(row)  # 長さは保つ（len(row) > 25 の判定を変えない）
    for key in _RETAINED_ROW_KEYS:
        if key < len(row):
            slim[key] = row[key]
    if len(slim) > 1 and isinstance(slim[1], str):
        slim[1] = shared_names.setdefault(slim[1], slim[1])  # 同じファイル名 str を共有
    rt = _rt_of(row)
    if rt is not None:
        slim[15] = [[1, [rt]]]
    if len(row) > 37 and isinstance(row[37], list):
        slim[37] = row[37][:2]
    return slim


def _convert_group(group: list, index: int, block_index: int, local_index: int) -> dict | None:
    if not (isinstance(group, list) and group
            and all(_is_alignment_peak_row(row) for row in group)):
        return None
    rep = group[0]
    times = _convert_to_times(rep[15]) if len(rep) > 15 and isinstance(rep[15], list) else {}
    return {
        "MasterAlignmentID": index,
        "AlignmentID": index,
        "RT": times.get("rt"),
        "MassCenter": rep[22] if len(rep) > 22 else None,
        "IonMode": IonMode.from_int(rep[23]) if len(rep) > 23 else IonMode.Unknown.name,
        "Name": _decode(rep[24]) if len(rep) > 24 else None,
        "HeightAverage": rep[18] if len(rep) > 18 else None,
        "AlignedPeakProperties": group,
        "SourceBlockIndex": block_index,
        "SourceLocalIndex": local_index,
    }


# --------------------------------------------------------------------------
# 公開 API
# --------------------------------------------------------------------------

def deserialize(data: bytes, slim: bool = True) -> list[dict]:
    """``.arf`` のバイト列をスポット辞書のリストにする。

    各スポットは ``MasterAlignmentID`` / ``RT`` / ``MassCenter`` / ``IonMode`` /
    ``Name`` / ``HeightAverage`` と、生の行リスト ``AlignedPeakProperties`` を持つ。
    ``slim=True``（既定）は行から peak_row が読まない Key を落としてメモリを抑える。
    生の行を全部残したいとき（未知の Key を調べるとき）は ``slim=False``。
    """
    results = []
    rows = 0
    shared_names: dict[str, str] = {}
    try:
        for index, (block_index, local_index, group) in enumerate(iter_peak_groups(data)):
            spot = _convert_group(group, index, block_index, local_index)
            if spot is None:
                continue
            if slim:
                spot["AlignedPeakProperties"] = [_slim_row(r, shared_names) for r in group]
            results.append(spot)
            rows += len(group)
    except MemoryError as exc:
        raise ArfMemoryError(spots=len(results), rows=rows) from exc
    return results


def read_arf(path, slim: bool = True) -> list[dict]:
    """``.arf`` ファイルを読んでスポット辞書のリストを返す。"""
    with open(path, "rb") as f:
        return deserialize(f.read(), slim=slim)


def peak_table(spots: list[dict]) -> list[dict]:
    """スポット × サンプルの縦持ち表（1 行 = 1 スポットの 1 サンプル）を返す。"""
    table = []
    for spot in spots:
        aligned = spot.get("AlignedPeakProperties")
        if not isinstance(aligned, list):
            continue
        for sample_index, row in enumerate(aligned):
            if not isinstance(row, list) or len(row) < 3:
                continue
            p = peak_row(row)
            table.append({
                "MasterAlignmentID": spot.get("MasterAlignmentID"),
                "SpotRT": spot.get("RT"),
                "SpotMassCenter": spot.get("MassCenter"),
                "IonMode": spot.get("IonMode"),
                "CompoundName": spot.get("Name"),
                "SampleIndex": sample_index,
                "FileName": p.get("file_name") or f"Sample_{sample_index}",
                "FileID": p.get("file_id"),
                "PeakID": p.get("peak_id"),
                "MasterPeakID": p.get("master_peak_id"),
                "PeakHeight": p.get("height"),
                "PeakArea": p.get("area"),
                "PeakAreaAboveBaseline": p.get("area_above_baseline"),
                "PeakMZ": p.get("m_z"),
                "PeakRT": p.get("rt"),
                "SignalToNoise": p.get("signal_to_noise"),
                "IsMsms": p.get("is_msms"),
                "IsGapFilled": p.get("is_gap_filled"),
            })
    return table


def intensity_matrix(spots: list[dict], prop: str = "height", include_gap_filled: bool = True):
    """サンプル × スポットの強度行列を返す: ``(matrix, sample_names, spot_ids)``。

    ``prop`` は peak_row のキー（``height`` / ``area`` / ``area_above_baseline`` など）。
    値の無いセルは NaN。``include_gap_filled=False`` で gap-fill のセルも NaN にする。
    補完・足切り・正規化はしない（ここは読むだけ）。
    """
    try:
        import numpy as np
    except ImportError as exc:
        raise ImportError(
            "intensity_matrix() requires NumPy. Install it with "
            '`pip install "msdial-parsers[matrix]"`. '
            "When installing from a source checkout, use "
            '`pip install ".[matrix]"`.'
        ) from exc

    sample_index: dict[str, int] = {}
    columns: list[dict[str, float]] = []
    spot_ids = []
    for spot in spots:
        values: dict[str, float] = {}
        for i, row in enumerate(spot.get("AlignedPeakProperties") or []):
            p = peak_row(row)
            if not p:
                continue
            name = p["file_name"] or f"Sample_{i}"
            sample_index.setdefault(name, len(sample_index))
            if p["is_gap_filled"] and not include_gap_filled:
                continue
            v = p.get(prop)
            if isinstance(v, (int, float)):
                values[name] = float(v)
        columns.append(values)
        spot_ids.append(spot.get("MasterAlignmentID"))

    samples = list(sample_index)
    matrix = np.full((len(samples), len(columns)), np.nan)
    for j, values in enumerate(columns):
        for name, v in values.items():
            matrix[sample_index[name], j] = v
    return matrix, samples, spot_ids


def summarize(spots: list[dict]) -> dict:
    rts = [s["RT"] for s in spots if isinstance(s.get("RT"), (int, float))]
    mzs = [s["MassCenter"] for s in spots if isinstance(s.get("MassCenter"), (int, float))]
    modes: dict[str, int] = {}
    for s in spots:
        modes[s.get("IonMode")] = modes.get(s.get("IonMode"), 0) + 1
    samples = {peak_row(r).get("file_name") for s in spots[:1]
               for r in s.get("AlignedPeakProperties") or []}
    return {
        "total_spots": len(spots),
        "samples": len(samples - {None}),
        "rt_range": (min(rts), max(rts)) if rts else None,
        "mz_range": (min(mzs), max(mzs)) if mzs else None,
        "ion_modes": modes,
        "named_spots": sum(1 for s in spots if s.get("Name")),
    }


# Preserve the low-level API exposed by the original standalone script.
from .raw.arf_parse import decode, decompress, describe, iter_objects, write_json





