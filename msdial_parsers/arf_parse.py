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

    from msdial_parsers import read_arf
    spots = read_arf("AlignmentResult_..._PeakProperties.arf")

使い方（CLI）::

    python -m msdial_parsers arf <path.arf> [--csv peaks.csv]

注意:
  - gap-fill（``is_gap_filled=True``）は未検出サンプルに MS-DIAL が補間した値で、
    実検出ではない。検出率を数えるときは除く。
  - スポットの ``Name`` は代表サンプルの行から取ったもの。同定の正準は ``.arf2``。

依存: msgpack, lz4。切り出し元: metabolomix/arf/reader.py。
"""

from __future__ import annotations

import io
import sys
from enum import Enum

import lz4.block
import msgpack

from ._sources import open_binary_source


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


def _iter_peak_groups(data: bytes):
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


def _sample_peak(row: list, *, slim: bool) -> dict:
    """Convert one AlignmentChromPeakFeature row to named Python fields."""
    times = _convert_to_times(row[15]) if len(row) > 15 and isinstance(row[15], list) else {}

    def optional_time(name):
        value = times.get(name)
        return None if value is None or value == -1 else value

    record = {
        "file_id": row[0] if len(row) > 0 else None,
        "file_name": _decode(row[1]) if len(row) > 1 else None,
        "master_peak_id": row[2] if len(row) > 2 else None,
        "peak_id": row[3] if len(row) > 3 else None,
        "rt": optional_time("rt"),
        "ri": optional_time("ri"),
        "drift_time": optional_time("dt"),
        "m/z": row[22] if len(row) > 22 else None,
        "peak_height_top": row[18] if len(row) > 18 else None,
        "peak_area_above_zero": row[20] if len(row) > 20 else None,
        "peak_area_above_baseline": row[21] if len(row) > 21 else None,
        "ion_mode": IonMode.from_int(row[23]) if len(row) > 23 else IonMode.Unknown.name,
        "name": _decode(row[24]) if len(row) > 24 else None,
        "formula": _decode(row[25][0]) if len(row) > 25 and isinstance(row[25], list) and row[25] else None,
    }
    if slim:
        return record

    record.update({
        "ms2_raw_spectrum_id": row[9] if len(row) > 9 else None,
        "ms2_raw_spectrum_id_2ce": row[10] if len(row) > 10 else None,
        "ontology": _decode(row[26]) if len(row) > 26 else None,
        "smiles": _decode(row[27]) if len(row) > 27 else None,
        "inchi_key": _decode(row[28]) if len(row) > 28 else None,
        "adduct_type": row[36] if len(row) > 36 else None,
        "peak_shape": row[37] if len(row) > 37 else None,
    })
    return record


def _convert_group(group: list, index: int) -> dict | None:
    if not (isinstance(group, list) and group
            and all(_is_alignment_peak_row(row) for row in group)):
        return None
    return {
        "spot_index": index,
        "sample_peaks": [_sample_peak(row, slim=True) for row in group],
    }


def _deserialize(data: bytes, slim: bool = True) -> list[dict]:
    """Convert each spot group into a spot index and named sample peak rows."""
    results = []
    rows = 0
    try:
        for index, (_block_index, _local_index, group) in enumerate(_iter_peak_groups(data)):
            spot = _convert_group(group, index)
            if spot is None:
                continue
            if not slim:
                spot["sample_peaks"] = [_sample_peak(row, slim=False) for row in group]
            results.append(spot)
            rows += len(group)
    except MemoryError as exc:
        raise ArfMemoryError(spots=len(results), rows=rows) from exc
    return results


def read_arf(source, slim: bool = True) -> list[dict]:
    """``.arf`` ファイルを読んでスポット辞書のリストを返す。"""
    with open_binary_source(source) as f:
        return _deserialize(f.read(), slim=slim)
