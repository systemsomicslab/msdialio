"""MS-DIAL ``.arf2``（アラインメントのスポット代表カタログ）の単体パーサ。

``AlignmentResult_<timestamp>.arf2`` は、各スポットの代表値（m/z・RT・平均強度）と
同定情報（名前・組成式・Ontology・SMILES・InChIKey・付加体）と品質指標を持つ。
**サンプル別の強度は持たない**——サンプル間比較には ``.arf``（arf_parse.py）を使い、
``MasterAlignmentID`` で突き合わせる。

ファイル構造:
  ``[header, payload]`` の MessagePack 1 個。payload は「展開後サイズ ＋ LZ4 ブロック」。
  展開後は ``AlignmentSpotProperty`` 配列の並び（入れ子のリストに包まれていることもある）。

Key 番号（AlignmentSpotProperty）:
  0 MasterAlignmentID / 1 AlignmentID / 4 TimesCenter / 5 MassCenter / 11 IonMode /
  12 Name / 13 Formula / 14 Ontology / 15 SMILES / 16 InChIKey / 31 HeightAverage /
  32-33 HeightMin/Max / 34 PeakWidthAverage / 35-37 SignalToNoise Ave/Max/Min /
  43-44 MassMin/Max / 49 FillPercentage / 51 MonoIsotopicPercentage /
  53 IsotopicPeaks / 54 AdductType

使い方（Python）::

    from msdial_parsers.arf2_parse import read_arf2
    spots = read_arf2("AlignmentResult_....arf2")

使い方（CLI）::

    python -m msdial_parsers arf2 <path.arf2> [--csv spots.csv]

注意: ``Name`` が付いていても確定同定とは限らない（``Unknown`` / ``no MS2:`` /
``low score:`` などの接頭辞を見る）。

依存: msgpack, lz4。切り出し元: metabolomix/arf2/reader.py。
"""

from __future__ import annotations

import io
from enum import Enum

import lz4.block
import msgpack


class IonMode(Enum):
    Positive = 0
    Negative = 1
    Both = 2
    Unknown = -1

    @classmethod
    def from_int(cls, val) -> str:
        try:
            return cls(val).name
        except ValueError:
            return cls.Unknown.name


def _decode(val) -> str:
    if isinstance(val, bytes):
        return val.decode("utf-8", errors="ignore")
    if isinstance(val, str):
        return val
    return ""


def _to_float(val, default: float = 0.0) -> float:
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


def _nested_at(val, idx: int, default=None):
    if isinstance(val, (list, tuple)) and len(val) > idx:
        return val[idx]
    return default


def _extract_times(data) -> dict:
    result = {}
    if not isinstance(data, list):
        return result
    for t in data[:4]:
        if isinstance(t, list) and len(t) >= 2 and isinstance(t[1], (list, tuple)) and len(t[1]) > 0:
            if t[0] == 1:
                result["rt"] = float(t[1][0])
            elif t[0] == 2:
                result["ri"] = float(t[1][0])
            elif t[0] == 3:
                result["m_z"] = float(t[1][0])
            elif t[0] == 4:
                result["dt"] = float(t[1][0])
    return result


def deserialize_lz4_packed_msgpack(data: bytes) -> list:
    """``[header, payload]`` を展開し、中の MessagePack オブジェクトを並びのまま返す。"""
    _header, compressed = msgpack.unpackb(data, raw=False)
    unpacker = msgpack.Unpacker(io.BytesIO(compressed), raw=False)
    size = next(unpacker)
    decompressed = lz4.block.decompress(compressed[unpacker.tell():], uncompressed_size=size)
    return list(msgpack.Unpacker(io.BytesIO(decompressed), raw=False, strict_map_key=False))


def iter_raw_spots(datas: list):
    """展開結果から AlignmentSpotProperty の生配列を順に返す。

    スポットを直接並べる形と、リストの入れ子に包む形の両方がある。
    """
    for d in datas:
        if isinstance(d, list) and len(d) < 14:
            for item in d:
                if isinstance(item, list) and len(item) > 0:
                    for spot_raw in item:
                        if isinstance(spot_raw, list) and len(spot_raw) >= 10:
                            yield spot_raw
            continue
        if isinstance(d, list) and len(d) >= 10:
            yield d


def spot_record(data: list) -> dict | None:
    """生の AlignmentSpotProperty 1 個をカタログ辞書にする。"""
    if not isinstance(data, list) or len(data) < 10:
        return None

    def get(index, default=None):
        return data[index] if len(data) > index else default

    times = _extract_times(get(4))
    return {
        "MasterAlignmentID": int(_to_float(get(0))),
        "AlignmentID": int(_to_float(get(1))),
        "RT": float(times.get("rt", 0.0)),
        "MassCenter": _to_float(get(5)),
        "IonMode": IonMode.from_int(get(11)),
        "Name": _decode(get(12)) or "Unknown",
        "HeightAverage": _to_float(get(31)),
        # 同定 / 化学情報
        "Formula": _decode(_nested_at(get(13), 0, "")) or "",   # Key13 = [組成式, 質量, ...]
        "Ontology": _decode(get(14)) or "",
        "SMILES": _decode(get(15)) or "",
        "InChIKey": _decode(get(16)) or "",
        "AdductType": _decode(_nested_at(get(54), 2, "")) or "",  # Key54 の index2 が '[M-H]-' 等
        # 強度 / 品質
        "HeightMin": _to_float(get(32)),
        "HeightMax": _to_float(get(33)),
        "PeakWidthAverage": _to_float(get(34)),
        "SignalToNoiseAve": _to_float(get(35)),
        "SignalToNoiseMax": _to_float(get(36)),
        "SignalToNoiseMin": _to_float(get(37)),
        "MassMin": _to_float(get(43)),
        "MassMax": _to_float(get(44)),
        "FillPercentage": _to_float(get(49)),      # 検出されたサンプルの割合
        "MonoIsotopicPercentage": _to_float(get(51)),
    }


def deserialize(data: bytes) -> list[dict]:
    """``.arf2`` のバイト列をスポット辞書のリストにする。"""
    return [r for r in (spot_record(s) for s in iter_raw_spots(deserialize_lz4_packed_msgpack(data))) if r]


def read_arf2(path) -> list[dict]:
    with open(path, "rb") as f:
        return deserialize(f.read())


def read_raw_spots(path) -> list[list]:
    """生の AlignmentSpotProperty 配列の一覧（Key 番号で自分で引くとき用）。"""
    with open(path, "rb") as f:
        return list(iter_raw_spots(deserialize_lz4_packed_msgpack(f.read())))


def read_isotopic_peaks(path) -> dict:
    """Key53 ``IsotopicPeaks``（代表試料の同位体ピークの生配列）を MasterAlignmentID ごとに返す。"""
    return {int(_to_float(raw[0])): (raw[53] if len(raw) > 53 else None)
            for raw in read_raw_spots(path)}


def _is_annotated(name: str) -> bool:
    s = (name or "").strip().lower()
    return bool(s) and s != "unknown" and not s.startswith(("no ms2", "low score"))


def summarize(spots: list[dict]) -> dict:
    if not spots:
        return {"total_spots": 0}
    rts = [s["RT"] for s in spots]
    mzs = [s["MassCenter"] for s in spots]
    modes: dict[str, int] = {}
    ontology: dict[str, int] = {}
    for s in spots:
        modes[s["IonMode"]] = modes.get(s["IonMode"], 0) + 1
        if s["Ontology"].strip():
            ontology[s["Ontology"]] = ontology.get(s["Ontology"], 0) + 1
    annotated = sum(1 for s in spots if _is_annotated(s["Name"]))
    return {
        "total_spots": len(spots),
        "annotated": annotated,
        "annotation_rate_pct": round(annotated / len(spots) * 100, 1),
        "rt_range": (min(rts), max(rts)),
        "mz_range": (min(mzs), max(mzs)),
        "ion_modes": modes,
        "ontology_top": dict(sorted(ontology.items(), key=lambda kv: -kv[1])[:10]),
    }




