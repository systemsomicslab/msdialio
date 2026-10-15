"""MS-DIAL ``.pai2``（1 測定ファイル分のピーク一覧）の単体パーサ。

``<測定名>_<timestamp>.pai2`` は、1 サンプルを個別にピーク検出した結果
（``ChromatogramPeakFeature`` の並び）。**サンプル間比較はできない**——
複数サンプルの比較はアライン後の ``.arf``（arf_parse.py）を使う。

MS/MS の実スペクトルは同名の ``.dcl`` にある。``.pai2`` の ``has_msms`` は
「取得の参照がある」ことしか示さない。実スペクトルを付けるには
``dcl_parse.attach_msms_to_features(peaks, dcl_parse.read_dcl(<同名.dcl>))``。
``.dcl`` の並び順（dcl_index）が ``.pai2`` のピーク順に対応する。

ファイル構造:
  ``[header, payload]`` の MessagePack 1 個。payload は「展開後サイズ ＋ LZ4 ブロック」。
  展開後は ChromatogramPeakFeature 配列のリスト。

Key 番号（ChromatogramPeakFeature）:
  3/4/5 ChromXsLeft/Top/Right / 6/7/8 PeakHeightLeft/Top/Right / 9 PeakAreaAboveZero /
  10 PeakAreaAboveBaseline / 11 PeakID / 18 MS2RawSpectrumID / 19 MS2RawSpectrumID2CE /
  22 IonMode / 24 Spectrum / 25 Name / 26 Formula / 27 Ontology / 28 SMILES / 29 InChIKey /
  30 AdductType / 31 CollisionCrossSection / 38 Comment / 40 PeakShape（[1] = S/N）/ 43 Mass

使い方（Python）::

    from msdial_parsers.pai2_parse import read_pai2
    peaks = read_pai2("sample_202605141515.pai2")

使い方（CLI）::

    python -m msdial_parsers pai2 <path.pai2> [--csv peaks.csv]

依存: msgpack, lz4。切り出し元: metabolomix/pai2/reader.py。
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


def deserialize_lz4_packed_msgpack(data: bytes):
    _header, compressed = msgpack.unpackb(data, raw=False)
    unpacker = msgpack.Unpacker(io.BytesIO(compressed), raw=False)
    size = next(unpacker)
    decompressed = lz4.block.decompress(compressed[unpacker.tell():], uncompressed_size=size)
    return msgpack.unpackb(decompressed, raw=False, strict_map_key=False)


def _convert_to_times(data) -> dict:
    """ChromXs を ``{"rt", "ri", "m/z", "dt"}`` にする（負値 = 未設定は落とす）。"""
    result = {}
    if not isinstance(data, list):
        return result
    for t in data[:4]:
        if isinstance(t, (list, tuple)) and len(t) >= 2:
            vals = t[1]
            if isinstance(vals, (list, tuple)) and len(vals) > 0 and vals[0] >= 0:
                key = {1: "rt", 2: "ri", 3: "m/z", 4: "dt"}.get(t[0])
                if key:
                    result[key] = vals[0]
    return result


def peak_record(data: list) -> dict:
    """生の ChromatogramPeakFeature 1 個をピーク辞書にする。"""
    formula_data = data[26]
    formula = formula_data[0] if isinstance(formula_data, (list, tuple)) and formula_data else "Unknown"
    adduct_data = data[30]
    adduct = adduct_data[2] if isinstance(adduct_data, (list, tuple)) and len(adduct_data) > 2 else "Unknown"

    # MS/MS: Key18 MS2RawSpectrumID / Key19 MS2RawSpectrumID2CE / Key24 Spectrum。
    # 実スペクトル本体は .dcl にあり、.pai2 の Key24 は通常空。
    ms2_id = data[18] if len(data) > 18 and isinstance(data[18], int) else -1
    ce_map = data[19] if len(data) > 19 and isinstance(data[19], dict) else {}
    spectrum = data[24] if len(data) > 24 and isinstance(data[24], list) else []
    msms_spectrum = [[p[0], p[1]] for p in spectrum if isinstance(p, (list, tuple)) and len(p) >= 2]

    return {
        "id": data[11],
        "name": data[25],
        "m/z": data[43],
        "time": _convert_to_times(data[4]),
        "time_left": _convert_to_times(data[3]),
        "time_right": _convert_to_times(data[5]),
        "peak_height": data[7],
        "peak_height_left": data[6],
        "peak_height_right": data[8],
        "peak_area": data[9],
        "peak_area_above_baseline": data[10],
        "S/N": data[40][1],
        "ion_mode": IonMode.from_int(data[22]),
        "formula": formula,
        "ontology": data[27],
        "smiles": data[28],
        "inchikey": data[29],
        "adduct": adduct,
        "collision_cross_section": data[31],
        "comment": data[38],
        "has_msms": bool(ce_map) or ms2_id >= 0,
        "ms2_raw_id": ms2_id,
        "collision_energies": sorted({float(v) for v in ce_map.values()}) if ce_map else [],
        "msms_peak_count": len(msms_spectrum),
        "msms_spectrum": msms_spectrum,
    }


def deserialize(data: bytes) -> list[dict]:
    return [peak_record(d) for d in deserialize_lz4_packed_msgpack(data)]


def read_pai2(path) -> list[dict]:
    with open(path, "rb") as f:
        return deserialize(f.read())


def flat_record(peak: dict) -> dict:
    """CSV 向けに 1 段へ平らにする（time を rt/ri へ展開し、スペクトル列は外す）。"""
    out = {k: v for k, v in peak.items()
           if k not in ("time", "time_left", "time_right", "msms_spectrum", "collision_energies")}
    out["rt"] = peak["time"].get("rt")
    out["rt_left"] = peak["time_left"].get("rt")
    out["rt_right"] = peak["time_right"].get("rt")
    out["collision_energies"] = ";".join(str(e) for e in peak["collision_energies"])
    return out


def _is_annotated(name) -> bool:
    if not isinstance(name, str):
        return False
    s = name.strip().lower()
    return bool(s) and s != "unknown" and not s.startswith(("no ms2", "low score"))


def summarize(peaks: list[dict]) -> dict:
    if not peaks:
        return {"total_peaks": 0}
    mzs = [p["m/z"] for p in peaks if isinstance(p["m/z"], (int, float))]
    rts = [p["time"]["rt"] for p in peaks if isinstance(p["time"].get("rt"), (int, float))]
    modes: dict[str, int] = {}
    for p in peaks:
        modes[p["ion_mode"]] = modes.get(p["ion_mode"], 0) + 1
    annotated = sum(1 for p in peaks if _is_annotated(p["name"]))
    return {
        "total_peaks": len(peaks),
        "annotated": annotated,
        "with_msms_reference": sum(1 for p in peaks if p["has_msms"]),
        "mz_range": (min(mzs), max(mzs)) if mzs else None,
        "rt_range": (min(rts), max(rts)) if rts else None,
        "ion_modes": modes,
    }


# Preserve the raw decoder API from the original standalone script.
from .raw.pai2_parse import decode, decompress, describe, iter_objects, write_json




