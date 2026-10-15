"""MS-DIAL ``.EIC.aef``（アラインメントの EIC クロマトグラム, CSS1 形式）の単体パーサ。

``AlignmentResult_<timestamp>.EIC.aef`` は、各スポット × 各サンプルの抽出イオン
クロマトグラム（横軸 RT・縦軸 強度の点列）を持つ独自バイナリ（リトルエンディアン）。
スポット順は ``.arf`` / ``.arf2`` の MasterAlignmentID 順に対応する。

バイト構造:
  ヘッダ: "CSS1"(4) + 6B（未使用）+ スポット数 int32  → 計 14B
  ポインタ表: スポット数 × int64（各スポットの開始位置）
  スポット: RT f32, RI f32, Mass f32, Drift f32, MainType i8, サンプル数 i32（21B）
    サンプル: FileID i32, 点数 i32, PeakTop f32, PeakLeft f32, PeakRight f32（20B）
      点列: 点数 × (横軸 f32, 強度 f32)

注意: ``peak_top`` / ``peak_left`` / ``peak_right`` は**横軸（RT）座標**で強度ではない。
強度で並べるときは点列の最大値（``max_intensity``）を使う。

使い方（Python）::

    from msdial_parsers.aef_parse import read_aef, read_aef_spots
    spots = read_aef("AlignmentResult_....EIC.aef")              # 全スポットの要約（点列なし）
    traces = read_aef_spots("AlignmentResult_....EIC.aef", [12, 345], file_ids=[0, 1])

使い方（CLI）::

    python -m msdial_parsers aef <path.EIC.aef> [--csv spots.csv] [--spot 12 --trace-csv eic.csv]

依存: 標準ライブラリのみ。切り出し元: metabolomix/eic/reader.py。
"""

from __future__ import annotations

import os
import struct

MAGIC = b"CSS1"
_HEADER_SIZE = 14  # magic 4 + 未使用 6 + スポット数 4


def read_aef(path, include_chromatogram: bool = False) -> list[dict]:
    """全スポットを順に読む。

    各サンプルは ``file_id`` / ``peak_top`` / ``num_peaks``（点数）/ ``mean_intensity`` /
    ``max_intensity``。``include_chromatogram=True`` で ``chromatogram``（(rt, 強度) の列）
    も付ける——全スポット分だとメモリを大きく食うので、点列が要るなら read_aef_spots を使う。
    """
    results = []
    with open(path, "rb") as f:
        magic = f.read(4)
        if magic != MAGIC:
            raise ValueError(f"Not a CSS1 file. Magic: {magic}")
        f.seek(10)
        num_spots = struct.unpack("<i", f.read(4))[0]
        pointers = struct.unpack(f"<{num_spots}q", f.read(8 * num_spots))

        for spot_id, ptr in enumerate(pointers):
            f.seek(ptr)
            rt, ri, mass, drift, main_type, num_samples = struct.unpack("<ffffbi", f.read(21))
            samples = []
            for _ in range(num_samples):
                file_id, num_points, top, _left, _right = struct.unpack("<iifff", f.read(20))
                raw = f.read(8 * num_points)
                points = list(struct.iter_unpack("<ff", raw))
                intensities = [p[1] for p in points]
                sample = {
                    "file_id": file_id,
                    "peak_top": top,
                    "num_peaks": num_points,
                    "mean_intensity": sum(intensities) / num_points if num_points else 0.0,
                    "max_intensity": max(intensities, default=0.0),
                }
                if include_chromatogram:
                    sample["chromatogram"] = points
                samples.append(sample)
            results.append({
                "spot_id": spot_id,
                "rt": rt,
                "ri": ri,
                "mz": mass,
                "drift": drift,
                "main_type": main_type,
                "num_samples": num_samples,
                "samples": samples,
            })
    return results


def _read_exact(stream, size, context):
    data = stream.read(size)
    if len(data) != size:
        raise ValueError(f"Unexpected end of EIC/AEF while reading {context}: "
                         f"expected {size} bytes, got {len(data)}")
    return data


def read_aef_spots(path, spot_ids, file_ids=None, *, max_total_points: int = 2_000_000,
                   strict: bool = False) -> list[dict]:
    """指定スポットだけをポインタ表から直接読み、点列（chromatogram）付きで返す。

    file_ids を渡すとそのサンプルだけを読み、他は点列を読み飛ばす。
    strict=True なら範囲外の spot_id と、見つからない file_id で例外にする。
    strict=False なら範囲外のスポットは結果から省き、file_id が無いスポットは
    ``samples`` を空にして返す（両者を区別できる）。
    """
    requested_set = set(file_ids) if file_ids is not None else None
    spots = []
    total_points = 0
    with open(path, "rb") as stream:
        file_size = os.fstat(stream.fileno()).st_size
        magic = _read_exact(stream, 4, "CSS1 magic")
        if magic != MAGIC:
            raise ValueError(f"Not a CSS1 file. Magic: {magic}")
        stream.seek(10)
        num_spots = struct.unpack("<i", _read_exact(stream, 4, "spot count"))[0]
        table_end = _HEADER_SIZE + 8 * num_spots

        for spot_id in sorted(set(spot_ids)):
            if not 0 <= spot_id < num_spots:
                if strict:
                    raise ValueError(f"spot_id {spot_id} is out of range for {num_spots} EIC spots")
                continue
            stream.seek(_HEADER_SIZE + 8 * spot_id)
            pointer = struct.unpack("<q", _read_exact(stream, 8, "spot pointer"))[0]
            if pointer < table_end or pointer >= file_size:
                raise ValueError(f"Invalid EIC spot pointer for spot_id {spot_id}: {pointer}")
            stream.seek(pointer)
            rt, ri, mass, drift, main_type, num_samples = struct.unpack(
                "<ffffbi", _read_exact(stream, 21, f"spot {spot_id} header"))

            samples = []
            found = set()
            for sample_index in range(num_samples):
                file_id, num_points, top, left, right = struct.unpack(
                    "<iifff", _read_exact(stream, 20, f"spot {spot_id} sample {sample_index} header"))
                point_bytes = 8 * num_points
                if requested_set is not None and file_id not in requested_set:
                    stream.seek(point_bytes, os.SEEK_CUR)
                    continue
                if total_points + num_points > max_total_points:
                    raise ValueError(f"Selected EIC traces exceed the {max_total_points} point limit; "
                                     "select fewer spots or file_ids")
                raw = _read_exact(stream, point_bytes, f"spot {spot_id} FileID {file_id} chromatogram")
                chromatogram = [[float(x), float(y)] for x, y in struct.iter_unpack("<ff", raw)]
                intensities = [p[1] for p in chromatogram]
                total_points += num_points
                found.add(file_id)
                samples.append({
                    "file_id": file_id,
                    "peak_left": float(left),
                    "peak_top": float(top),
                    "peak_right": float(right),
                    "num_points": num_points,
                    "mean_intensity": sum(intensities) / len(intensities) if intensities else 0.0,
                    "max_intensity": max(intensities, default=0.0),
                    "chromatogram": chromatogram,
                })

            if strict and requested_set is not None:
                missing = [i for i in file_ids if i not in found]
                if missing:
                    raise ValueError(f"Requested FileID values were not found in spot_id {spot_id}: {missing}")

            spots.append({
                "spot_id": spot_id,
                "rt": float(rt),
                "ri": float(ri),
                "mz": float(mass),
                "drift": float(drift),
                "main_type": int(main_type),
                "num_samples": num_samples,
                "samples": samples,
            })
    return spots


def search_by_mz_range(spots: list[dict], min_mz: float, max_mz: float) -> list[dict]:
    return [s for s in spots if min_mz <= s.get("mz", 0.0) <= max_mz]


def search_by_rt_range(spots: list[dict], min_rt: float, max_rt: float) -> list[dict]:
    return [s for s in spots if min_rt <= s.get("rt", 0.0) <= max_rt]


def top_by_max_intensity(spots: list[dict], top_n: int = 20) -> list[dict]:
    """サンプル間で最大のクロマトグラム最大強度で降順に並べる（peak_top は RT 座標なので使わない）。"""
    ranked = [{
        "spot_id": s["spot_id"], "rt": s["rt"], "mz": s["mz"], "num_samples": s["num_samples"],
        "max_intensity": max((x.get("max_intensity", 0.0) for x in s["samples"]), default=0.0),
    } for s in spots]
    return sorted(ranked, key=lambda x: x["max_intensity"], reverse=True)[:top_n]


def summarize(spots: list[dict]) -> dict:
    if not spots:
        return {"total_spots": 0}
    samples = [x for s in spots for x in s["samples"]]
    return {
        "total_spots": len(spots),
        "rt_range": (min(s["rt"] for s in spots), max(s["rt"] for s in spots)),
        "mz_range": (min(s["mz"] for s in spots), max(s["mz"] for s in spots)),
        "total_traces": len(samples),
        "total_points": sum(x["num_peaks"] for x in samples),
        "unique_file_ids": sorted({x["file_id"] for x in samples}),
    }


# Preserve the original raw, uppercase-key API beside ``read_aef``.
from .raw.aef_parse import decode, decode_spots




