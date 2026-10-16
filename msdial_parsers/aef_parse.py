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

``peak_top`` はピーク頂点の横軸（RT）座標で、強度ではない。保存された ``peak_left`` /
``peak_right`` は返り値に含めない。

使い方（Python）::

    from msdial_parsers.aef_parse import read_aef, read_aef_spots
    spots = read_aef("AlignmentResult_....EIC.aef")
    traces = read_aef_spots("AlignmentResult_....EIC.aef", [12, 345], file_ids=[0, 1])

使い方（CLI）::

    python -m msdial_parsers aef <path.EIC.aef> [--csv spots.csv] [--spot 12 --trace-csv eic.csv]

依存: 標準ライブラリのみ。切り出し元: metabolomix/eic/reader.py。
"""

from __future__ import annotations

import os
import struct

from ._sources import open_binary_source

MAGIC = b"CSS1"
_HEADER_SIZE = 14  # magic 4 + 未使用 6 + スポット数 4


def read_aef(source) -> list[dict]:
    """全スポットを順に読む。

    各サンプルは ``file_id`` / ``peak_top`` / ``chromatogram``（(RT, 強度) の点列）を持つ。
    全スポットの点列が必要でない場合は ``read_aef_spots`` で対象を絞る。
    path または seek 可能なバイナリストリームを受け取り、ストリームは閉じない。
    """
    results = []
    with open_binary_source(source) as f:
        base_offset = f.tell()
        magic = f.read(4)
        if magic != MAGIC:
            raise ValueError(f"Not a CSS1 file. Magic: {magic}")
        f.seek(base_offset + 10)
        num_spots = struct.unpack("<i", f.read(4))[0]
        pointers = struct.unpack(f"<{num_spots}q", f.read(8 * num_spots))

        for spot_id, ptr in enumerate(pointers):
            f.seek(base_offset + ptr)
            rt, ri, mass, drift, main_type, num_samples = struct.unpack("<ffffbi", f.read(21))
            samples = []
            for _ in range(num_samples):
                file_id, num_points, top, _left, _right = struct.unpack("<iifff", f.read(20))
                raw = f.read(8 * num_points)
                points = [[float(rt), float(intensity)]
                          for rt, intensity in struct.iter_unpack("<ff", raw)]
                samples.append({
                    "file_id": file_id,
                    "peak_top": float(top),
                    "chromatogram": points,
                })
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


def read_aef_spots(source, spot_ids, file_ids=None, *, max_total_points: int = 2_000_000,
                   strict: bool = False) -> list[dict]:
    """指定スポットだけをポインタ表から直接読み、点列（chromatogram）付きで返す。

    file_ids を渡すとそのサンプルだけを読み、他は点列を読み飛ばす。
    strict=True なら範囲外の spot_id と、見つからない file_id で例外にする。
    strict=False なら範囲外のスポットは結果から省き、file_id が無いスポットは
    ``samples`` を空にして返す（両者を区別できる）。
    path または seek 可能なバイナリストリームを受け取り、ストリームは閉じない。
    """
    requested_set = set(file_ids) if file_ids is not None else None
    spots = []
    total_points = 0
    with open_binary_source(source) as stream:
        base_offset = stream.tell()
        stream.seek(0, os.SEEK_END)
        file_size = stream.tell() - base_offset
        stream.seek(base_offset)
        magic = _read_exact(stream, 4, "CSS1 magic")
        if magic != MAGIC:
            raise ValueError(f"Not a CSS1 file. Magic: {magic}")
        stream.seek(base_offset + 10)
        num_spots = struct.unpack("<i", _read_exact(stream, 4, "spot count"))[0]
        table_end = _HEADER_SIZE + 8 * num_spots

        for spot_id in sorted(set(spot_ids)):
            if not 0 <= spot_id < num_spots:
                if strict:
                    raise ValueError(f"spot_id {spot_id} is out of range for {num_spots} EIC spots")
                continue
            stream.seek(base_offset + _HEADER_SIZE + 8 * spot_id)
            pointer = struct.unpack("<q", _read_exact(stream, 8, "spot pointer"))[0]
            if pointer < table_end or pointer >= file_size:
                raise ValueError(f"Invalid EIC spot pointer for spot_id {spot_id}: {pointer}")
            stream.seek(base_offset + pointer)
            rt, ri, mass, drift, main_type, num_samples = struct.unpack(
                "<ffffbi", _read_exact(stream, 21, f"spot {spot_id} header"))

            samples = []
            found = set()
            for sample_index in range(num_samples):
                file_id, num_points, top, _left, _right = struct.unpack(
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
                total_points += num_points
                found.add(file_id)
                samples.append({
                    "file_id": file_id,
                    "peak_top": float(top),
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


