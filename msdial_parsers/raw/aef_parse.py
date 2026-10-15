"""Raw decoder for MS-DIAL ``.EIC.aef`` files (alignment EIC chromatograms, "CSS1").

The file stores, for every alignment spot and every sample, the extracted ion
chromatogram as a list of (x, intensity) points. It is a custom little-endian binary
format, not MessagePack. This module decodes every field and does not interpret them
further. Spots follow the alignment spot order of the ``.arf`` / ``.arf2`` files.

Binary layout
-------------
Header (14 bytes)::

    "CSS1" (4) | Reserved (6, zero in all files seen) | SpotCount int32

Pointer table: ``SpotCount`` x int64, the byte offset of each spot.

Each spot, starting at its pointer::

    RT f32 | RI f32 | Mass f32 | Drift f32 | MainType i8 | SampleCount i32      (21 B)
    then SampleCount x sample:
        FileID i32 | PointCount i32 | PeakTop f32 | PeakLeft f32 | PeakRight f32 (20 B)
        then PointCount x (X f32, Intensity f32)                                  (8 B each)

Note: PeakTop / PeakLeft / PeakRight are positions on the x axis (retention time),
not intensities.

Usage (Python)::

    from aef_parse import decode, decode_spots
    eic = decode("AlignmentResult_....EIC.aef")                    # every spot
    some = decode_spots("AlignmentResult_....EIC.aef", [0, 42])     # selected spots only

Usage (CLI)::

    python aef_parse.py <file.EIC.aef> [--json out.json] [--no-points] [--spot N ...]

Requires: Python standard library only.
"""

from __future__ import annotations

import argparse
import json
import os
import struct

MAGIC = b"CSS1"
HEADER_SIZE = 14

_SPOT = struct.Struct("<ffffbi")
_SAMPLE = struct.Struct("<iifff")
_POINT = struct.Struct("<ff")


def _read_exact(f, size: int, what: str) -> bytes:
    data = f.read(size)
    if len(data) != size:
        raise ValueError(f"unexpected end of file while reading {what}: "
                         f"expected {size} bytes, got {len(data)}")
    return data


def _read_header(f) -> tuple[dict, list[int]]:
    head = _read_exact(f, HEADER_SIZE, "header")
    if head[:4] != MAGIC:
        raise ValueError(f"not a CSS1 file (magic {head[:4]!r})")
    count = struct.unpack_from("<i", head, 10)[0]
    pointers = list(struct.unpack(f"<{count}q", _read_exact(f, 8 * count, "pointer table")))
    header = {"Magic": MAGIC.decode(), "Reserved": head[4:10].hex(), "SpotCount": count}
    return header, pointers


def _read_spot(f, pointer: int, include_points: bool) -> dict:
    f.seek(pointer)
    rt, ri, mass, drift, main_type, sample_count = _SPOT.unpack(_read_exact(f, _SPOT.size, "spot"))
    samples = []
    for _ in range(sample_count):
        file_id, point_count, top, left, right = _SAMPLE.unpack(_read_exact(f, _SAMPLE.size, "sample"))
        sample = {"FileID": file_id, "PointCount": point_count,
                  "PeakTop": top, "PeakLeft": left, "PeakRight": right}
        if include_points:
            raw = _read_exact(f, point_count * _POINT.size, "points")
            sample["Points"] = [list(p) for p in _POINT.iter_unpack(raw)]
        else:
            f.seek(point_count * _POINT.size, os.SEEK_CUR)
        samples.append(sample)
    return {"RT": rt, "RI": ri, "Mass": mass, "Drift": drift, "MainType": main_type,
            "SampleCount": sample_count, "Samples": samples}


def decode(path, include_points: bool = True) -> dict:
    """Decode the whole file. ``include_points=False`` skips the chromatogram points."""
    with open(path, "rb") as f:
        header, pointers = _read_header(f)
        spots = [_read_spot(f, p, include_points) for p in pointers]
    return {"header": header, "pointers": pointers, "spots": spots}


def decode_spots(path, spot_indices, include_points: bool = True) -> dict:
    """Decode only the spots at the given indices (direct access via the pointer table)."""
    with open(path, "rb") as f:
        header, pointers = _read_header(f)
        spots = {}
        for i in spot_indices:
            if not 0 <= i < len(pointers):
                raise IndexError(f"spot index {i} out of range (0..{len(pointers) - 1})")
            spots[i] = _read_spot(f, pointers[i], include_points)
    return {"header": header, "spots": spots}


def main(argv=None):
    ap = argparse.ArgumentParser(description="Decode an MS-DIAL .EIC.aef (CSS1) file.")
    ap.add_argument("path", help="path to the .EIC.aef file")
    ap.add_argument("--json", help="write the decoded content as JSON to this path")
    ap.add_argument("--no-points", action="store_true", help="skip chromatogram points")
    ap.add_argument("--spot", type=int, nargs="+", help="decode only these spot indices")
    args = ap.parse_args(argv)

    include_points = not args.no_points
    if args.spot:
        eic = decode_spots(args.path, args.spot, include_points)
        spots = list(eic["spots"].items())
    else:
        eic = decode(args.path, include_points)
        spots = list(enumerate(eic["spots"]))
    print(f"{args.path}: {eic['header']}")
    for i, s in spots[:5]:
        points = sum(x["PointCount"] for x in s["Samples"])
        print(f"  spot {i}: Mass={s['Mass']:.4f} RT={s['RT']:.3f} "
              f"samples={s['SampleCount']} points={points}")
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(eic, f)
        print(f"wrote JSON: {args.json}")


if __name__ == "__main__":
    main()
