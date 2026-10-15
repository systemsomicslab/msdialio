"""Raw decoder for MS-DIAL ``.dcl`` files (deconvoluted MS/MS, ``MSDecResult``).

Unlike ``.arf`` / ``.arf2`` / ``.pai2``, a ``.dcl`` file is not MessagePack. It is written
by C# ``BinaryWriter`` / ``BitConverter`` (little-endian). This module decodes every block
of every record into plain Python values and does not interpret them further. Records
of a per-measurement ``.dcl`` follow the peak order of the ``.pai2`` with the same base
name; records of an alignment ``.dcl`` follow the alignment spot order.

Binary layout (``MsdecResultsReader.ReadMSDecResultVer1`` in MS-DIAL)
---------------------------------------------------------------------
Header (11 bytes)::

    "DC" (2) | Version int32 | IsAnnotationInfoIncluded bool (1) | RecordCount int32

Seek pointer table: ``RecordCount`` x int64, the byte offset of each record.

Each record, starting at its seek pointer::

    Scan     (60 B): SeekPoint i64, ScanID i32, RawSpectrumID i32, PrecursorMz f64,
                     IonMode i32, RT f64, RI f64, Drift f64, Mz f64
    Quant    (40 B): ModelPeakMz, ModelPeakHeight, ModelPeakArea,
                     IntegratedHeight, IntegratedArea (f64 x 5)
    Score    (20 B): AmplitudeScore, ModelPeakPurity, ModelPeakQuality,
                     SignalNoiseRatio, EstimatedNoise (f32 x 5)
    Counts   (12 B): SpectrumCount i32, ChromatogramCount i32, ModelMassCount i32
    Spectrum       : SpectrumCount x (Mass f64, Intensity f64, PeakQuality i32)        20 B each
    Chromatogram   : ChromatogramCount x (ID i32, ChromX f64, Mass f64, Intensity f64)  28 B each
    ModelMasses    : ModelMassCount x f64                                               8 B each

The Chromatogram and ModelMasses layouts were confirmed against real files: record
lengths match exactly, ChromX runs around the record RT, and Mass equals ModelPeakMz.
When the header flag ``IsAnnotationInfoIncluded`` is set, an annotation block follows;
its layout is not decoded and its bytes are returned as ``TrailingBytes``. Any other
bytes left before the next record are returned the same way.

Usage (Python)::

    from dcl_parse import decode
    dcl = decode("sample.dcl")
    dcl["header"], dcl["records"][0]["Spectrum"]

Usage (CLI)::

    python dcl_parse.py <file.dcl> [--json out.json] [--no-arrays]

Requires: Python standard library only.
"""

from __future__ import annotations

import argparse
import base64
import json
import struct

MAGIC = b"DC"
HEADER_SIZE = 11

_SCAN = struct.Struct("<qiididddd")
_QUANT = struct.Struct("<5d")
_SCORE = struct.Struct("<5f")
_COUNTS = struct.Struct("<3i")
_SPECTRUM_PEAK = struct.Struct("<ddi")
_CHROM_PEAK = struct.Struct("<iddd")
_MODEL_MASS = struct.Struct("<d")


def decode_bytes(raw: bytes, include_arrays: bool = True) -> dict:
    """Decode the content of a ``.dcl`` file.

    ``include_arrays=False`` skips the Spectrum / Chromatogram / ModelMasses arrays
    (their counts are still returned), which is much lighter for large files.
    """
    if raw[:2] != MAGIC:
        raise ValueError(f"not a .dcl file (magic {raw[:2]!r}, expected {MAGIC!r})")
    version, has_annotation, count = struct.unpack_from("<i?i", raw, 2)
    pointers = list(struct.unpack_from(f"<{count}q", raw, HEADER_SIZE))

    # A record ends where the next record (in file order) starts.
    ordered = sorted(set(pointers)) + [len(raw)]
    record_end = {start: ordered[i + 1] for i, start in enumerate(ordered[:-1])}

    records = []
    for pointer in pointers:
        o = pointer
        (seek_point, scan_id, raw_spectrum_id, precursor_mz, ion_mode,
         rt, ri, drift, mz) = _SCAN.unpack_from(raw, o)
        o += _SCAN.size
        model_peak_mz, model_peak_height, model_peak_area, integrated_height, integrated_area = \
            _QUANT.unpack_from(raw, o)
        o += _QUANT.size
        amplitude_score, model_peak_purity, model_peak_quality, signal_noise_ratio, estimated_noise = \
            _SCORE.unpack_from(raw, o)
        o += _SCORE.size
        spectrum_count, chromatogram_count, model_mass_count = _COUNTS.unpack_from(raw, o)
        o += _COUNTS.size

        record = {
            "SeekPoint": seek_point,
            "ScanID": scan_id,
            "RawSpectrumID": raw_spectrum_id,
            "PrecursorMz": precursor_mz,
            "IonMode": ion_mode,
            "RT": rt,
            "RI": ri,
            "Drift": drift,
            "Mz": mz,
            "ModelPeakMz": model_peak_mz,
            "ModelPeakHeight": model_peak_height,
            "ModelPeakArea": model_peak_area,
            "IntegratedHeight": integrated_height,
            "IntegratedArea": integrated_area,
            "AmplitudeScore": amplitude_score,
            "ModelPeakPurity": model_peak_purity,
            "ModelPeakQuality": model_peak_quality,
            "SignalNoiseRatio": signal_noise_ratio,
            "EstimatedNoise": estimated_noise,
            "SpectrumCount": spectrum_count,
            "ChromatogramCount": chromatogram_count,
            "ModelMassCount": model_mass_count,
        }

        arrays = (
            ("Spectrum", _SPECTRUM_PEAK, spectrum_count),
            ("Chromatogram", _CHROM_PEAK, chromatogram_count),
            ("ModelMasses", _MODEL_MASS, model_mass_count),
        )
        for name, layout, n in arrays:
            if include_arrays:
                values = [list(layout.unpack_from(raw, o + i * layout.size)) for i in range(n)]
                record[name] = [v[0] for v in values] if layout is _MODEL_MASS else values
            o += n * layout.size

        end = record_end[pointer]
        record["TrailingBytes"] = raw[o:end] if o < end else b""
        records.append(record)

    return {
        "header": {"Magic": MAGIC.decode(), "Version": version,
                   "IsAnnotationInfoIncluded": has_annotation, "RecordCount": count},
        "pointers": pointers,
        "records": records,
    }


def decode(path, include_arrays: bool = True) -> dict:
    with open(path, "rb") as f:
        return decode_bytes(f.read(), include_arrays=include_arrays)


def _json_default(obj):
    if isinstance(obj, (bytes, bytearray)):
        return {"$bytes": base64.b64encode(obj).decode("ascii")}
    raise TypeError(f"cannot serialize {type(obj).__name__}")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Decode an MS-DIAL .dcl file (MSDecResult records).")
    ap.add_argument("path", help="path to the .dcl file")
    ap.add_argument("--json", help="write the decoded content as JSON to this path")
    ap.add_argument("--no-arrays", action="store_true",
                    help="skip Spectrum / Chromatogram / ModelMasses arrays")
    ap.add_argument("--show", type=int, default=5, help="number of records to print")
    args = ap.parse_args(argv)

    dcl = decode(args.path, include_arrays=not args.no_arrays)
    print(f"{args.path}: {dcl['header']}")
    for r in dcl["records"][:args.show]:
        print(f"  PrecursorMz={r['PrecursorMz']:.4f} RT={r['RT']:.3f} "
              f"Spectrum={r['SpectrumCount']} Chromatogram={r['ChromatogramCount']} "
              f"ModelMasses={r['ModelMassCount']} TrailingBytes={len(r['TrailingBytes'])}")
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(dcl, f, default=_json_default)
        print(f"wrote JSON: {args.json}")


if __name__ == "__main__":
    main()
