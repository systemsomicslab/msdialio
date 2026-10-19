"""MS-DIAL ``.dcl``（デコンボリューション済み MS/MS = MSDecResult）の単体パーサ。

``.dcl`` は他の出力（.arf/.arf2/.pai2）と違い MessagePack ではなく、C# の
BinaryWriter が書く独自バイナリ（リトルエンディアン）。個別測定の ``.dcl`` は同名
``.pai2`` のピーク順に、アラインメントの ``.dcl`` はスポット順に対応する。

バイト構造（MS-DIAL MsdecResultsReader.ReadMSDecResultVer1 準拠）:
  ヘッダ(11B): "DC"(2) + version int32 + annotation フラグ bool(1) + 結果数 int32
  seekpointer 配列: 結果数 × int64
  各 MSDecResult（seekpointer の位置から）:
    Scan(60B)   : SeekPoint i64, ScanID i32, RawSpectrumID i32, PrecursorMz f64,
                  IonMode i32, RT f64, RI f64, Drift f64, Mz f64
    Quant(40B)  : ModelPeakMz, ModelPeakHeight, ModelPeakArea, IntegratedHeight, IntegratedArea (f64×5)
    Scoring(20B): AmplitudeScore, ModelPeakPurity, ModelPeakQuality, SignalNoiseRatio, EstimatedNoise (f32×5)
    Counts(12B) : spectraNumber i32, datapointNumber i32, modelMassNumber i32
    Spectrum    : spectraNumber × (Mass f64, Intensity f64, PeakQuality i32) = 20B/peak
    （以降の Chromatogram / ModelMasses / Annotation ブロックは読まない）

使い方（Python）::

    from msdial_parsers import read_dcl
    results = read_dcl("sample.dcl")
    selected = read_dcl("aligned.dcl", dcl_indices=[12, 345])

使い方（CLI）::

    python -m msdial_parsers dcl <path.dcl> [--csv msdec.csv]

依存: 標準ライブラリのみ。切り出し元: metabolomix/dcl/reader.py。
"""

from __future__ import annotations

import io
import struct
from collections.abc import Iterable

from ._sources import open_binary_source

DCL_MAGIC = b"DC"


def _read_exact(stream, size: int) -> bytes:
    chunks = []
    remaining = size
    while remaining:
        chunk = stream.read(remaining)
        if not chunk:
            raise ValueError("Unexpected end of DCL file")
        chunks.append(chunk)
        remaining -= len(chunk)
    return b"".join(chunks)


def _read_records(stream, dcl_indices: Iterable[int] | None = None) -> list[dict]:
    base_offset = stream.tell()
    magic = _read_exact(stream, 2)
    if magic != DCL_MAGIC:
        raise ValueError("先頭が 'DC' でないため .dcl 形式ではありません。")
    header = magic + _read_exact(stream, 9)

    count = struct.unpack_from("<i", header, 7)[0]
    if count < 0:
        raise ValueError(f"Invalid DCL record count: {count}")
    table = _read_exact(stream, count * 8)
    seekpoints = struct.unpack(f"<{count}q", table) if count else ()

    if dcl_indices is None:
        selected = range(count)
    else:
        requested = list(dcl_indices)
        for index in requested:
            if not isinstance(index, int):
                raise TypeError("dcl_indices must contain integers")
            if index < 0 or index >= count:
                raise IndexError(f"dcl_index {index} is out of range for {count} records")
        selected = sorted(set(requested))

    results = []
    for idx in selected:
        offset = seekpoints[idx]
        if offset < 0:
            raise ValueError(f"Invalid DCL record offset for dcl_index {idx}: {offset}")
        stream.seek(base_offset + offset)
        scan = _read_exact(stream, 60)
        scan_id, raw_spec_id = struct.unpack_from("<2i", scan, 8)
        precursor_mz = struct.unpack_from("<d", scan, 16)[0]
        ion_mode = struct.unpack_from("<i", scan, 24)[0]
        rt = struct.unpack_from("<d", scan, 28)[0]
        quant = _read_exact(stream, 40)
        _model_mz, model_h, _model_a, _int_h, _int_a = struct.unpack("<5d", quant)
        scoring = _read_exact(stream, 20)
        _amp, _purity, _quality, sn, est_noise = struct.unpack("<5f", scoring)
        counts = _read_exact(stream, 12)
        spectra_n, _datapoint_n, _model_n = struct.unpack("<3i", counts)
        if spectra_n < 0:
            raise ValueError(f"Invalid MS/MS peak count at dcl_index {idx}: {spectra_n}")
        raw_spectrum = _read_exact(stream, spectra_n * 20)
        spectrum = [
            list(struct.unpack_from("<2d", raw_spectrum, peak_offset))
            for peak_offset in range(0, len(raw_spectrum), 20)
        ]

        results.append({
            "dcl_index": idx,          # .dcl 内の順序（= .pai2 のピーク順 / アラインメントのスポット順）
            "scan_id": scan_id,
            "raw_spec_id": raw_spec_id,
            "precursor_mz": precursor_mz,
            "ion_mode": ion_mode,      # 0 = Positive, 1 = Negative
            "rt": rt,
            "model_peak_height": model_h,
            "signal_to_noise": sn,
            "estimated_noise": est_noise,
            "n_msms_peaks": spectra_n,
            "msms_spectrum": spectrum,
        })
    return results


def _deserialize(raw: bytes) -> list[dict]:
    """Deserialize a DCL byte string."""
    return _read_records(io.BytesIO(raw))


def read_dcl(source, dcl_indices: Iterable[int] | None = None) -> list[dict]:
    """Read all DCL records, or only the selected record indices.

    ``dcl_indices`` refers to the record order stored in the DCL file. In aligned
    DCL files this corresponds to spot order; in individual DCL files it
    corresponds to PAI2 peak order. Selected records are returned in file order;
    duplicate indices are collapsed. The source must be a path or seekable
    binary stream. Caller-owned streams are left open.
    """
    with open_binary_source(source) as f:
        return _read_records(f, dcl_indices)



