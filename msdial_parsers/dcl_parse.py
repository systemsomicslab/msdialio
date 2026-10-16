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
    results = read_dcl("sample.dcl", top_n_peaks=50)

使い方（CLI）::

    python -m msdial_parsers dcl <path.dcl> [--top-peaks 20] [--csv msdec.csv]

依存: 標準ライブラリのみ。切り出し元: metabolomix/dcl/reader.py。
"""

from __future__ import annotations

import struct

from ._sources import open_binary_source

DCL_MAGIC = b"DC"


def _deserialize(raw: bytes, include_spectrum: bool = True,
                 top_n_peaks: int | None = None) -> list[dict]:
    """``.dcl`` のバイト列を MSDecResult 辞書のリストにする。

    include_spectrum: False なら MS/MS ピーク列を読まない（メタのみ・軽量）。
    top_n_peaks: 各スペクトルを強度上位 N 本に間引く（m/z 昇順で返す）。None で全件。
    """
    if raw[:2] != DCL_MAGIC:
        raise ValueError("先頭が 'DC' でないため .dcl 形式ではありません。")

    count = struct.unpack_from("<i", raw, 7)[0]
    seekpoints = struct.unpack_from(f"<{count}q", raw, 11)

    results = []
    for idx, o in enumerate(seekpoints):
        scan_id, raw_spec_id = struct.unpack_from("<2i", raw, o + 8)
        precursor_mz = struct.unpack_from("<d", raw, o + 16)[0]
        ion_mode = struct.unpack_from("<i", raw, o + 24)[0]
        rt = struct.unpack_from("<d", raw, o + 28)[0]
        o += 60
        _model_mz, model_h, _model_a, _int_h, _int_a = struct.unpack_from("<5d", raw, o)
        o += 40
        _amp, _purity, _quality, sn, est_noise = struct.unpack_from("<5f", raw, o)
        o += 20
        spectra_n, _datapoint_n, _model_n = struct.unpack_from("<3i", raw, o)
        o += 12

        spectrum = []
        if include_spectrum and spectra_n > 0:
            for _ in range(spectra_n):
                mass, inten = struct.unpack_from("<2d", raw, o)
                o += 20  # mass f64 + intensity f64 + peakquality i32
                spectrum.append([mass, inten])
            if top_n_peaks is not None and len(spectrum) > top_n_peaks:
                spectrum = sorted(spectrum, key=lambda p: p[1], reverse=True)[:top_n_peaks]
                spectrum.sort(key=lambda p: p[0])

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


def read_dcl(source, include_spectrum: bool = True, top_n_peaks: int | None = None) -> list[dict]:
    """Read a DCL file from a path or binary stream; caller-owned streams stay open."""
    with open_binary_source(source) as f:
        return _deserialize(f.read(), include_spectrum=include_spectrum, top_n_peaks=top_n_peaks)



