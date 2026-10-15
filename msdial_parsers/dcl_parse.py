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

    from msdial_parsers.dcl_parse import read_dcl, find_msms
    results = read_dcl("sample.dcl", top_n_peaks=50)
    hits = find_msms(results, precursor_mz=760.585, rt=12.3)

使い方（CLI）::

    python -m msdial_parsers dcl <path.dcl> [--top-peaks 20] [--csv msdec.csv]

依存: 標準ライブラリのみ。切り出し元: metabolomix/dcl/reader.py。
"""

from __future__ import annotations

import os
import struct

DCL_MAGIC = b"DC"


def deserialize(raw: bytes, include_spectrum: bool = True,
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


def read_dcl(path, include_spectrum: bool = True, top_n_peaks: int | None = None) -> list[dict]:
    with open(path, "rb") as f:
        return deserialize(f.read(), include_spectrum=include_spectrum, top_n_peaks=top_n_peaks)


def find_msms(results: list[dict], precursor_mz: float, tol: float = 0.01,
              rt: float | None = None, rt_tol: float = 0.2) -> list[dict]:
    """precursor m/z（任意で RT）が一致し、MS/MS を持つ結果を返す。"""
    return [r for r in results
            if abs(r["precursor_mz"] - precursor_mz) <= tol
            and (rt is None or abs(r["rt"] - rt) <= rt_tol)
            and r["n_msms_peaks"] > 0]


def find_dcl_for_pai2(pai2_path: str) -> str | None:
    """``.pai2`` と同じ basename の ``.dcl`` を同じフォルダから探す。"""
    base = pai2_path[:-len(".pai2")] if pai2_path.endswith(".pai2") else pai2_path
    cand = base + ".dcl"
    return cand if os.path.exists(cand) else None


def attach_msms_to_features(features: list[dict], dcl_results: list[dict],
                            mz_tol: float = 0.01) -> int:
    """``.pai2`` のピーク（pai2_parse の辞書）に ``.dcl`` の MS/MS を索引対応で付ける。

    precursor m/z が ``mz_tol`` を超えて食い違う組は、対応が崩れているとみなして飛ばす。
    付けられた（MS/MS ピークが 1 本以上ある）件数を返す。
    """
    attached = 0
    for feat, res in zip(features, dcl_results):
        feat_mz = feat.get("m/z")
        if isinstance(feat_mz, (int, float)) and abs(feat_mz - res["precursor_mz"]) > mz_tol:
            continue
        feat["msms_spectrum"] = res["msms_spectrum"]
        feat["n_msms_peaks"] = res["n_msms_peaks"]
        feat["msms_peak_count"] = res["n_msms_peaks"]
        if res["n_msms_peaks"] > 0:
            attached += 1
    return attached


def summarize(results: list[dict]) -> dict:
    if not results:
        return {"total_results": 0}
    counts = sorted(r["n_msms_peaks"] for r in results if r["n_msms_peaks"] > 0)
    mzs = [r["precursor_mz"] for r in results if r["precursor_mz"] > 0]
    rts = [r["rt"] for r in results]
    return {
        "total_results": len(results),
        "with_msms": len(counts),
        "msms_rate_pct": round(len(counts) / len(results) * 100, 1),
        "msms_peak_count_median": counts[len(counts) // 2] if counts else 0,
        "msms_peak_count_max": counts[-1] if counts else 0,
        "precursor_mz_range": (round(min(mzs), 4), round(max(mzs), 4)) if mzs else None,
        "rt_range": (round(min(rts), 2), round(max(rts), 2)),
    }




