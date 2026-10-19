"""Read individual-sample peak features from PAI2 files."""

from __future__ import annotations

from enum import Enum

import msgpack

from ._msgpack_lz4 import decompress_packed_msgpack
from ._sources import open_binary_source


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


def _peak_record(data: list) -> dict:
    """Convert one raw ChromatogramPeakFeature to a named peak record."""
    formula_data = data[26]
    formula = formula_data[0] if isinstance(formula_data, (list, tuple)) and formula_data else "Unknown"
    adduct_data = data[30]
    adduct = adduct_data[2] if isinstance(adduct_data, (list, tuple)) and len(adduct_data) > 2 else "Unknown"

    # MS/MS: Key18 MS2RawSpectrumID / Key19 MS2RawSpectrumID2CE / Key24 Spectrum.
    # The actual spectrum is usually stored in the corresponding .dcl file.
    ms2_id = data[18] if len(data) > 18 and isinstance(data[18], int) else -1
    ce_map = data[19] if len(data) > 19 and isinstance(data[19], dict) else {}
    spectrum = data[24] if len(data) > 24 and isinstance(data[24], list) else []
    msms_spectrum = [[p[0], p[1]] for p in spectrum if isinstance(p, (list, tuple)) and len(p) >= 2]
    center = _convert_to_times(data[4])
    left = _convert_to_times(data[3])
    right = _convert_to_times(data[5])

    return {
        "PeakID": data[11],
        "Name": data[25],
        "Mz": data[43],
        "RT": center.get("rt"),
        "RI": center.get("ri"),
        "DriftTime": center.get("dt"),
        "PeakHeightTop": data[7],
        "PeakHeightLeft": data[6],
        "PeakHeightRight": data[8],
        "PeakAreaAboveZero": data[9],
        "PeakAreaAboveBaseline": data[10],
        "RTLeft": left.get("rt"),
        "RILeft": left.get("ri"),
        "DriftTimeLeft": left.get("dt"),
        "RTRight": right.get("rt"),
        "RIRight": right.get("ri"),
        "DriftTimeRight": right.get("dt"),
        "SignalToNoise": data[40][1],
        "IonMode": IonMode.from_int(data[22]),
        "Formula": formula,
        "Ontology": data[27],
        "SMILES": data[28],
        "InChIKey": data[29],
        "AdductType": adduct,
        "CollisionCrossSection": data[31],
        "Comment": data[38],
        "HasMSMS": bool(ce_map) or ms2_id >= 0,
        "MS2RawSpectrumID": ms2_id,
        "CollisionEnergies": sorted({float(v) for v in ce_map.values()}) if ce_map else [],
        "MSMSPeakCount": len(msms_spectrum),
        "MSMSSpectrum": msms_spectrum,
    }


def _deserialize(data: bytes, slim: bool = True) -> list[dict]:
    unpacked = msgpack.unpackb(
        decompress_packed_msgpack(data), raw=False, strict_map_key=False
    )
    peaks = []
    for raw_peak in unpacked:
        peak = _peak_record(raw_peak)
        peaks.append({key: peak[key] for key in _SLIM_FIELDS} if slim else peak)
    return peaks


_SLIM_FIELDS = (
    "PeakID", "Name", "Mz", "RT", "RI", "DriftTime", "PeakHeightTop",
    "PeakAreaAboveZero", "IonMode", "HasMSMS", "MS2RawSpectrumID", "MSMSPeakCount",
)


def read_pai2(source, slim: bool = True) -> list[dict]:
    """Read PAI2 peak records.

    ``slim=True`` keeps peak identifiers, RT/RI/drift time, primary
    measurements, ion mode, and MS/MS availability. Set ``slim=False`` to
    include boundaries, annotations, and the embedded spectrum. Accepts a path
    or binary stream; caller-owned streams stay open.
    """
    with open_binary_source(source) as f:
        return _deserialize(f.read(), slim=slim)
