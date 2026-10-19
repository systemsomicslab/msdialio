"""Read stored alignment spot representatives from ARF2 files.

Use ARF for sample-level peaks. Validate file ordinals against stored alignment
identifiers before joining records across files.
"""

from __future__ import annotations

import io
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


def _unpack_lz4_packed_msgpack(data: bytes) -> list:
    """展開後の複数MessagePackオブジェクトをリストにする。"""
    decompressed = decompress_packed_msgpack(data)
    return list(msgpack.Unpacker(io.BytesIO(decompressed), raw=False, strict_map_key=False))


def _iter_raw_spots(datas: list):
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


def _spot_record(data: list) -> dict | None:
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


def _deserialize(data: bytes, slim: bool = True) -> list[dict]:
    """``.arf2`` のバイト列をスポット辞書のリストにする。"""
    raw_spots = _iter_raw_spots(_unpack_lz4_packed_msgpack(data))
    records = []
    for raw in raw_spots:
        record = _spot_record(raw)
        if record is None:
            continue
        records.append({key: record[key] for key in _SLIM_FIELDS} if slim else record)
    return records


_SLIM_FIELDS = (
    "MasterAlignmentID", "AlignmentID", "RT", "MassCenter", "IonMode", "Name", "HeightAverage",
)


def read_arf2(source, slim: bool = True) -> list[dict]:
    """Read ARF2 spot records.

    ``slim=True`` keeps identifiers, primary spot measurements, ion mode, and
    name. Set ``slim=False`` to include identification and quality metadata.
    Accepts a path or binary stream; caller-owned streams stay open.
    """
    with open_binary_source(source) as f:
        return _deserialize(f.read(), slim=slim)


