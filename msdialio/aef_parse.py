"""MS-DIAL ``.EIC.aef`` (CSS1) reader.

AEF stores extracted ion chromatograms for each alignment spot and sample.
``peak_top`` is the RT coordinate of the peak apex. Stored peak bounds are
ignored because they are known to be incorrect in current files.

Examples::

    from msdialio import read_aef

    spots = read_aef("AlignmentResult.EIC.aef")
    selected = read_aef("AlignmentResult.EIC.aef", spot_ids=[12, 345], file_ids=[0, 1])

``spot_ids`` and ``file_ids`` restrict which chromatogram points are read from
the file. Omit them to read all spots and samples.
"""

from __future__ import annotations

import os
import struct
from collections.abc import Iterable

from ._sources import open_binary_source

MAGIC = b"CSS1"
_HEADER_SIZE = 14  # magic 4 + unused 6 + spot count 4


def _read_exact(stream, size: int, context: str) -> bytes:
    """Read exactly size bytes, including from streams that return short chunks."""
    chunks = []
    remaining = size
    while remaining:
        chunk = stream.read(remaining)
        if not chunk:
            received = size - remaining
            raise ValueError(
                f"Unexpected end of EIC/AEF while reading {context}: "
                f"expected {size} bytes, got {received}"
            )
        chunks.append(chunk)
        remaining -= len(chunk)
    return b"".join(chunks)


def read_aef(
    source,
    *,
    spot_ids: Iterable[int] | None = None,
    file_ids: Iterable[int] | None = None,
) -> list[dict]:
    """Read all AEF spots, or selected spots and sample chromatograms.

    ``spot_ids`` selects alignment spot indices and ``file_ids`` selects sample
    FileIDs within each selected spot. The reader seeks over unselected point
    data, so these options reduce the amount of chromatogram data read. Selected
    spots always include their metadata and return a ``samples`` list; samples
    include ``file_id``, ``peak_top`` (RT coordinate), and ``chromatogram``.

    Raises ``ValueError`` when a requested spot or FileID is missing. Accepts a
    path or seekable binary stream. A caller-provided stream remains open and is
    read from its current position.
    """
    requested_spots = None if spot_ids is None else sorted(set(spot_ids))
    requested_file_ids = None if file_ids is None else list(file_ids)
    requested_file_set = None if requested_file_ids is None else set(requested_file_ids)
    spots = []
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
        if num_spots < 0:
            raise ValueError(f"Invalid EIC spot count: {num_spots}")
        table_end = _HEADER_SIZE + 8 * num_spots

        if requested_spots is None:
            selected_spots = range(num_spots)
        else:
            selected_spots = requested_spots

        for spot_id in selected_spots:
            if not 0 <= spot_id < num_spots:
                raise ValueError(f"spot_id {spot_id} is out of range for {num_spots} EIC spots")

            stream.seek(base_offset + _HEADER_SIZE + 8 * spot_id)
            pointer = struct.unpack("<q", _read_exact(stream, 8, "spot pointer"))[0]
            if pointer < table_end or pointer >= file_size:
                raise ValueError(f"Invalid EIC spot pointer for spot_id {spot_id}: {pointer}")
            stream.seek(base_offset + pointer)
            rt, ri, mass, drift, main_type, num_samples = struct.unpack(
                "<ffffbi", _read_exact(stream, 21, f"spot {spot_id} header")
            )
            if num_samples < 0:
                raise ValueError(f"Invalid sample count for EIC spot_id {spot_id}: {num_samples}")

            samples = []
            found_file_ids = set()
            for sample_index in range(num_samples):
                file_id, num_points, top, _left, _right = struct.unpack(
                    "<iifff", _read_exact(stream, 20, f"spot {spot_id} sample {sample_index} header")
                )
                if num_points < 0:
                    raise ValueError(
                        f"Invalid point count for EIC spot_id {spot_id}, FileID {file_id}: {num_points}"
                    )
                point_bytes = 8 * num_points
                if requested_file_set is not None and file_id not in requested_file_set:
                    stream.seek(point_bytes, os.SEEK_CUR)
                    continue
                raw = _read_exact(stream, point_bytes, f"spot {spot_id} FileID {file_id} chromatogram")
                chromatogram = [[float(x), float(y)] for x, y in struct.iter_unpack("<ff", raw)]
                found_file_ids.add(file_id)
                samples.append({
                    "file_id": file_id,
                    "peak_top": float(top),
                    "chromatogram": chromatogram,
                })

            if requested_file_set is not None:
                missing = [file_id for file_id in requested_file_ids if file_id not in found_file_ids]
                if missing:
                    raise ValueError(
                        f"Requested FileID values were not found in spot_id {spot_id}: {missing}"
                    )

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
