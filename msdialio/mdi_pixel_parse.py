"""Selective readers for MS-DIAL imaging intensity files (``.mdi_pixel``)."""

from __future__ import annotations

from array import array
from collections.abc import Iterable
import operator
import os
import struct
import sys

from ._sources import open_binary_source

_HEADER = struct.Struct("<ddi")
_CHUNK_BYTES = 1024 * 1024


def _read_exact(stream, size):
    result = bytearray()
    while len(result) < size:
        chunk = stream.read(size - len(result))
        if not chunk:
            raise ValueError("Unexpected end of mdi_pixel file")
        result.extend(chunk)
    return result


def _layout(stream):
    base = stream.tell()
    stream.seek(0, os.SEEK_END)
    size = stream.tell() - base
    stream.seek(base)
    count = struct.unpack("<i", _read_exact(stream, 4))[0]
    if count < 0:
        raise ValueError("Invalid mdi_pixel image count")
    pixels = 0
    if count:
        _, _, pixels = _HEADER.unpack(_read_exact(stream, _HEADER.size))
        if pixels < 0:
            raise ValueError("Invalid mdi_pixel pixel count")
    stride = _HEADER.size + pixels * 8
    if size != 4 + count * stride:
        raise ValueError("mdi_pixel file size does not match image and pixel counts")
    return base, count, pixels, stride


def _metadata(stream, base, index, pixels, stride):
    stream.seek(base + 4 + index * stride)
    mz, drift, actual_pixels = _HEADER.unpack(_read_exact(stream, _HEADER.size))
    if actual_pixels != pixels:
        raise ValueError(f"Inconsistent mdi_pixel pixel count at image {index}")
    return {"image_index": index, "mz": mz, "drift": drift, "pixel_count": pixels}


def list_mdi_pixel(source) -> list[dict]:
    """List image metadata without reading intensity arrays.

    Accept a path or seekable binary stream, starting at its current position.
    Caller-owned streams remain open. All images must have the same pixel count,
    as required by the MS-DIAL reader. Indices are zero-based file positions,
    not peak IDs from another file. The format has no magic or version marker.
    """
    with open_binary_source(source) as stream:
        base, count, pixels, stride = _layout(stream)
        return [_metadata(stream, base, i, pixels, stride) for i in range(count)]


def read_mdi_pixel(source, *, image_indices: Iterable[int]) -> list[dict]:
    """Read only explicitly selected images, in requested order.

    ``image_indices`` is required to avoid accidental loading of an entire file.
    Duplicate indices are returned once, at their first requested position.
    Each result contains ``image_index``, ``mz``, ``drift``, ``pixel_count``, and
    ``intensities`` (an ``array('d')`` of float64 values in saved pixel order).
    No XY coordinates or image dimensions are stored in this file; callers must
    supply matching coordinates separately to reconstruct a spatial image.

    Only selected intensity arrays are read. Memory scales with their total
    pixel count, plus a bounded temporary buffer. Accept a path or seekable
    binary stream; caller-owned streams remain open and start at their current
    position. Invalid indices raise ValueError, noninteger indices TypeError.
    """
    indices = []
    seen = set()
    for value in image_indices:
        if isinstance(value, bool):
            raise TypeError("image_indices must contain integers, not bool")
        index = operator.index(value)
        if index not in seen:
            indices.append(index)
            seen.add(index)

    with open_binary_source(source) as stream:
        base, count, pixels, stride = _layout(stream)
        for index in indices:
            if not 0 <= index < count:
                raise ValueError(f"image_index {index} is out of range for {count} images")
        images = []
        for index in indices:
            image = _metadata(stream, base, index, pixels, stride)
            intensities = array("d")
            remaining = pixels * 8
            while remaining:
                size = min(remaining, _CHUNK_BYTES)
                intensities.frombytes(_read_exact(stream, size))
                remaining -= size
            if sys.byteorder != "little":
                intensities.byteswap()
            image["intensities"] = intensities
            images.append(image)
        return images
