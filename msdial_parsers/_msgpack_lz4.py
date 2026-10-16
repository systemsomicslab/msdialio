"""Internal helpers for MS-DIAL's LZ4-compressed MessagePack payloads."""

from __future__ import annotations

import io

import lz4.block
import msgpack


def decompress_packed_msgpack(data: bytes) -> bytes:
    """Return the decompressed payload from an MS-DIAL MessagePack envelope."""
    _header, compressed = msgpack.unpackb(data, raw=False)
    unpacker = msgpack.Unpacker(io.BytesIO(compressed), raw=False)
    size = next(unpacker)
    return lz4.block.decompress(compressed[unpacker.tell():], uncompressed_size=size)
