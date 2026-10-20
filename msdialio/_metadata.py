"""Shared decoding for dataset metadata and reference database archives."""

import io

import lz4.block
import msgpack

from ._metadata_schema import SCHEMAS
from ._messagepack_objects import _objects


def unpack(data):
    return msgpack.unpackb(data, raw=False, strict_map_key=False,
                           unicode_errors="surrogateescape")


def decompress(extension):
    if extension.code != 99:
        raise ValueError(f"Unsupported MessagePack extension: {extension.code}")
    stream = msgpack.Unpacker(io.BytesIO(extension.data), raw=False)
    try:
        size = next(stream)
    except StopIteration as exc:
        raise ValueError("Missing LZ4 uncompressed size") from exc
    if type(size) is not int or size < 0:
        raise ValueError("Invalid LZ4 uncompressed size")
    result = lz4.block.decompress(extension.data[stream.tell():], uncompressed_size=size)
    if len(result) != size:
        raise ValueError("LZ4 uncompressed size mismatch")
    return result


def decode(data):
    """Decode the v1 LZ4 extension or an uncompressed MessagePack object."""
    value = unpack(data)
    return unpack(decompress(value)) if isinstance(value, msgpack.ExtType) else value


def large_list(data):
    """Read LargeListMessagePack's concatenated, padded array blocks.

    C# reserves five bytes for each array header, even when the actual header
    takes only one or three bytes. Ordinary MessagePack array decoding would
    interpret the padding as records and drop real rows at the end.
    """
    if data == b"\xc0":
        return []
    rows = []
    for block in _objects(data):
        if not isinstance(block, msgpack.ExtType):
            raise ValueError("Expected a LargeListMessagePack LZ4 block")
        payload = decompress(block)
        stream = msgpack.Unpacker(io.BytesIO(payload), raw=False, strict_map_key=False)
        count = stream.read_array_header()
        if len(payload) < 5:
            raise ValueError("Incomplete LargeListMessagePack array header")
        records = _objects(payload[5:])
        if len(records) != count:
            raise ValueError("LargeListMessagePack record count mismatch")
        rows.extend(records)
    return rows


def named(value, schema):
    """Name known stored fields, retaining unknown numeric keys and values."""
    if value is None:
        return None
    if not isinstance(value, (list, dict)):
        raise ValueError(f"Expected an array or map for {schema}")
    fields = SCHEMAS[schema]
    result = {}
    for key, item in enumerate(value) if isinstance(value, list) else value.items():
        name, kind = fields.get(key, (key, None))
        if kind in SCHEMAS:
            item = named(item, kind)
        result[name] = item
    return result
