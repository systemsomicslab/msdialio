"""Shared decoding for dataset metadata and reference database archives."""

import io

import lz4.block
import msgpack

from ._metadata_schema import SCHEMAS


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
    """Read all LargeListMessagePack records from bytes."""
    return list(iter_large_list(io.BytesIO(data)))


def _read_exact(stream, size):
    chunks = []
    remaining = size
    while remaining:
        chunk = stream.read(remaining)
        if not chunk:
            raise ValueError("Incomplete LargeListMessagePack block")
        chunks.append(chunk)
        remaining -= len(chunk)
    return b"".join(chunks)


def _large_list_payloads(stream):
    """Consume one MessagePack extension at a time, without reading ahead."""
    while prefix := stream.read(1):
        tag = prefix[0]
        if tag == 0xC0:
            if stream.read(1):
                raise ValueError("Unexpected data after nil LargeListMessagePack")
            return
        if tag in (0xC7, 0xC8, 0xC9):
            length_size = {0xC7: 1, 0xC8: 2, 0xC9: 4}[tag]
            length = int.from_bytes(_read_exact(stream, length_size), "big")
        elif 0xD4 <= tag <= 0xD8:
            length = 1 << (tag - 0xD4)
        else:
            raise ValueError("Expected a LargeListMessagePack LZ4 block")
        code = int.from_bytes(_read_exact(stream, 1), "big", signed=True)
        if code != 99:
            raise ValueError(f"Unsupported MessagePack extension: {code}")
        yield decompress(msgpack.ExtType(code, _read_exact(stream, length)))


def iter_large_list(stream):
    """Yield records while retaining only one decompressed block.

    C# reserves five bytes for each array header, even when the actual header
    takes only one or three bytes. Ordinary MessagePack array decoding would
    interpret the padding as records and drop real rows at the end.
    """
    for payload in _large_list_payloads(stream):
        if len(payload) < 5:
            raise ValueError("Incomplete LargeListMessagePack array header")
        buffer = io.BytesIO(payload)
        header = msgpack.Unpacker(buffer, raw=False)
        count = header.read_array_header()
        buffer.seek(5)
        records = msgpack.Unpacker(buffer, raw=False, strict_map_key=False,
                                  unicode_errors="surrogateescape",
                                  max_buffer_size=max(len(payload), 1))
        for _ in range(count):
            try:
                row = next(records)
            except StopIteration as exc:
                raise ValueError("LargeListMessagePack record count mismatch") from exc
            yield row
        if records.tell() != len(payload) - 5:
            raise ValueError("LargeListMessagePack record count mismatch")
        del payload, buffer, header, records


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
