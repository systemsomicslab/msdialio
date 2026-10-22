"""Strict MessagePack stream decoding shared by raw readers."""

import io

import msgpack


class Float32(float):
    """A float whose MessagePack representation must use the float32 tag.

    Raw readers create these automatically. Use Float32(new_value) when
    replacing a single-precision field; ordinary Python floats encode as float64.
    """


def _read_value(unpacker, data):
    position = unpacker.tell()
    if position >= len(data):
        raise ValueError("Incomplete MessagePack object")
    tag = data[position]
    if 0x90 <= tag <= 0x9f or tag in (0xdc, 0xdd):
        return [_read_value(unpacker, data) for _ in range(unpacker.read_array_header())]
    if 0x80 <= tag <= 0x8f or tag in (0xde, 0xdf):
        result = {}
        for _ in range(unpacker.read_map_header()):
            key = _read_value(unpacker, data)
            result[key] = _read_value(unpacker, data)
        return result
    value = unpacker.unpack()
    return Float32(value) if tag == 0xca else value


def _objects(data):
    unpacker = msgpack.Unpacker(
        io.BytesIO(data), raw=False, strict_map_key=False,
        unicode_errors="surrogateescape", max_buffer_size=max(len(data), 1),
    )
    objects = []
    try:
        while unpacker.tell() < len(data):
            objects.append(_read_value(unpacker, data))
    except msgpack.OutOfData as exc:
        raise ValueError("Incomplete MessagePack object") from exc
    return objects


def _one_object(data):
    objects = _objects(data)
    if len(objects) != 1:
        raise ValueError("Expected exactly one MessagePack object")
    return objects[0]

