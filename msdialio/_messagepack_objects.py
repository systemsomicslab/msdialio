"""Strict MessagePack stream decoding shared by raw readers."""

import io

import msgpack


def _objects(data):
    unpacker = msgpack.Unpacker(
        io.BytesIO(data), raw=False, strict_map_key=False,
        unicode_errors="surrogateescape", max_buffer_size=max(len(data), 1),
    )
    objects = []
    end = 0
    for obj in unpacker:
        objects.append(obj)
        end = unpacker.tell()
    if end != len(data):
        raise ValueError("Incomplete MessagePack object")
    return objects

