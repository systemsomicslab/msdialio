"""Decode MS-DIAL MessagePack payloads without interpreting their fields.

Arrays, maps, unknown fields and sentinel values are retained. MessagePack
strings become str and binary values remain bytes. Invalid UTF-8 is preserved
with surrogateescape. These functions accept paths or binary streams and leave
caller-owned streams open. Raw readers are deliberately not re-exported at the
package root. DCL and AEF use different binary formats, not MessagePack.
"""

import io

import msgpack

from ._msgpack_lz4 import decompress_packed_msgpack
from ._sources import open_binary_source
from .arf_parse import _decompress_lz4_payload, _payload_from_top_level_object

__all__ = ["read_raw_arf", "read_raw_arf2", "read_raw_pai2"]


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


def read_raw_arf(source) -> list[list]:
    """Return blocks, each containing its decoded MessagePack objects.

    Only compression envelopes are removed. Block order, object order and
    payload nesting are retained, including non-peak metadata.
    """
    with open_binary_source(source) as stream:
        envelopes = _objects(stream.read())
    blocks = []
    for envelope in envelopes:
        payload = _payload_from_top_level_object(envelope)
        if payload is None:
            raise ValueError("Expected an ARF LZ4 MessagePack envelope")
        blocks.append(_objects(_decompress_lz4_payload(payload)))
    return blocks


def read_raw_arf2(source) -> list:
    """Return decoded payload objects in order, preserving their nesting.

    The outer list represents the MessagePack stream; even a single payload
    object occupies one element. It is not a flattened list of spots.
    """
    with open_binary_source(source) as stream:
        return _objects(decompress_packed_msgpack(stream.read()))


def read_raw_pai2(source):
    """Return the single decoded payload object (normally a peak array)."""
    with open_binary_source(source) as stream:
        return msgpack.unpackb(
            decompress_packed_msgpack(stream.read()), raw=False,
            strict_map_key=False, unicode_errors="surrogateescape",
        )
