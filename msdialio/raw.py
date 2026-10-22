"""Decode MS-DIAL MessagePack payloads without interpreting their fields.

Arrays, maps, unknown fields and sentinel values are retained. MessagePack
strings become str and binary values remain bytes. Invalid UTF-8 is preserved
with surrogateescape. These functions accept paths or binary streams and leave
caller-owned streams open. Raw readers are deliberately not re-exported at the
package root. DCL and AEF use different binary formats, not MessagePack.
Float32 values use a float subclass to retain their wire width for writers.
"""

from ._msgpack_lz4 import decompress_packed_msgpack
from ._sources import open_binary_source
from ._messagepack_objects import Float32, _objects, _one_object
from .mddata_parse import read_raw_mddata
from .dbs_parse import read_raw_dbs
from ._raw_write import (write_raw_arf, write_raw_arf2, write_raw_pai2,
                         write_raw_mddata, write_raw_dbs)
from .arf_parse import _decompress_lz4_payload, _payload_from_top_level_object

__all__ = ["read_raw_arf", "read_raw_arf2", "read_raw_pai2", "read_raw_mddata", "read_raw_dbs"]
__all__ += ["write_raw_arf", "write_raw_arf2", "write_raw_pai2", "write_raw_mddata", "write_raw_dbs"]
__all__ += ["Float32"]


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
        return _one_object(decompress_packed_msgpack(stream.read()))
