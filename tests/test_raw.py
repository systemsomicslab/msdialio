from io import BytesIO
from pathlib import Path

import lz4.block
import msgpack
import pytest

import msdialio
from msdialio.raw import read_raw_arf, read_raw_arf2, read_raw_pai2


def envelope(payload):
    compressed = msgpack.packb(len(payload)) + lz4.block.compress(payload, store_size=False)
    return msgpack.packb([0, compressed])


def test_raw_preserves_unknown_fields_nesting_and_sentinels():
    value = [[[1, -1, None, {42: b"binary", "FutureField": [123]}]]]
    packed = msgpack.packb(value, use_bin_type=True)
    for reader, expected in [(read_raw_pai2, value), (read_raw_arf2, [value])]:
        stream = BytesIO(envelope(packed))
        assert reader(stream) == expected
        assert not stream.closed
    assert not hasattr(msdialio, "read_raw_arf")


def test_arf_retains_blocks_objects_and_metadata():
    data = envelope(msgpack.packb([1, [2]]) + msgpack.packb({7: "metadata"}))
    data += envelope(msgpack.packb([3]))
    assert read_raw_arf(BytesIO(data)) == [[[1, [2]], {7: "metadata"}], [[3]]]


def test_raw_rejects_truncated_payload_and_unknown_envelope():
    with pytest.raises(ValueError):
        read_raw_arf2(BytesIO(envelope(b"\x92\x01")))
    with pytest.raises(ValueError):
        read_raw_arf(BytesIO(msgpack.packb([1, 2, 3])))


@pytest.mark.parametrize(("reader", "filename"), [
    (read_raw_arf, "alignment_peaks.arf"),
    (read_raw_arf2, "alignment_spots.arf2"),
    (read_raw_pai2, "sample_peaks.pai2"),
])
def test_raw_real_files(reader, filename):
    assert reader(Path(__file__).parent / "data" / filename)
