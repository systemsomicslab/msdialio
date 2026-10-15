import struct

import pytest

from msdial_parsers.aef_parse import read_aef, read_aef_spots


def test_reads_empty_css1_file(tmp_path):
    path = tmp_path / "empty.EIC.aef"
    path.write_bytes(b"CSS1" + bytes(6) + struct.pack("<i", 0))

    assert read_aef(path) == []
    assert read_aef_spots(path, []) == []


def test_rejects_non_css1_file(tmp_path):
    path = tmp_path / "invalid.aef"
    path.write_bytes(b"NOPE" + bytes(10))

    with pytest.raises(ValueError, match="Not a CSS1 file"):
        read_aef(path)
