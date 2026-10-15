"""Compatibility checks for the original standalone decoder APIs."""

from pathlib import Path

from msdial_parsers.aef_parse import decode as decode_aef, decode_spots
from msdial_parsers.arf2_parse import decode as decode_arf2
from msdial_parsers.arf_parse import decode as decode_arf, iter_objects
from msdial_parsers.dcl_parse import decode as decode_dcl, decode_bytes
from msdial_parsers.pai2_parse import decode as decode_pai2


DATA = Path(__file__).parent / "data"


def test_original_arf_raw_decoder_api_is_available():
    path = DATA / "alignment_peaks.arf"
    objects = decode_arf(path)

    assert objects == list(iter_objects(path))
    assert len(objects) == 3
    assert len(objects[0]) == 187


def test_original_arf2_raw_decoder_api_is_available():
    objects = decode_arf2(DATA / "alignment_spots.arf2")

    assert len(objects) == 1
    assert len(objects[0]) == 5
    assert len(objects[0][3]) == 187


def test_original_pai2_raw_decoder_api_is_available():
    objects = decode_pai2(DATA / "sample_peaks.pai2")

    assert len(objects) == 1
    assert len(objects[0]) == 364


def test_original_dcl_full_layout_api_is_available():
    path = DATA / "sample_msdec.dcl"
    raw = path.read_bytes()
    decoded = decode_dcl(path)

    assert decoded == decode_bytes(raw)
    assert decoded["header"]["RecordCount"] == 364
    assert "pointers" in decoded and "records" in decoded


def test_original_aef_raw_api_is_available():
    path = DATA / "alignment_eic.aef"
    decoded = decode_aef(path, include_points=False)
    selected = decode_spots(path, [0], include_points=False)

    assert decoded["header"]["SpotCount"] == 187
    assert set(selected["spots"]) == {0}
