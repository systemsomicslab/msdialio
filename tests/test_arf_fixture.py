from pathlib import Path

from msdial_parsers.arf_parse import read_arf


DATA = Path(__file__).parent / "data"


def test_reads_peak_properties_fixture():
    spots = read_arf(DATA / "alignment_peaks.arf")

    assert len(spots) == 187
    assert spots[0]["MasterAlignmentID"] == 0
    assert "m/z" in spots[0]
    assert "RT" in spots[0]
    assert "RI" in spots[0]
    assert "DriftTime" in spots[0]
    assert "MassCenter" not in spots[0]
    assert "HeightAverage" not in spots[0]
