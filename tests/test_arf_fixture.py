from pathlib import Path

from msdial_parsers.arf_parse import read_arf


DATA = Path(__file__).parent / "data"


def test_reads_peak_properties_fixture():
    spots = read_arf(DATA / "alignment_peaks.arf")

    assert len(spots) == 187
    assert spots[0]["MasterAlignmentID"] == 0
