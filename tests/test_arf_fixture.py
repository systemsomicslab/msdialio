from pathlib import Path

from msdial_parsers.arf_parse import read_arf


DATA = Path(__file__).parent / "data"


def test_reads_peak_properties_fixture():
    spots = read_arf(DATA / "alignment_peaks.arf")

    assert len(spots) == 187
    first = spots[0]
    assert first["spot_index"] == 0
    assert "MasterAlignmentID" not in first
    assert len(first["sample_peaks"]) == 1
    peak = first["sample_peaks"][0]
    assert peak["file_id"] == 0
    assert peak["file_name"].startswith("Nega_Swath_")
    assert peak["m/z"] == 139.00817
    assert peak["peak_height_top"] == 2581.4375
    assert peak["peak_area_above_zero"] == 13872.236593749998
    assert peak["ri"] is None
    assert peak["drift_time"] is None

    full = read_arf(DATA / "alignment_peaks.arf", slim=False)
    assert "peak_shape" in full[0]["sample_peaks"][0]
