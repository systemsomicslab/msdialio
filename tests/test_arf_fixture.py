from pathlib import Path

from msdialio.arf_parse import read_arf


DATA = Path(__file__).parent / "data"


def test_reads_peak_properties_fixture():
    spots = read_arf(DATA / "alignment_peaks.arf")

    assert len(spots) == 187
    first = spots[0]
    assert first["SpotIndex"] == 0
    assert "MasterAlignmentID" not in first
    assert len(first["SamplePeaks"]) == 1
    peak = first["SamplePeaks"][0]
    assert peak["FileID"] == 0
    assert peak["FileName"].startswith("Nega_Swath_")
    assert peak["Mz"] == 139.00817
    assert peak["PeakHeightTop"] == 2581.4375
    assert peak["PeakAreaAboveZero"] == 13872.236593749998
    assert peak["RT"] == 2.9854
    assert peak["RI"] is None
    assert peak["DriftTime"] is None

    full = read_arf(DATA / "alignment_peaks.arf", slim=False)
    assert "PeakShape" in full[0]["SamplePeaks"][0]
