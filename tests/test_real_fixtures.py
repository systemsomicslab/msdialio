from pathlib import Path

import pytest

from msdialio import read_aef, read_arf2, read_dcl, read_pai2


DATA = Path(__file__).parent / "data"


def test_reads_arf2_fixture():
    spots = read_arf2(DATA / "alignment_spots.arf2")

    assert len(spots) == 187
    assert spots[0]["MasterAlignmentID"] == 0
    assert spots[0]["MassCenter"] == pytest.approx(139.00817)
    assert spots[0]["IonMode"] == "Negative"


def test_reads_pai2_fixture():
    peaks = read_pai2(DATA / "sample_peaks.pai2")

    assert len(peaks) == 364
    assert peaks[0]["Mz"] == pytest.approx(157.86403)
    assert peaks[0]["IonMode"] == "Negative"
    assert "PeakHeightTop" in peaks[0] and "PeakAreaAboveZero" in peaks[0]
    assert "MS2RawSpectrumID" in peaks[0]
    assert "RT" in peaks[0] and "RI" in peaks[0] and "DriftTime" in peaks[0]
    assert "time" not in peaks[0]


@pytest.mark.parametrize(
    ("filename", "expected_count", "expected_mz"),
    [
        ("alignment_msdec.dcl", 187, 139.00817),
        ("sample_msdec.dcl", 364, 157.86403),
    ],
)
def test_reads_dcl_fixtures(filename, expected_count, expected_mz):
    results = read_dcl(DATA / filename)

    assert len(results) == expected_count
    assert results[0]["precursor_mz"] == pytest.approx(expected_mz)
    assert results[0]["n_msms_peaks"] > 0


def test_reads_aef_fixture_and_selected_trace():
    path = DATA / "alignment_eic.aef"
    spots = read_aef(path)
    traces = read_aef(path, spot_ids=[0], file_ids=[0])

    assert len(spots) == 187
    assert sum(len(sample["chromatogram"]) for spot in spots for sample in spot["samples"]) == 14684
    assert traces[0]["spot_id"] == 0
    assert len(traces[0]["samples"][0]["chromatogram"]) == 41
