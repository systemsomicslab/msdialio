from io import BytesIO
from pathlib import Path
import struct
import zipfile

import lz4.block
import msgpack
import pytest

from msdialio import raw, read_arf, read_arf2, read_pai2, read_mddata, read_dbs
from test_metadata import dataset, dbs_bytes


def float_widths(value):
    if isinstance(value, float):
        return 32 if isinstance(value, raw.Float32) else 64
    if isinstance(value, list):
        return [float_widths(item) for item in value]
    if isinstance(value, dict):
        return {key: float_widths(item) for key, item in value.items()}
    return None


@pytest.mark.parametrize("kind,filename", [
    ("mddata", "Dataset_2025_06_18_22_25_39.mddata"),
    ("dbs", "Dataset_2025_06_18_22_25_39_Loaded.msp2.dbs"),
])
def test_local_project_roundtrip(kind, filename):
    # Large user-provided originals stay local, outside the Git repository.
    source = Path(__file__).parent / "data" / "local-project" / filename
    if not source.exists():
        pytest.skip("Local project fixture is not installed")
    read = getattr(raw, f"read_raw_{kind}")
    original = read(source)
    stream = BytesIO()
    getattr(raw, f"write_raw_{kind}")(stream, original)
    stream.seek(0)
    saved = read(stream)
    if kind == "dbs":
        assert {k: v for k, v in saved["Entries"].items() if k != "Storage"} == {
            k: v for k, v in original["Entries"].items() if k != "Storage"}
        original, saved = original["Storage"], saved["Storage"]
    assert saved == original
    assert float_widths(saved) == float_widths(original)


@pytest.mark.parametrize("kind", ["pai2", "arf2", "arf", "mddata", "dbs"])
def test_roundtrip_preserves_mixed_float_widths(kind):
    # C# MessagePackBinary.ReadSingle rejects float64 (code 203). Python's
    # equality hides that incompatibility, so assert the uncompressed wire data.
    from msdialio._metadata import decompress
    from test_metadata import compressed, archive
    single = b"\xca" + struct.pack(">f", 1.25)
    double = b"\xcb" + struct.pack(">d", 1.2345678901234567)
    payload = b"\x82\xa6nested\x92" + single + double + b"\xa3key" + single
    original = compressed(payload)
    if kind == "dbs":
        original = archive({"Storage": original, "opaque": b"unchanged"})
    value = getattr(raw, f"read_raw_{kind}")(BytesIO(original))
    stream = BytesIO()
    getattr(raw, f"write_raw_{kind}")(stream, value)
    encoded = stream.getvalue()
    if kind == "dbs":
        with zipfile.ZipFile(BytesIO(encoded)) as z:
            assert z.read("opaque") == b"unchanged"
            encoded = z.read("Storage")
    assert decompress(msgpack.unpackb(encoded)) == payload


def test_float32_edits_and_special_values():
    import copy
    from msdialio._metadata import decompress
    # Retaining width must survive copying and explicitly typed field edits.
    payload = b"\x93\xca" + struct.pack(">f", -0.0)
    payload += b"\xca" + struct.pack(">f", float("inf"))
    payload += b"\xcb" + struct.pack(">d", 1.2345678901234567)
    value = raw.read_raw_mddata(BytesIO(payload))
    assert isinstance(value[0], float)
    assert isinstance(value[0], raw.Float32)
    assert type(value[2]) is float
    value = copy.deepcopy(value)
    value[1] = raw.Float32(2.5)
    stream = BytesIO()
    raw.write_raw_mddata(stream, value)
    expected = payload[:6] + b"\xca" + struct.pack(">f", 2.5) + payload[11:]
    assert decompress(msgpack.unpackb(stream.getvalue())) == expected


@pytest.mark.parametrize("kind,filename,reader", [
    ("pai2", "sample_peaks.pai2", read_pai2),
    ("arf2", "alignment_spots.arf2", read_arf2),
    ("arf", "alignment_peaks.arf", read_arf),
])
def test_real_files_roundtrip(kind, filename, reader):
    source = Path(__file__).parent / "data" / filename
    value = getattr(raw, f"read_raw_{kind}")(source)
    stream = BytesIO()
    getattr(raw, f"write_raw_{kind}")(stream, value)
    assert not stream.closed
    stream.seek(0)
    saved = getattr(raw, f"read_raw_{kind}")(stream)
    assert saved == value
    assert float_widths(saved) == float_widths(value)
    stream.seek(0)
    assert reader(stream, slim=False) == reader(source, slim=False)


@pytest.mark.parametrize("kind", ["pai2", "arf2", "arf", "mddata"])
def test_preserves_types_and_unknown_values(kind):
    value = [None, -1, {42: b"binary", "future": msgpack.ExtType(12, b"x")},
             "bad\udcff", 1.25, [], {}]
    if kind == "arf":
        value = [[value, {9: "metadata"}], [None]]
    stream = BytesIO()
    getattr(raw, f"write_raw_{kind}")(stream, value)
    stream.seek(0)
    assert getattr(raw, f"read_raw_{kind}")(stream) == value


def test_mddata_edit():
    value = dataset()
    value[0][0][1] = "edited"
    stream = BytesIO()
    raw.write_raw_mddata(stream, value)
    stream.seek(0)
    assert read_mddata(stream)["AnalysisFiles"][0]["AnalysisFileName"] == "edited"


@pytest.mark.parametrize("legacy", [False, True])
@pytest.mark.parametrize("windows", [False, True])
def test_dbs_preserves_payloads_and_updates_storage(legacy, windows):
    original = dbs_bytes([[0, 100.0, None, 0, [], "molecule"]], legacy=legacy, windows=windows)
    value = raw.read_raw_dbs(BytesIO(original))
    value["Storage"]["Future"] = {42: b"changed"}
    stream = BytesIO()
    raw.write_raw_dbs(stream, value)
    assert not stream.closed
    stream.seek(0)
    saved = raw.read_raw_dbs(stream)
    assert saved["Storage"] == value["Storage"]
    assert {k: v for k, v in saved["Entries"].items() if k != "Storage"} == {
        k: v for k, v in value["Entries"].items() if k != "Storage"}
    stream.seek(0)
    assert read_dbs(stream) == read_dbs(BytesIO(original))


def test_output_paths_atomic_and_overwrite(tmp_path, monkeypatch):
    path = tmp_path / "test.pai2"
    raw.write_raw_pai2(path, [1])
    with pytest.raises(FileExistsError):
        raw.write_raw_pai2(path, [2])
    assert raw.read_raw_pai2(path) == [1]
    raw.write_raw_pai2(path, [3], overwrite=True)
    assert raw.read_raw_pai2(path) == [3]

    def fail(*args, **kwargs):
        raise OSError("disk failure")
    monkeypatch.setattr(zipfile.ZipFile, "writestr", fail)
    with pytest.raises(OSError, match="disk failure"):
        raw.write_raw_dbs(path, {"Storage": {}, "Entries": {}}, overwrite=True)
    assert raw.read_raw_pai2(path) == [3]
    assert list(tmp_path.iterdir()) == [path]


def test_fixed_wire_header():
    stream = BytesIO()
    raw.write_raw_pai2(stream, [1, "value"])
    data = stream.getvalue()
    assert data[0] == 0xC9 and data[5:7] == b"\x63\xd2"
    assert struct.unpack(">I", data[1:5])[0] == len(data) - 6
    size = struct.unpack(">i", data[7:11])[0]
    assert msgpack.unpackb(lz4.block.decompress(data[11:], uncompressed_size=size)) == [1, "value"]


def test_arf_large_list_layout_matches_csharp_reader():
    # C# reads the array count then jumps to byte 5, ignoring header padding.
    path = Path(__file__).parent / "data" / "alignment_peaks.arf"
    stream = BytesIO()
    raw.write_raw_arf(stream, raw.read_raw_arf(path))
    def rows(data):
        ext = msgpack.unpackb(data)
        size = struct.unpack(">i", ext.data[1:5])[0]
        payload = lz4.block.decompress(ext.data[5:], uncompressed_size=size)
        header = msgpack.Unpacker(BytesIO(payload))
        count = header.read_array_header()
        records = msgpack.Unpacker(BytesIO(payload[5:]), raw=False, strict_map_key=False)
        result = [next(records) for _ in range(count)]
        assert records.tell() == len(payload) - 5
        return result
    assert rows(stream.getvalue()) == rows(path.read_bytes())


def test_invalid_input_does_not_touch_destination(tmp_path):
    path = tmp_path / "invalid.dbs"
    with pytest.raises(ValueError, match="duplicate"):
        raw.write_raw_dbs(path, {"Storage": {}, "Entries": {"a/b": b"a", "a\\b": b"b"}})
    assert not path.exists()
    with pytest.raises(TypeError):
        raw.write_raw_arf(path, [1])
    assert not path.exists()


def test_dbs_emits_windows_entry_names_in_zip_headers():
    # Python's ZipInfo constructor normalizes backslashes on Windows. Check
    # actual archive bytes: MS-DIAL uses Path.Combine + GetEntry with '\\'.
    stream = BytesIO()
    raw.write_raw_dbs(stream, {"Storage": {}, "Entries": {
        "MetabolomicsDB/lib/DataBase": b"payload"}})
    data = stream.getvalue()
    name = b"MetabolomicsDB\\lib\\DataBase"
    name_size = struct.unpack_from("<H", data, 26)[0]
    assert data[30:30 + name_size] == name
    central = data.index(b"PK\x01\x02")
    name_size = struct.unpack_from("<H", data, central + 28)[0]
    assert data[central + 46:central + 46 + name_size] == name
