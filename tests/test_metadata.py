"""MS-DIAL v1 LZ4 and LargeListMessagePack wire-format regression tests."""

from io import BytesIO
import json
import struct
import warnings
import zipfile

import lz4.block
import msgpack
import pytest

from msdialio import read_dbs, read_mddata
from msdialio.cli import main
from msdialio.raw import read_raw_dbs, read_raw_mddata


def compressed(payload):
    body = b"\xd2" + struct.pack(">i", len(payload))
    body += lz4.block.compress(payload, store_size=False)
    return msgpack.packb(msgpack.ExtType(99, body))


def packed(value):
    return compressed(msgpack.packb(value, use_bin_type=True))


def large_block(rows):
    count = len(rows)
    header = msgpack.packb(rows)[:1] if count < 16 else b"\xdc" + struct.pack(">H", count)
    payload = header.ljust(5, b"\0") + b"".join(msgpack.packb(row) for row in rows)
    return compressed(payload)


def archive(entries):
    stream = BytesIO()
    with zipfile.ZipFile(stream, "w") as z:
        for name, payload in entries.items():
            z.writestr(name, payload)
    return stream.getvalue()


def dataset():
    project = [None, None, "5.5", "/original", "project.mddata"]
    analysis = ["raw.mzML", "sample", 0, "Control", 1, 42, True,
                "sample.dcl", [], "sample.pai2"]
    return [[analysis], [[7, "aligned.arf2", "aligned"]], [], [], [], None,
            [project], None, None, {42: b"future"}]


def test_mddata_names_and_retains_unknown_values():
    value = dataset()
    stream = BytesIO(packed(value))
    result = read_mddata(stream)
    assert not stream.closed
    assert result["AnalysisFiles"][0]["AnalysisFileId"] == 42
    assert result["AlignmentFiles"][0]["FilePath"] == "aligned.arf2"
    assert result["Parameter"]["ProjectParam"]["MsdialVersionNumber"] == "5.5"
    assert result[9] == {42: b"future"}
    assert read_raw_mddata(BytesIO(packed(value))) == value
    assert read_mddata(BytesIO(msgpack.packb(value))) == result


@pytest.mark.parametrize("value", [[], None, {}, "not a dataset"])
def test_mddata_rejects_wrong_shape(value):
    with pytest.raises(ValueError):
        read_mddata(BytesIO(packed(value)))


def dbs_bytes(rows, *, legacy=False, windows=False):
    storage = {"MetabolomicsDataBases": [
        {"DataBase": ["lib", 4, 1, "original.msp"], "Pairs": [[0, {"Future": -1}]]}
    ], "ProteomicsDataBases": [], "EadLipidomicsDatabases": []}
    payload = large_block(rows[:1]) + large_block(rows[1:])
    key = "MetabolomicsDB/lib"
    if legacy:
        payload = archive({"DataBase": payload})
    else:
        key += "/DataBase"
    if windows:
        key = key.replace("/", "\\")
    return archive({"Storage": packed(storage), key: payload, "future": b"opaque"})


@pytest.mark.parametrize("legacy", [False, True])
@pytest.mark.parametrize("windows", [False, True])
def test_dbs_loads_all_blocks_and_both_zip_layouts(legacy, windows):
    rows = [[i, 100.5, None, 1, [[50.1, 123.0, "fragment"]], f"compound{i}"]
            for i in range(20)]
    data = dbs_bytes(rows, legacy=legacy, windows=windows)
    stream = BytesIO(data)
    result = read_dbs(stream)
    assert not stream.closed
    db = result[0]
    assert db["DataBaseID"] == "lib"
    assert db["DataBase"]["DataBaseSourceFilePath"] == "original.msp"
    assert db["Pairs"] == [[0, {"Future": -1}]]
    assert len(db["References"]) == 20
    assert db["References"][-1]["ScanID"] == 19
    assert db["References"][0]["Spectrum"] == [[50.1, 123.0, "fragment"]]
    raw = read_raw_dbs(BytesIO(data))
    assert raw["Entries"]["future"] == b"opaque"
    assert raw["Storage"]["MetabolomicsDataBases"][0]["DataBase"][0] == "lib"


def test_dbs_rejects_missing_payload_and_storage():
    with pytest.raises(ValueError, match="Storage"):
        read_dbs(BytesIO(archive({"other": b""})))
    raw = read_raw_dbs(BytesIO(dbs_bytes([[0]])))
    with pytest.raises(ValueError, match="Missing DBS"):
        read_dbs(BytesIO(archive({"Storage": raw["Entries"]["Storage"]})))


def test_dbs_empty_and_unsupported_database_payloads():
    assert read_dbs(BytesIO(archive({"Storage": packed({})}))) == []
    for group, prefix, metadata in (
        ("ProteomicsDataBases", "ProteomicsDB", [None, "library.fasta", "lib", 4]),
        ("EadLipidomicsDatabases", "EadLipidomicsDB", {"Id": "lib", "DatabaseFormat": 2, "Source": 8}),
    ):
        raw = read_raw_dbs(BytesIO(dbs_bytes([[0, 100.0, None, 0, [], "molecule"]])))
        storage = raw["Storage"]
        storage[group] = [{"DataBase": metadata, "Pairs": []}]
        entries = raw["Entries"]
        entries["Storage"] = packed(storage)
        entries[prefix + "/lib/DataBase"] = b"opaque"
        result = read_dbs(BytesIO(archive(entries)))
        assert result[0]["References"][0]["Name"] == "molecule"
        assert result[1]["DataBaseID"] == "lib"
        assert result[1]["DataBase"] == metadata
        assert result[1]["References"] is None
        assert result[1]["Data"] == b"opaque"


def test_dbs_rejects_ambiguous_zip_names():
    # Windows' zipfile writer already normalizes the two spellings, whereas
    # Unix retains both. Either archive must be rejected by the reader.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        data = archive({"Storage": packed({}), "a/b": b"first", "a\\b": b"second"})
    with pytest.raises(ValueError, match="Duplicate"):
        read_raw_dbs(BytesIO(data))


def test_truncated_and_unknown_extensions():
    from msdialio._metadata import large_list

    with pytest.raises(ValueError):
        read_raw_mddata(BytesIO(packed(dataset())[:-1]))
    with pytest.raises(ValueError, match="Unsupported"):
        read_raw_mddata(BytesIO(msgpack.packb(msgpack.ExtType(98, b"abc"))))
    with pytest.raises(ValueError):
        large_list(large_block([[1]])[:-1])
    with pytest.raises(ValueError, match="count"):
        large_list(compressed(b"\x92\0\0\0\0\x91\x01"))


def test_cli_paths_json_and_csv(tmp_path, capsys):
    mddata = tmp_path / "sample.mddata"
    mddata.write_bytes(packed(dataset()))
    output = tmp_path / "output.json"
    assert main(["mddata", str(mddata), "--json", str(output)]) == 0
    assert json.loads(output.read_text(encoding="utf-8"))["AnalysisFiles"][0]["AnalysisFileId"] == 42
    assert "1 records" in capsys.readouterr().out
    csv = tmp_path / "output.csv"
    assert main(["mddata", str(mddata), "--csv", str(csv)]) == 0
    assert "AnalysisFiles" in csv.read_text(encoding="utf-8")
    dbs = tmp_path / "sample.dbs"
    dbs.write_bytes(dbs_bytes([[0, 100.5, None, 0, [], "compound"]]))
    assert main(["dbs", str(dbs), "--json", str(output)]) == 0
    assert json.loads(output.read_text(encoding="utf-8"))[0]["References"][0]["Name"] == "compound"
