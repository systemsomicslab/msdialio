"""Read MS-DIAL 5 database storage ZIP archives without extracting files."""

import io
import zipfile

from ._metadata import decode, large_list, named
from ._sources import open_binary_source

_GROUPS = {"MetabolomicsDataBases": "MetabolomicsDB",
           "ProteomicsDataBases": "ProteomicsDB",
           "EadLipidomicsDatabases": "EadLipidomicsDB"}


def _entries(source):
    with open_binary_source(source) as stream, zipfile.ZipFile(stream) as archive:
        entries = {}
        for info in archive.infolist():
            if info.is_dir():
                continue
            name = info.filename.replace("\\", "/")
            if name in entries:
                raise ValueError(f"Duplicate DBS entry: {name}")
            entries[name] = archive.read(info)
    if "Storage" not in entries:
        raise ValueError("DBS archive has no Storage entry")
    return entries


def read_raw_dbs(source) -> dict:
    """Return decoded Storage and all original ZIP entry bytes in Entries.

    ZIP separators are normalized to '/'. Unknown payloads are retained as bytes.
    The source must be seekable; a caller-owned stream is left open.
    """
    entries = _entries(source)
    return {"Storage": decode(entries["Storage"]), "Entries": entries}


def read_dbs(source) -> list[dict]:
    """Return database metadata, annotator pairs and molecule references.

    References retain C# field names; Spectrum contains stored peak arrays with
    mass and intensity at positions 0 and 1. Proteomics and EAD payloads are
    returned as Data bytes with References=None, without interpreting them.
    Both current flat and legacy nested ZIP layouts are supported.
    """
    raw = read_raw_dbs(source)
    storage, entries = raw["Storage"], raw["Entries"]
    if not isinstance(storage, dict):
        raise ValueError("Expected a DBS Storage map")
    records = []
    for group, prefix in _GROUPS.items():
        for item in storage.get(group) or []:
            metadata = item["DataBase"]
            if group == "EadLipidomicsDatabases" and isinstance(metadata, dict):
                dbid = metadata.get("Id")
            elif isinstance(metadata, list):
                id_key = 2 if group == "ProteomicsDataBases" else 0
                dbid = metadata[id_key] if len(metadata) > id_key else None
            else:
                dbid = None
            if not isinstance(dbid, str):
                raise ValueError("Invalid DBS database metadata")
            path = f"{prefix}/{dbid}"
            if path + "/DataBase" in entries:
                payload = entries[path + "/DataBase"]
            elif path in entries:
                with zipfile.ZipFile(io.BytesIO(entries[path])) as nested:
                    payload = nested.read("DataBase")
            else:
                raise ValueError(f"Missing DBS database payload: {path}/DataBase")
            record = {"DataBaseID": dbid, "Kind": prefix,
                      "DataBase": named(metadata, "MoleculeDataBase") if group == "MetabolomicsDataBases" else metadata,
                      "Pairs": item.get("Pairs", []), "References": None}
            if group == "MetabolomicsDataBases":
                record["References"] = [named(row, "MoleculeMsReference") for row in large_list(payload)]
            else:
                record["Data"] = payload
            records.append(record)
    return records
