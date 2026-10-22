"""Read MS-DIAL 5 database storage ZIP archives without extracting files."""

from contextlib import contextmanager
import zipfile

from ._metadata import decode, iter_large_list, named
from ._sources import open_binary_source

_GROUPS = {"MetabolomicsDataBases": "MetabolomicsDB",
           "ProteomicsDataBases": "ProteomicsDB",
           "EadLipidomicsDatabases": "EadLipidomicsDB"}


def _entry_index(archive):
    entries = {}
    for info in archive.infolist():
        if info.is_dir():
            continue
        name = info.filename.replace("\\", "/")
        if name in entries:
            raise ValueError(f"Duplicate DBS entry: {name}")
        entries[name] = info
    return entries


def _storage(archive, entries):
    if "Storage" not in entries:
        raise ValueError("DBS archive has no Storage entry")
    storage = decode(archive.read(entries["Storage"]))
    if not isinstance(storage, dict):
        raise ValueError("Expected a DBS Storage map")
    return storage


def _entries(source):
    with open_binary_source(source) as stream, zipfile.ZipFile(stream) as archive:
        entries = _entry_index(archive)
        if "Storage" not in entries:
            raise ValueError("DBS archive has no Storage entry")
        return {name: archive.read(info) for name, info in entries.items()}


def read_raw_dbs(source) -> dict:
    """Return decoded Storage and all original ZIP entry bytes in Entries.

    ZIP separators are normalized to '/'. Unknown payloads are retained as bytes.
    The source must be seekable; a caller-owned stream is left open.
    """
    entries = _entries(source)
    return {"Storage": decode(entries["Storage"], preserve_float_width=True), "Entries": entries}


def _databases(storage, database_ids):
    if isinstance(database_ids, str):
        database_ids = {database_ids}
    elif database_ids is not None:
        database_ids = set(database_ids)
    databases = []
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
            if database_ids is None or dbid in database_ids:
                databases.append({
                    "DataBaseID": dbid, "Kind": prefix,
                    "DataBase": named(metadata, "MoleculeDataBase") if prefix == "MetabolomicsDB" else metadata,
                    "Pairs": item.get("Pairs", []),
                })
    if database_ids is not None:
        missing = database_ids - {db["DataBaseID"] for db in databases}
        if missing:
            raise KeyError(f"Unknown DBS database IDs: {sorted(missing)}")
    return databases


def list_dbs(source) -> list[dict]:
    """List database metadata without opening any database payloads.

    Returns DataBaseID, Kind, DataBase and Pairs using read_dbs's metadata
    schema and ordering. References, Data and compound counts are not included.
    Only the ZIP directory and Storage entry are read. Payload presence and
    integrity are not checked. The source must be seekable; caller-owned
    streams remain open.
    """
    with open_binary_source(source) as stream, zipfile.ZipFile(stream) as archive:
        entries = _entry_index(archive)
        return _databases(_storage(archive, entries), None)


@contextmanager
def _database_stream(archive, entries, database):
    path = f"{database['Kind']}/{database['DataBaseID']}"
    if path + "/DataBase" in entries:
        with archive.open(entries[path + "/DataBase"]) as stream:
            yield stream
    elif path in entries:
        with archive.open(entries[path]) as stream, zipfile.ZipFile(stream) as nested:
            nested_entries = _entry_index(nested)
            if "DataBase" not in nested_entries:
                raise ValueError(f"Missing DBS database payload: {path}/DataBase")
            with nested.open(nested_entries["DataBase"]) as payload:
                yield payload
    else:
        raise ValueError(f"Missing DBS database payload: {path}/DataBase")


def read_dbs(source, *, database_ids=None) -> list[dict]:
    """Return selected databases, or all databases when database_ids is None.

    database_ids accepts one ID string or an iterable of IDs. Unknown IDs raise
    KeyError; an empty iterable selects nothing. Matching IDs in multiple kinds
    are all returned, in storage order. Only selected database entries are read.
    References retain C# names and stored Spectrum peak arrays. Proteomics/EAD
    payloads are Data bytes with References=None. Streams must be seekable and
    remain open when caller-owned. Flat and legacy nested ZIPs are supported.
    """
    with open_binary_source(source) as stream, zipfile.ZipFile(stream) as archive:
        entries = _entry_index(archive)
        databases = _databases(_storage(archive, entries), database_ids)
        for database in databases:
            database["References"] = None
            with _database_stream(archive, entries, database) as payload:
                if database["Kind"] == "MetabolomicsDB":
                    database["References"] = [named(row, "MoleculeMsReference")
                                              for row in iter_large_list(payload)]
                else:
                    database["Data"] = payload.read()
        return databases


def iter_dbs_references(source, *, database_id):
    """Yield references from one molecule database without collecting all rows.

    Only Storage and the selected database payload are read. LZ4 decompression
    still needs a complete block; a large single block can use substantial RAM.
    Rows have the same schema as read_dbs()[i]["References"][j]. IDs and ordering
    are unchanged. Missing IDs raise KeyError; non-molecule-only IDs raise
    ValueError. Errors are raised during iteration, and the unread remainder is
    not validated on early termination. Exhaust or close the generator to release
    resources; caller-owned source streams remain open.
    """
    if not isinstance(database_id, str):
        raise TypeError("database_id must be a string")
    with open_binary_source(source) as stream, zipfile.ZipFile(stream) as archive:
        entries = _entry_index(archive)
        databases = _databases(_storage(archive, entries), [database_id])
        molecules = [db for db in databases if db["Kind"] == "MetabolomicsDB"]
        if len(molecules) != 1:
            raise ValueError("Expected one molecule database for database_id")
        with _database_stream(archive, entries, molecules[0]) as payload:
            for row in iter_large_list(payload):
                yield named(row, "MoleculeMsReference")
