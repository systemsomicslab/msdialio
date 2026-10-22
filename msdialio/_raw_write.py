"""Writers for the decoded representations exposed by msdialio.raw."""

import os
from pathlib import Path
import struct
import tempfile
import zipfile

import lz4.block
import msgpack

from ._messagepack_objects import Float32


def _pack(value):
    packer = msgpack.Packer(use_bin_type=True, unicode_errors="surrogateescape")
    result = bytearray()

    def append(item):
        if isinstance(item, Float32):
            result.extend(b"\xca" + struct.pack(">f", item))
        elif isinstance(item, msgpack.ExtType):
            result.extend(packer.pack(item))
        elif isinstance(item, (list, tuple)):
            result.extend(packer.pack_array_header(len(item)))
            for child in item:
                append(child)
        elif isinstance(item, dict):
            result.extend(packer.pack_map_header(len(item)))
            for key, child in item.items():
                append(key)
                append(child)
        else:
            result.extend(packer.pack(item))

    append(value)
    return bytes(result)


def _compressed(payload):
    # MS-DIAL's LargeList reader requires ext32 and a fixed-width int32 size.
    body = b"\xd2" + struct.pack(">i", len(payload))
    body += lz4.block.compress(payload, store_size=False)
    return b"\xc9" + struct.pack(">I", len(body)) + b"\x63" + body


def _save(destination, emit, overwrite):
    if hasattr(destination, "write"):
        emit(destination)
        return
    path = Path(destination)
    if not overwrite and path.exists():
        raise FileExistsError(path)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            emit(stream)
        if overwrite:
            os.replace(temporary, path)
        elif os.name == "nt":
            # Unlike POSIX rename, Windows rename refuses an existing target.
            os.rename(temporary, path)
        else:
            # An atomic no-clobber publication, including a concurrent creator.
            os.link(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _write_bytes(destination, data, overwrite):
    def emit(stream):
        remaining = memoryview(data)
        while remaining:
            written = stream.write(remaining)
            if written is None or written <= 0:
                raise OSError("Output stream did not accept data")
            remaining = remaining[written:]
    _save(destination, emit, overwrite)


def write_raw_pai2(destination, value, *, overwrite=False):
    """Write the single object returned by read_raw_pai2."""
    _write_bytes(destination, _compressed(_pack(value)), overwrite)


def write_raw_mddata(destination, value, *, overwrite=False):
    """Write the storage object returned by read_raw_mddata."""
    _write_bytes(destination, _compressed(_pack(value)), overwrite)


def write_raw_arf2(destination, objects, *, overwrite=False):
    """Write the ordered payload objects returned by read_raw_arf2."""
    if not isinstance(objects, list):
        raise TypeError("ARF2 objects must be a list from read_raw_arf2")
    _write_bytes(destination, _compressed(b"".join(_pack(obj) for obj in objects)), overwrite)


def write_raw_arf(destination, blocks, *, overwrite=False):
    """Write the block/object structure returned by read_raw_arf.

    Retain the existing nesting and padding objects when editing field values.
    This interface does not rebuild the LargeList layout after inserting rows.
    """
    if not isinstance(blocks, list) or any(not isinstance(block, list) for block in blocks):
        raise TypeError("ARF blocks must be lists of objects from read_raw_arf")
    data = b"".join(_compressed(b"".join(_pack(obj) for obj in block)) for block in blocks)
    _write_bytes(destination, data, overwrite)


def write_raw_dbs(destination, value, *, overwrite=False):
    """Write a raw DBS archive, regenerating Storage from value['Storage'].

    All other Entries are copied as bytes, including unknown/legacy payloads.
    ZIP metadata and compression are not preserved. Entries use Windows
    separators on disk for MS-DIAL's Path.Combine/GetEntry lookup.
    """
    if not isinstance(value, dict) or not isinstance(value.get("Storage"), dict):
        raise TypeError("DBS requires a Storage map")
    if not isinstance(value.get("Entries"), dict):
        raise TypeError("DBS requires an Entries map")
    entries = {}
    for name, payload in value["Entries"].items():
        if not isinstance(name, str) or not name or "\x00" in name:
            raise ValueError("Invalid DBS entry name")
        name = name.replace("\\", "/")
        if name.endswith("/") or name in entries:
            raise ValueError(f"Invalid or duplicate DBS entry: {name}")
        if not isinstance(payload, bytes):
            raise TypeError(f"DBS entry must contain bytes: {name}")
        entries[name] = payload
    entries["Storage"] = _compressed(_pack(value["Storage"]))

    def emit(stream):
        with zipfile.ZipFile(stream, "w") as archive:
            for name, payload in entries.items():
                info = zipfile.ZipInfo(name)
                # Assign after construction: ZipInfo normalizes os.sep to '/'
                # on Windows, but MS-DIAL looks up literal backslash names.
                info.filename = name.replace("/", "\\")
                archive.writestr(info, payload)
    _save(destination, emit, overwrite)
