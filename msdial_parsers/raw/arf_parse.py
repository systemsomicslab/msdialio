"""Raw decoder for MS-DIAL ``.arf`` files (aligned per-sample peaks).

This only undoes the binary encoding. It does not interpret any field: records come
back exactly as MessagePack stores them, i.e. as nested lists where the list index is
the ``[Key(N)]`` number of the MS-DIAL C# class.

Binary layout
-------------
The file is a sequence of MessagePack ``ExtType(99)`` objects (MessagePack-CSharp
"Lz4Block" compression). Each extension's payload is::

    <uncompressed length: MessagePack int> <LZ4 block>

Decompressing every block and concatenating the results gives a plain MessagePack
stream. For ``*_PeakProperties.arf`` that stream holds the alignment spots: a container
array ``[int, int, group, group, ...]`` (some trailing groups may also appear as
standalone top-level arrays). A group is one alignment spot and is a list of per-sample
rows; each row is an ``AlignmentChromPeakFeature`` array. ``*_DriftSopts.arf`` has the
same encoding and is normally almost empty.

Usage (Python)::

    from arf_parse import decode, decompress, iter_objects
    objects = decode("AlignmentResult_..._PeakProperties.arf")   # list of decoded objects
    raw = decompress("AlignmentResult_..._PeakProperties.arf")   # decompressed MessagePack bytes

Usage (CLI)::

    python arf_parse.py <file.arf> [--json out.json] [--msgpack out.msgpack]

Without options the CLI prints an outline of the decoded structure. Large ``.arf`` files
can decompress to around 1 GB, and their fully decoded Python objects need several times
that much RAM; write ``--msgpack`` and process it in a streaming way when that is a problem.

Requires: msgpack, lz4.
"""

from __future__ import annotations

import argparse
import base64
import io
import json

import lz4.block
import msgpack

#: MessagePack-CSharp extension type code for "Lz4Block" compression.
LZ4_BLOCK_EXT = 99


def _unpacker(data: bytes) -> msgpack.Unpacker:
    return msgpack.Unpacker(
        io.BytesIO(data),
        raw=False,
        # Some positive-mode files contain invalid UTF-8 in header strings; strict
        # decoding would abort the whole file on a single bad byte.
        unicode_errors="ignore",
        strict_map_key=False,
        max_buffer_size=2**31 - 1,
    )


def _decompress_lz4_block(payload: bytes) -> bytes:
    """Decompress one Lz4Block payload: a MessagePack int (length) followed by LZ4 data."""
    header = _unpacker(payload)
    size = header.unpack()
    return lz4.block.decompress(payload[header.tell():], uncompressed_size=size)


def _iter_top_level(data: bytes):
    """Yield ``(is_compressed, value)`` for each top-level object of the file."""
    for obj in _unpacker(data):
        if isinstance(obj, msgpack.ExtType) and obj.code == LZ4_BLOCK_EXT:
            yield True, _decompress_lz4_block(obj.data)
        else:
            yield False, obj  # stored without compression


def iter_objects(path):
    """Yield every decoded MessagePack object in the file, decompressing on the fly."""
    with open(path, "rb") as f:
        data = f.read()
    for compressed, value in _iter_top_level(data):
        if compressed:
            yield from _unpacker(value)
        else:
            yield value


def decode(path) -> list:
    """Return all decoded MessagePack objects in the file as a list."""
    return list(iter_objects(path))


def decompress(path) -> bytes:
    """Return the decompressed content as one plain MessagePack byte stream."""
    with open(path, "rb") as f:
        data = f.read()
    return b"".join(value if compressed else msgpack.packb(value, use_bin_type=True)
                    for compressed, value in _iter_top_level(data))


def _json_default(obj):
    """Encode the MessagePack types that JSON cannot represent directly."""
    if isinstance(obj, (bytes, bytearray)):
        return {"$bytes": base64.b64encode(obj).decode("ascii")}
    if isinstance(obj, msgpack.ExtType):
        return {"$ext": obj.code, "data": base64.b64encode(obj.data).decode("ascii")}
    raise TypeError(f"cannot serialize {type(obj).__name__}")


def write_json(objects: list, out_path: str) -> None:
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(objects, f, default=_json_default, ensure_ascii=False)


def describe(obj, depth: int = 2, width: int = 3) -> str:
    """Short one-line outline of a decoded object (types and lengths only)."""
    if isinstance(obj, list):
        if depth <= 0 or not obj:
            return f"list[{len(obj)}]"
        head = ", ".join(describe(x, depth - 1, width) for x in obj[:width])
        return f"list[{len(obj)}]({head}{', ...' if len(obj) > width else ''})"
    if isinstance(obj, dict):
        return f"map[{len(obj)}]"
    if isinstance(obj, (bytes, bytearray)):
        return f"bin[{len(obj)}]"
    if isinstance(obj, str):
        return repr(obj if len(obj) <= 24 else obj[:24] + "...")
    return repr(obj)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Decode an MS-DIAL .arf file into raw MessagePack objects.")
    ap.add_argument("path", help="path to the .arf file")
    ap.add_argument("--json", help="write the decoded objects as JSON to this path")
    ap.add_argument("--msgpack", help="write the decompressed MessagePack stream to this path")
    args = ap.parse_args(argv)

    if args.msgpack:
        with open(args.msgpack, "wb") as f:
            f.write(decompress(args.path))
        print(f"wrote decompressed MessagePack: {args.msgpack}")

    objects = decode(args.path)
    print(f"{args.path}: {len(objects)} top-level object(s)")
    for i, obj in enumerate(objects[:10]):
        print(f"  [{i}] {describe(obj)}")
    if args.json:
        write_json(objects, args.json)
        print(f"wrote JSON: {args.json}")


if __name__ == "__main__":
    main()
