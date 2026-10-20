"""Command-line interface for the MS-DIAL readers."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from . import aef_parse, arf2_parse, arf_parse, dcl_parse, pai2_parse, mddata_parse, dbs_parse


_READERS = {
    "aef": aef_parse.read_aef,
    "arf": arf_parse.read_arf,
    "arf2": arf2_parse.read_arf2,
    "dcl": dcl_parse.read_dcl,
    "pai2": pai2_parse.read_pai2,
    "mddata": mddata_parse.read_mddata,
    "dbs": dbs_parse.read_dbs,
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Read an MS-DIAL output file")
    parser.add_argument("format", choices=sorted(_READERS))
    parser.add_argument("path", type=Path)
    parser.add_argument("--database-id", dest="database_ids", action="append",
                        help="read only this DBS database ID (repeatable; DBS only)")
    output = parser.add_mutually_exclusive_group()
    output.add_argument("--json", metavar="PATH", help="write parsed records as JSON")
    output.add_argument("--csv", metavar="PATH", help="write top-level record fields as CSV")
    args = parser.parse_args(argv)
    if args.database_ids is not None and args.format != "dbs":
        parser.error("--database-id is only supported for dbs")

    read = _READERS[args.format]
    records = read(args.path, database_ids=args.database_ids) if args.format == "dbs" else read(args.path)
    print(f"{args.path}: {1 if isinstance(records, dict) else len(records)} records")

    if args.json:
        with open(args.json, "w", encoding="utf-8") as stream:
            json.dump(records, stream, ensure_ascii=False, indent=2, default=str)
            stream.write("\n")
    elif args.csv:
        _write_csv(args.csv, [records] if isinstance(records, dict) else records)
    return 0


def _write_csv(path: str, records: list[dict]) -> None:
    if not records:
        Path(path).write_text("", encoding="utf-8")
        return
    fields = list(dict.fromkeys(key for record in records for key in record))
    with open(path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for record in records:
            writer.writerow({
                key: json.dumps(value, ensure_ascii=False, default=str)
                if isinstance(value, (dict, list, tuple)) else value
                for key, value in record.items()
            })


if __name__ == "__main__":
    sys.exit(main())
