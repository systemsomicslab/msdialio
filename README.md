# msdialio

[English](README.md) | [日本語](README.ja.md)

A Python package for reading MS-DIAL `.arf`, `.arf2`, `.pai2`, `.dcl`, `.EIC.aef`, `.mddata`, and `.dbs` output files.

## Installation

### Clone the repository

Install the package and its required dependencies with pip:

```bash
python -m pip install .
```

Or use uv from the repository root:

```bash
uv sync
```

### Install directly from GitHub

Replace `OWNER/REPO` with the GitHub repository location.

Standard install with pip:

```bash
python -m pip install "msdialio @ git+https://github.com/OWNER/REPO.git"
```

With uv, install into an existing virtual environment:

```bash
uv pip install "msdialio @ git+https://github.com/OWNER/REPO.git"
```

## Python API

```python
from msdialio import read_arf

spots = read_arf("AlignmentResult_PeakProperties.arf")
```

The package-level API consists of readers: `read_arf`, `read_arf2`, `read_pai2`, `read_dcl`, `read_aef`, `read_mddata`, and `read_dbs`. `read_aef` accepts optional `spot_ids` and `file_ids` to read selected chromatograms directly.

## Command line

```bash
python -m msdialio arf input.arf --json peaks.json
python -m msdialio aef input.EIC.aef --csv spots.csv
```

Choose one of `aef`, `arf`, `arf2`, `dcl`, `pai2`, `mddata`, or `dbs`. Run `python -m msdialio --help` for details.

## Development

```bash
python -m pip install -e ".[test]"
python -m pytest
```

## License

MIT

## Reader contract

This is an alpha package; the import and distribution name is `msdialio`.
Raw writers support ARF, ARF2, PAI2, MDDATA and DBS (see below).
The CLI can export JSON and CSV.

All readers accept a path or a binary stream. Caller-owned streams remain open.
AEF and DCL streams must support seeking. Results are lists of dictionaries, except `read_mddata`, which returns one dictionary.

| Reader | One result represents | Options |
| --- | --- | --- |
| `read_arf2` | An alignment spot's stored representative measurements | `slim=True` |
| `read_arf` | A spot ordinal (`SpotIndex`) and its `SamplePeaks` | `slim=True` |
| `read_pai2` | A peak in one measurement | `slim=True` |
| `read_dcl` | A deconvoluted spectrum | `dcl_indices=None` |
| `read_aef` | A spot and its sample chromatograms | `spot_ids=None`, `file_ids=None` |

For ARF/ARF2/PAI2, `slim=False` includes additional annotations and metadata;
it does not expose every serialized field. Public fields use PascalCase, with
`RT`, `RI`, `DriftTime`, `Mz`, `PeakHeightTop` and `PeakAreaAboveZero` in sample
peak records. ARF2 retains spot-specific names such as `MassCenter` and
`HeightAverage`; its current time field is `RT` only. DCL and AEF currently retain
their existing snake_case fields. Unset ARF/PAI2 time coordinates become `None`.

ARF does not manufacture representative values or alignment identifiers.
`SpotIndex` is the file ordinal; use ARF2 for stored representative measurements.
Match ordinals and identifiers only after checking the relationship for your data.
MS/MS spectra are normally in the corresponding DCL file, even when PAI2 has
MS/MS references. Spectrum pairs are `[mz, intensity]`.

```python
from msdialio import read_arf, read_pai2

spots = read_arf("alignment_PeakProperties.arf")
first_sample = spots[0]["SamplePeaks"][0]
peaks = read_pai2("sample.pai2", slim=False)
```

## Raw MessagePack API

```python
from msdialio.raw import read_raw_arf, read_raw_arf2, read_raw_pai2

blocks = read_raw_arf("alignment_PeakProperties.arf")
objects = read_raw_arf2("alignment.arf2")
peaks = read_raw_pai2("sample.pai2")
```

Raw readers remove compression envelopes and decode MessagePack, retaining
payload nesting, unknown fields, numeric map keys and sentinel values. They do
not rename fields, select peak rows, normalize missing values or accept `slim`.
They are only exported from `msdialio.raw`.

- `read_raw_arf`: a list of blocks, each a list of decoded payload objects.
- `read_raw_arf2`: a list of decoded payload objects, even when there is only one.
- `read_raw_pai2`: the single decoded payload object, normally a peak array.

MessagePack strings become Python strings; binary values remain bytes. Invalid
UTF-8 uses `surrogateescape` so bytes can be recovered by encoding with that
error handler. Incomplete MessagePack objects raise an error. Raw reads retain
all decoded data in memory. DCL/AEF are not MessagePack and have no Raw API here.

## MS-DIAL 5 dataset and database files

```python
from msdialio import read_mddata, read_dbs

dataset = read_mddata("Dataset.mddata")
files = dataset["AnalysisFiles"]
parameters = dataset["Parameter"]
project = parameters["ProjectParam"]

databases = read_dbs("Dataset_Loaded.msp2.dbs")
for database in databases:
    for reference in database["References"] or []:
        name = reference["Name"]
        precursor_mz = reference["PrecursorMz"]
        spectrum = reference["Spectrum"]  # stored peaks: [mass, intensity, ...]
```

`read_mddata` returns `AnalysisFiles`, `AlignmentFiles`, embedded reference lists,
`Parameter`, and stored database metadata. File records and base parameter sections
use the C# serialized field names, including `AnalysisFileId`, `AnalysisFileClass`,
`PeakAreaBeanInformationFilePath`, `ProjectParam`, and `PeakPickBaseParam`.
It reads only the selected file: stored paths are neither followed nor rewritten,
and companion DBS/MSP2 libraries are not loaded automatically. A saved dataset's
embedded `MspDB` can be empty even when a companion library contains references.

`read_dbs` returns one dictionary per database with `DataBaseID`, `Kind`,
`DataBase`, `Pairs`, and `References`. Molecule database metadata includes `Id`,
`SourceType`, `DataBaseSource`, and `DataBaseSourceFilePath`. Molecule reference
records include `ScanID`, `PrecursorMz`, `Name`, `Formula`, `SMILES`, `InChIKey`,
`ChromXs`, `AdductType`, and `Spectrum`, where present in the source.
`ScanID` is scoped to its database; retain `DataBaseID` when joining annotations.
Both current flat ZIP entries and older nested ZIP entries are accepted, including
Windows backslash separators. No files are extracted to disk.

These readers target **MS-DIAL 5**. Numeric enums, nested chemical properties
(`Formula`, `AdductType`, `ChromXs`), annotator pairs, and unrecognized nested
parameter types retain their serialized representation. Unknown numeric fields
are retained with numeric keys; specialized parameters beyond the base contract
are not assigned guessed names. Proteomics and EAD lipid database payloads are
preserved in `Data` as bytes with `References=None`; their records are not decoded.
DBS input streams must support seeking. `read_dbs` collects selected databases in
memory; `read_mddata` collects the dataset in memory.
Raw writers for these formats are available from `msdialio.raw`.

### Select databases and iterate references

Use `list_dbs` to discover IDs before selecting a database. It reads only the
ZIP directory and `Storage`, without opening database or annotator payloads.
Each result contains `DataBaseID`, `Kind`, `DataBase`, and `Pairs`, with the same
metadata schema and order as `read_dbs`. `References`, `Data`, and compound counts
are omitted. This is a metadata listing; missing or corrupt payloads are not
validated until they are read. Inputs can be paths or seekable binary streams;
caller-owned streams remain open.

```python
from msdialio import list_dbs, read_dbs

path = "Dataset_Loaded.msp2.dbs"
catalog = list_dbs(path)
for database in catalog:
    print(database["DataBaseID"], database["Kind"], database["DataBase"])

if catalog:
    selected = read_dbs(path, database_ids=[catalog[0]["DataBaseID"]])
```

```python
from contextlib import closing
from msdialio import read_dbs, iter_dbs_references

# Read one database or several; the default still reads all databases.
databases = read_dbs("Dataset_Loaded.msp2.dbs", database_ids=["MspDB"])

# Iterate one molecule database without accumulating all Python records.
with closing(iter_dbs_references("Dataset_Loaded.msp2.dbs", database_id="MspDB")) as references:
    for reference in references:
        print(reference["Name"], reference["PrecursorMz"])
```

`database_ids` accepts a single string or an iterable. Empty selections return
`[]`; unknown IDs raise `KeyError`. Output retains storage order. If an ID occurs
in different database kinds, `read_dbs` returns all matching kinds. Only `Storage`
and selected payload entries are read; unrelated databases and annotator payloads
are not opened. The raw DBS reader continues to read every entry.

`iter_dbs_references` requires one molecule database ID and yields the same
reference dictionaries as `read_dbs`. Other database kinds are not decoded by
this iterator. Errors occur when advancing the iterator. Exhaust it or use
`contextlib.closing` when stopping early to release its open archive resources;
caller-owned source streams remain open.

Iteration retains one decompressed LZ4 block and decodes records one at a time.
It avoids accumulating the full reference list, but cannot seek directly to a
compound. A database stored as a single large block still needs that entire
block decompressed in memory, so iteration may offer only a modest memory saving.
Stopping early skips subsequent blocks and does not validate the unread remainder.

```bash
python -m msdialio dbs Dataset_Loaded.msp2.dbs --database-id MspDB --json selected.json
# Repeat --database-id to select several databases.
```

```python
from msdialio.raw import read_raw_mddata, read_raw_dbs

storage = read_raw_mddata("Dataset.mddata")  # original array / numeric keys
archive = read_raw_dbs("Dataset_Loaded.msp2.dbs")
metadata = archive["Storage"]               # decoded MessagePack map
entries = archive["Entries"]                # normalized ZIP names -> original bytes
```

The raw DBS API retains every ZIP entry, including unknown and annotator payloads;
ZIP separators are normalized to `/`. The raw dataset API removes only the LZ4
envelope. Both retain unknown serialized fields.

```bash
python -m msdialio mddata Dataset.mddata --json dataset.json
python -m msdialio dbs Dataset_Loaded.msp2.dbs --json databases.json
```

Dataset CSV export writes one row with nested fields encoded as JSON. The existing
CLI JSON fallback converts binary and extension values to strings; use the raw
Python API when preserving those values exactly is required.

Wire contracts were checked against the MS-DIAL source:
[MsdialDataStorage](https://github.com/systemsomicslab/MsdialWorkbench/blob/master/src/MSDIAL5/MsdialCore/DataObj/MsdialDataStorage.cs),
[DataBaseStorage](https://github.com/systemsomicslab/MsdialWorkbench/blob/master/src/MSDIAL5/MsdialCore/DataObj/DataBaseStorage.cs),
[DataBaseItem](https://github.com/systemsomicslab/MsdialWorkbench/blob/master/src/MSDIAL5/MsdialCore/DataObj/DataBaseItem.cs), and
[LargeListMessagePack](https://github.com/systemsomicslab/MsdialWorkbench/blob/master/src/Common/CommonStandard/MessagePack/LargeListMessagePack.cs).

## Raw writers

`msdialio.raw` exports `write_raw_arf`, `write_raw_arf2`, `write_raw_pai2`,
`write_raw_mddata`, and `write_raw_dbs`. Pass the corresponding raw reader's
result, not the interpreted `read_*` dictionaries (even with `slim=False`).

```python
from msdialio.raw import read_raw_pai2, write_raw_pai2

peaks = read_raw_pai2("sample.pai2")
# Edit serialized field values as needed.
write_raw_pai2("sample_modified.pai2", peaks)
```

Every writer accepts `(destination, value, *, overwrite=False)` and returns
`None`. Destinations can be paths or binary streams. Path writes use a temporary
file in the destination directory and publish only after success; existing paths
require `overwrite=True`. On POSIX, no-clobber publication requires hard-link
support on the destination filesystem. Caller-owned streams remain open, write at their
current position, are not truncated, and may contain partial output on failure.
The `overwrite` option applies only to paths. Encoding currently holds the output
in memory for MessagePack formats.

Unknown fields, numeric keys, binary/extension values, sentinel values and
surrogate-escaped strings are retained. Raw readers preserve float32 values as
`msdialio.raw.Float32` (a `float` subclass); float64 remains a normal Python float.
Writers retain these widths because MS-DIAL's C# float32 reader rejects float64.
When replacing a float32 field, use `Float32(new_value)`; arithmetic or converting
to ordinary floats discards the width marker. New ordinary floats encode as
float64. Compression and MessagePack integer widths may change; byte-identical
output is not promised. These are low-level
writers: they do not validate C# object schemas, cross-file IDs, or stored paths.
Keep existing ARF block/object nesting and padding objects intact; inserting or
removing rows requires rebuilding its LargeList layout and is not provided here.
Noncanonical array-header widths from other producers are not retained.

For DBS, `value["Storage"]` is authoritative: `Entries["Storage"]` is regenerated.
All other entry bytes, including unknown entries, annotators and nested legacy
archives, are copied unchanged. Raw API entry separators are normalized to `/`;
writers emit Windows backslashes in ZIP entry names for MS-DIAL's lookup.
Collisions are rejected. ZIP timestamps and compression metadata are not preserved. Editing
decoded database reference records is not part of this API.

Tests cover real ARF/ARF2/PAI2 fixture round trips, synthetic and optional local
MDDATA/DBS round trips, float32/float64 preservation, ZIP entry names, and LZ4/ARF
wire layouts. Edited sample files were also checked in the MS-DIAL GUI; this is
sample-based verification, not a guarantee for every version or dataset. The
optional [C# interoperability probe](tools/msdial-interop/README.md) checks the
actual deserializers and LCMS dataset initialization. DCL/AEF writers remain
unsupported.
