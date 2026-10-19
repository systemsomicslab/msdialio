# msdialio

[English](README.md) | [日本語](README.ja.md)

A Python package for reading MS-DIAL `.arf`, `.arf2`, `.pai2`, `.dcl`, and `.EIC.aef` output files.

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

The package-level API consists of readers: `read_arf`, `read_arf2`, `read_pai2`, `read_dcl`, and `read_aef`. `read_aef` accepts optional `spot_ids` and `file_ids` to read selected chromatograms directly.

## Command line

```bash
python -m msdialio arf input.arf --json peaks.json
python -m msdialio aef input.EIC.aef --csv spots.csv
```

Choose one of `aef`, `arf`, `arf2`, `dcl`, or `pai2`. Run `python -m msdialio --help` for details.

## Development

```bash
python -m pip install -e ".[test]"
python -m pytest
```

## License

MIT

## Reader contract

This is an alpha package; the import and distribution name is `msdialio`.
Writing MS-DIAL files is not implemented yet. The CLI can export JSON and CSV.

All readers accept a path or a binary stream. Caller-owned streams remain open.
AEF and DCL streams must support seeking. Results are lists of dictionaries.

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
