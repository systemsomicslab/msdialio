# MS-DIAL Parsers

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
python -m pip install "msdial-parsers @ git+https://github.com/OWNER/REPO.git"
```

With uv, install into an existing virtual environment:

```bash
uv pip install "msdial-parsers @ git+https://github.com/OWNER/REPO.git"
```

## Python API

```python
from msdial_parsers import read_arf

spots = read_arf("AlignmentResult_PeakProperties.arf")
```

The package-level API consists of readers: `read_arf2`, `read_pai2`, `read_dcl`, `read_aef`, and `read_aef_spots` are also available. The `.arf2` reader module additionally provides `read_raw_spots` and `read_isotopic_peaks` for specialized data access.

## Command line

```bash
python -m msdial_parsers arf input.arf --json peaks.json
python -m msdial_parsers aef input.EIC.aef --csv spots.csv
```

Choose one of `aef`, `arf`, `arf2`, `dcl`, or `pai2`. Run `python -m msdial_parsers --help` for details.

## Development

```bash
python -m pip install -e ".[test]"
python -m pytest
```

## License

MIT
