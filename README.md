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

The `intensity_matrix` utility converts parsed `.arf` data into a numeric matrix and requires NumPy. Other reader APIs work without NumPy. If NumPy is missing when you call this utility, the error message explains how to install it.

With pip:

```bash
python -m pip install ".[matrix]"
```

With uv:

```bash
uv sync --extra matrix
```

### Install directly from GitHub

Replace `OWNER/REPO` with the GitHub repository location.

Standard install with pip:

```bash
python -m pip install "msdial-parsers @ git+https://github.com/OWNER/REPO.git"
```

To use `intensity_matrix`:

```bash
python -m pip install "msdial-parsers[matrix] @ git+https://github.com/OWNER/REPO.git"
```

With uv, install into an existing virtual environment:

```bash
uv pip install "msdial-parsers @ git+https://github.com/OWNER/REPO.git"
```

To use `intensity_matrix`:

```bash
uv pip install "msdial-parsers[matrix] @ git+https://github.com/OWNER/REPO.git"
```

## Python API

```python
from msdial_parsers import read_arf, peak_table

spots = read_arf("AlignmentResult_PeakProperties.arf")
rows = peak_table(spots)
```

The package also provides `read_arf2`, `read_pai2`, `read_dcl`, `read_aef`, and helper functions in each reader module.

The original low-level decoder APIs are also preserved. For example, `msdial_parsers.arf_parse.decode(path)` returns raw decoded MessagePack objects, while `read_arf(path)` returns interpreted spot records. Equivalent `decode` functions are available in `arf2_parse`, `pai2_parse`, `dcl_parse`, and `aef_parse`.

The original standalone command line decoders are available as modules under `msdial_parsers.raw`, for example `python -m msdial_parsers.raw.arf_parse input.arf`.

## Command line

```bash
python -m msdial_parsers arf input.arf --json peaks.json
python -m msdial_parsers aef input.EIC.aef --csv spots.csv
```

Choose one of `aef`, `arf`, `arf2`, `dcl`, or `pai2`. Run `python -m msdial_parsers --help` for details.

## Development

```bash
python -m pip install -e ".[test,matrix]"
python -m pytest
```

## License

MIT
