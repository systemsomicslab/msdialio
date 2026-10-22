# Test data

These compact binary files were copied from sample MS-DIAL outputs provided by the project maintainer. Their filenames were normalized for the test suite; the binary contents are unchanged.

| Fixture | Format | Expected records | Notes |
|---|---|---:|---|
| `alignment_peaks.arf` | `.arf` | 187 spots | Single-sample aligned peaks |
| `alignment_spots.arf2` | `.arf2` | 187 spots | Aligned spot catalog |
| `alignment_msdec.dcl` | `.dcl` | 187 results | Alignment-level MSDec results |
| `alignment_eic.aef` | `.EIC.aef` | 187 spots | EIC chromatograms |
| `sample_peaks.pai2` | `.pai2` | 364 peaks | Single-sample peaks |
| `sample_msdec.dcl` | `.dcl` | 364 results | Sample-level MSDec results paired with the `.pai2` file |

The parser tests assert record counts and representative values so accidental format or decoding changes are visible.

User-supplied MDDATA and DBS originals for local writer regression tests are
retained unchanged under `local-project/`. These large files are Git-ignored;
tests that need them skip if absent. The tests check float32/float64 preservation
as well as value equality and opaque DBS entry preservation.
