# Optional MS-DIAL interoperability check

`VerifyMsdial.cs` uses an installed Windows MS-DIAL build to check files with its
actual C# deserializers. It requires Windows PowerShell/.NET Framework and the
MS-DIAL assembly folder; it is not part of the portable pytest suite.

Compile once in Windows PowerShell:

```powershell
Add-Type -Path tools/msdial-interop/VerifyMsdial.cs -OutputAssembly tools/msdial-interop/VerifyMsdial.exe -OutputType ConsoleApplication
```

Pass the assembly folder, then one or more format/path pairs:

```powershell
./tools/msdial-interop/VerifyMsdial.exe C:/path/to/MSDIAL pai2 sample.pai2 arf2 alignment.arf2 arf alignment_PeakProperties.arf
./tools/msdial-interop/VerifyMsdial.exe C:/path/to/MSDIAL dataset-model Dataset.mddata
```

Supported modes:

| Mode | Check |
| --- | --- |
| `pai2`, `arf2`, `arf`, `mddata`, `msp2` | Deserialize the corresponding LCMS type |
| `dbs-storage` | Deserialize an extracted DBS `Storage` entry |
| `dataset` | Load the LCMS dataset with companion DBS and MSP2 files |
| `dataset-model` | Also construct the GUI DatasetModel and restore identification settings |

No application window is shown. Dataset modes require companion files beside
the MDDATA. They do not exercise raw-data loading, analysis views, or other GUI
interactions. Compatible `Common.dll`, `MsdialCore.dll`, `MsdialLcMsApi.dll` and
their dependencies are required; `dataset-model` also uses `MSDIAL.exe` and
`ReactiveProperty.dll`. Other MS-DIAL versions may use different assembly APIs.

This check caught two defects that Python value equality alone missed:

- Replacing float32 with float64 causes C# `ReadSingle` to reject code 203.
- Forward-slash DBS entry names leave annotation factories uninitialized when
  Windows MS-DIAL looks up backslash names, failing later during settings restore.

Portable regression coverage is in `tests/test_raw_write.py`, including checks
of float tags and ZIP header bytes. Large local MDDATA/DBS originals are retained
in the Git-ignored `tests/data/local-project/`; their tests skip if absent.
