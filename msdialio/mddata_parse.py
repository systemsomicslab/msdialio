"""Read MS-DIAL 5 dataset files (not the older MS-DIAL 4 project format)."""

from ._metadata import decode, named
from ._sources import open_binary_source

_FIELDS = ("AnalysisFiles", "AlignmentFiles", "MspDB", "TextDB", "IsotopeTextDB",
           "IupacDatabase", "Parameter", "DataBaseMapper", "DataBases")


def read_raw_mddata(source):
    """Return the decoded storage object, preserving numeric keys and nesting."""
    with open_binary_source(source) as stream:
        return decode(stream.read())


def read_mddata(source) -> dict:
    """Return dataset metadata with named file and base parameter fields.

    Does not follow stored paths or load companion files. Specialized parameter
    keys and nested types without a known schema retain their serialized values.
    Enums remain integers. Caller-owned binary streams remain open.
    """
    value = read_raw_mddata(source)
    if not isinstance(value, list) or len(value) < 7:
        raise ValueError("Expected an MS-DIAL 5 dataset storage array")
    result = {_FIELDS[i] if i < len(_FIELDS) else i: item for i, item in enumerate(value)}
    for field, schema in (("AnalysisFiles", "AnalysisFileBean"),
                          ("AlignmentFiles", "AlignmentFileBean")):
        result[field] = [named(row, schema) for row in result[field] or []]
    result["Parameter"] = named(result["Parameter"], "ParameterBase")
    for field in ("MspDB", "TextDB", "IsotopeTextDB"):
        result[field] = [named(row, "MoleculeMsReference") for row in result[field] or []]
    return result
