"""Readers for binary output files produced by MS-DIAL."""

from importlib import import_module


_EXPORTS = {
    "attach_msms_to_features": "dcl_parse",
    "find_msms": "dcl_parse",
    "intensity_matrix": "arf_parse",
    "peak_table": "arf_parse",
    "read_aef": "aef_parse",
    "read_aef_spots": "aef_parse",
    "read_arf": "arf_parse",
    "read_arf2": "arf2_parse",
    "read_dcl": "dcl_parse",
    "read_pai2": "pai2_parse",
}

__all__ = list(_EXPORTS)

__version__ = "0.1.0"


def __getattr__(name: str):
    """Load public readers on demand so raw modules can run with ``-m``."""
    module_name = _EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(f".{module_name}", __name__), name)
    globals()[name] = value
    return value
