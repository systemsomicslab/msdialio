import builtins

import pytest

from msdial_parsers.arf_parse import intensity_matrix


def test_intensity_matrix_explains_how_to_install_numpy(monkeypatch):
    original_import = builtins.__import__

    def import_without_numpy(name, *args, **kwargs):
        if name == "numpy":
            raise ModuleNotFoundError("No module named 'numpy'")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", import_without_numpy)

    with pytest.raises(ImportError, match=r'msdial-parsers\[matrix\]'):
        intensity_matrix([])
