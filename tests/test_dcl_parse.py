import pytest

from msdialio.dcl_parse import _deserialize


def test_rejects_invalid_magic():
    with pytest.raises(ValueError, match="先頭が 'DC'"):
        _deserialize(b"not a dcl")
