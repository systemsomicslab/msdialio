import pytest

from msdial_parsers.dcl_parse import deserialize


def test_rejects_invalid_magic():
    with pytest.raises(ValueError, match="先頭が 'DC'"):
        deserialize(b"not a dcl")


def test_empty_dcl_summary():
    from msdial_parsers.dcl_parse import summarize

    assert summarize([]) == {"total_results": 0}
