import io
from pathlib import Path
import struct

import pytest

from msdialio import read_aef
from msdialio.raw import read_raw_aef, write_raw_aef


def test_real_file_roundtrip_and_modified_points():
    original = (Path(__file__).parent / 'data' / 'alignment_eic.aef').read_bytes()
    value = read_raw_aef(io.BytesIO(original))
    result = io.BytesIO()
    write_raw_aef(result, value)
    assert result.getvalue() == original
    sample = value['spots'][0]['samples'][0]
    sample['chromatogram'].append([1.25, 42.0])
    result = io.BytesIO()
    write_raw_aef(result, value)
    actual = read_raw_aef(io.BytesIO(result.getvalue()))
    assert actual == value
    assert read_aef(io.BytesIO(result.getvalue()))[0]['samples'][0]['chromatogram'][-1] == [1.25, 42.0]
    assert len(read_aef(io.BytesIO(result.getvalue()), spot_ids=[1])) == 1


def test_empty_header_stream_and_path(tmp_path):
    value = {'header': b'CSS1extra!', 'spots': []}
    stream = io.BytesIO(b'prefix')
    stream.seek(6)
    write_raw_aef(stream, value)
    stream.seek(6)
    assert read_raw_aef(stream) == value
    assert not stream.closed
    path = tmp_path / 'out.aef'
    write_raw_aef(path, value)
    with pytest.raises(FileExistsError):
        write_raw_aef(path, value)
    write_raw_aef(path, value, overwrite=True)
    assert read_raw_aef(path) == value


@pytest.mark.parametrize('data', [
    b'bad!' + bytes(10),
    b'CSS1' + bytes(6) + struct.pack('<i', -1),
    b'CSS1' + bytes(6) + struct.pack('<i', 1) + struct.pack('<q', 0),
    b'CSS1' + bytes(6) + struct.pack('<i', 0) + b'extra',
    b'CSS1' + bytes(6) + struct.pack('<i', 1) + struct.pack('<q', 22) + bytes(5),
])
def test_rejects_malformed(data):
    with pytest.raises(ValueError):
        read_raw_aef(io.BytesIO(data))
