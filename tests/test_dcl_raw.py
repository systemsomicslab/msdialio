import io
import struct
from pathlib import Path

import pytest

from msdialio import read_dcl
from msdialio.raw import read_raw_dcl, write_raw_dcl

DATA = Path(__file__).parent / 'data'


@pytest.mark.parametrize('name', ['alignment_msdec.dcl', 'sample_msdec.dcl'])
def test_fixture_roundtrip_and_edit(name):
    original = (DATA / name).read_bytes()
    value = read_raw_dcl(io.BytesIO(original))
    out = io.BytesIO()
    write_raw_dcl(out, value)
    assert out.getvalue() == original
    value['records'][0]['spectrum'].append([321.25, 42.0, 3])
    out = io.BytesIO()
    write_raw_dcl(out, value)
    decoded = read_raw_dcl(io.BytesIO(out.getvalue()))
    assert decoded['records'][0]['spectrum'][-1] == [321.25, 42.0, 3]
    assert decoded['records'][1]['spectrum'] == value['records'][1]['spectrum']
    assert read_dcl(io.BytesIO(out.getvalue()))[0]['msms_spectrum'][-1] == [321.25, 42.0]


def test_annotation_empty_and_stream_offset(tmp_path):
    value = {'version': 1, 'annotation': True, 'records': [{
        'scan': bytes(60), 'quant': bytes(40), 'scoring': bytes(20),
        'spectrum': [], 'chromatogram': [[1, 2., 3., 4.]], 'model_masses': [5.],
        'annotation': bytes(73) + struct.pack('<iii', 2, 11, 12),
    }]}
    stream = io.BytesIO(b'prefix')
    stream.seek(6)
    write_raw_dcl(stream, value)
    stream.seek(6)
    result = read_raw_dcl(stream)
    assert result['records'][0]['annotation'] == value['records'][0]['annotation']
    assert result['records'][0]['chromatogram'] == value['records'][0]['chromatogram']
    assert not stream.closed
    path = tmp_path / 'out.dcl'
    write_raw_dcl(path, value)
    with pytest.raises(FileExistsError):
        write_raw_dcl(path, value)
    write_raw_dcl(path, value, overwrite=True)
    empty = {'version': 1, 'annotation': False, 'records': []}
    write_raw_dcl(path, empty, overwrite=True)
    assert read_raw_dcl(path) == empty


def test_rejects_invalid_layout():
    with pytest.raises(ValueError):
        read_raw_dcl(io.BytesIO(struct.pack('<2siBi', b'DC', 2, 0, 0)))
    with pytest.raises(ValueError):
        read_raw_dcl(io.BytesIO(struct.pack('<2siBi', b'DC', 1, 0, 1) + struct.pack('<q', 0)))
    with pytest.raises(ValueError):
        read_raw_dcl(io.BytesIO(struct.pack('<2siBi', b'DC', 1, 0, 0) + b'x'))
