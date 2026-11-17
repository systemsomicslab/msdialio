"""Editable CSS1 chromatograms, including fields omitted by read_aef."""

import io
import struct

from ._raw_write import _write_bytes
from ._sources import open_binary_source
from .aef_parse import _read_exact


def _take(stream, size):
    return _read_exact(stream, size, 'raw AEF data')


def read_raw_aef(source):
    """Read all CSS1 spots and samples, retaining peak bounds and header."""
    with open_binary_source(source) as stream:
        base = stream.tell()
        header = _take(stream, 10)
        if header[:4] != b'CSS1':
            raise ValueError('Not a CSS1 file')
        count = struct.unpack('<i', _take(stream, 4))[0]
        if count < 0:
            raise ValueError('Negative AEF spot count')
        pointers = [x[0] for x in struct.iter_unpack('<q', _take(stream, count * 8))]
        expected = 14 + count * 8
        spots = []
        for pointer in pointers:
            if pointer != expected:
                raise ValueError('Invalid AEF spot pointer')
            rt, ri, mz, drift, main_type, n = struct.unpack('<ffffBi', _take(stream, 21))
            if n < 0:
                raise ValueError('Negative AEF sample count')
            samples = []
            for _ in range(n):
                file_id, points, top, left, right = struct.unpack('<iifff', _take(stream, 20))
                if points < 0:
                    raise ValueError('Negative AEF point count')
                chromatogram = [list(p) for p in struct.iter_unpack('<ff', _take(stream, points * 8))]
                samples.append(dict(file_id=file_id, peak_top=top, peak_left=left,
                                    peak_right=right, chromatogram=chromatogram))
            spots.append(dict(rt=rt, ri=ri, mz=mz, drift=drift, main_type=main_type,
                              samples=samples))
            expected = stream.tell() - base
        if stream.read(1):
            raise ValueError('Unexpected trailing AEF data')
        return dict(header=header, spots=spots)


def write_raw_aef(destination, value, *, overwrite=False):
    """Write CSS1 spots and regenerate counts and seek pointers."""
    header = value['header']
    if not isinstance(header, bytes) or len(header) != 10 or header[:4] != b'CSS1':
        raise ValueError('Invalid CSS1 header')
    spots = value['spots']
    if not isinstance(spots, list):
        raise TypeError('AEF spots must be a list')
    output = io.BytesIO()
    output.write(header)
    output.write(struct.pack('<i', len(spots)))
    output.write(bytes(8 * len(spots)))
    pointers = []
    for spot in spots:
        pointers.append(output.tell())
        samples = spot['samples']
        output.write(struct.pack('<ffffBi', spot['rt'], spot['ri'], spot['mz'], spot['drift'],
                                 spot['main_type'], len(samples)))
        for sample in samples:
            points = sample['chromatogram']
            output.write(struct.pack('<iifff', sample['file_id'], len(points), sample['peak_top'],
                                     sample['peak_left'], sample['peak_right']))
            for point in points:
                output.write(struct.pack('<ff', *point))
    output.seek(14)
    for pointer in pointers:
        output.write(struct.pack('<q', pointer))
    _write_bytes(destination, output.getvalue(), overwrite)
