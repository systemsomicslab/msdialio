"""Lossless, editable representation of version 1 DCL records."""

import io
import struct

from ._sources import open_binary_source
from .dcl_parse import _read_exact
from ._raw_write import _write_bytes


def _take(stream, size):
    return _read_exact(stream, size)


def _counted(stream, count, fmt):
    if count < 0:
        raise ValueError("Negative DCL item count")
    size = struct.calcsize(fmt)
    return [list(item) for item in struct.iter_unpack(fmt, _take(stream, count * size))]


def read_raw_dcl(source):
    """Read a version 1 DCL, preserving fixed fields and annotation as bytes.

    Spectrum rows are [mass, intensity, quality]; chromatogram rows are
    [id, coordinate, mass, intensity]; model masses are floats.
    """
    with open_binary_source(source) as stream:
        base = stream.tell()
        magic, version, annotation, count = struct.unpack('<2siBi', _take(stream, 11))
        if magic != b'DC' or version != 1 or annotation not in (0, 1) or count < 0:
            raise ValueError('Invalid or unsupported DCL header')
        offsets = [x[0] for x in _counted(stream, count, '<q')]
        records = []
        expected = 11 + count * 8
        for offset in offsets:
            if offset != expected:
                raise ValueError('Invalid DCL record offset')
            stream.seek(base + offset)
            scan = _take(stream, 60)
            if struct.unpack_from('<q', scan)[0] != offset:
                raise ValueError('Invalid DCL record seek point')
            quant = _take(stream, 40)
            scoring = _take(stream, 20)
            ns, nc, nm = struct.unpack('<3i', _take(stream, 12))
            spectrum = _counted(stream, ns, '<ddi')
            chromatogram = _counted(stream, nc, '<iddd')
            model_masses = [x[0] for x in _counted(stream, nm, '<d')]
            annotation_data = None
            if annotation:
                fixed = _take(stream, 73)
                n = struct.unpack('<i', _take(stream, 4))[0]
                if n < 0:
                    raise ValueError('Negative DCL annotation count')
                annotation_data = fixed + struct.pack('<i', n) + _take(stream, n * 4)
            records.append(dict(scan=scan, quant=quant, scoring=scoring, spectrum=spectrum,
                                chromatogram=chromatogram, model_masses=model_masses,
                                annotation=annotation_data))
            expected = stream.tell() - base
        if stream.read(1):
            raise ValueError('Unexpected trailing DCL data')
        return dict(version=version, annotation=bool(annotation), records=records)


def write_raw_dcl(destination, value, *, overwrite=False):
    """Write a representation returned by read_raw_dcl; recompute offsets."""
    if value['version'] != 1 or type(value['annotation']) is not bool:
        raise ValueError('Invalid or unsupported DCL header')
    records = value['records']
    if not isinstance(records, list):
        raise TypeError('DCL records must be a list')
    data = io.BytesIO()
    data.write(struct.pack('<2siBi', b'DC', 1, value['annotation'], len(records)))
    data.write(b'\0' * (8 * len(records)))
    offsets = []
    for record in records:
        offset = data.tell()
        offsets.append(offset)
        for key, size in [('scan', 60), ('quant', 40), ('scoring', 20)]:
            block = record[key]
            if not isinstance(block, bytes) or len(block) != size:
                raise ValueError(f'DCL {key} must be {size} bytes')
            data.write(struct.pack('<q', offset) + block[8:] if key == 'scan' else block)
        spectrum, chrom, masses = (record[k] for k in ('spectrum', 'chromatogram', 'model_masses'))
        data.write(struct.pack('<3i', len(spectrum), len(chrom), len(masses)))
        for row in spectrum:
            data.write(struct.pack('<ddi', *row))
        for row in chrom:
            data.write(struct.pack('<iddd', *row))
        for mass in masses:
            data.write(struct.pack('<d', mass))
        annotation = record['annotation']
        if value['annotation']:
            if not isinstance(annotation, bytes) or len(annotation) < 77:
                raise ValueError('Invalid DCL annotation')
            n = struct.unpack_from('<i', annotation, 73)[0]
            if n < 0 or len(annotation) != 77 + n * 4:
                raise ValueError('Invalid DCL annotation count')
            data.write(annotation)
        elif annotation is not None:
            raise ValueError('DCL annotation must be None when disabled')
    data.seek(11)
    for offset in offsets:
        data.write(struct.pack('<q', offset))
    _write_bytes(destination, data.getvalue(), overwrite)
