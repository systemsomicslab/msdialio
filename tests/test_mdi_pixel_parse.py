import io
import struct
from array import array

import pytest

from msdialio import list_mdi_pixel, read_mdi_pixel


def pixel_file(rows=((100.5, 1.2, [1., 2., 3.]), (200.5, 2.4, [4., 5., 6.]))):
    return struct.pack('<i', len(rows)) + b''.join(
        struct.pack('<ddi', mz, drift, len(values))
        + struct.pack(f'<{len(values)}d', *values)
        for mz, drift, values in rows
    )


class TrackingStream(io.BytesIO):
    def __init__(self, data):
        super().__init__(data)
        self.reads = []

    def read(self, size=-1):
        self.reads.append((self.tell(), size))
        return super().read(min(size, 7))  # exercise short reads


def test_catalog_skips_all_intensities():
    stream = TrackingStream(pixel_file())
    catalog = list_mdi_pixel(stream)
    assert catalog == [
        dict(image_index=0, mz=100.5, drift=1.2, pixel_count=3),
        dict(image_index=1, mz=200.5, drift=2.4, pixel_count=3),
    ]
    assert not stream.closed
    assert all(start < 24 or 48 <= start < 68 for start, _ in stream.reads)


def test_selection_preserves_order_and_deduplicates(tmp_path):
    path = tmp_path / 'sample.mdi_pixel'
    path.write_bytes(pixel_file())
    images = read_mdi_pixel(path, image_indices=(i for i in [1, 0, 1]))
    assert [image['image_index'] for image in images] == [1, 0]
    assert images[0]['intensities'] == array('d', [4., 5., 6.])
    assert images[1]['intensities'] == array('d', [1., 2., 3.])
    assert images[0]['mz'] == 200.5


def test_selection_skips_unselected_intensities_and_borrows_stream():
    stream = TrackingStream(b'prefix' + pixel_file())
    stream.seek(6)
    image, = read_mdi_pixel(stream, image_indices=[1])
    assert list(image['intensities']) == [4., 5., 6.]
    assert not stream.closed
    assert all(not 30 <= start < 54 for start, _ in stream.reads)


def test_selection_is_required():
    with pytest.raises(TypeError):
        read_mdi_pixel(io.BytesIO(pixel_file()))


@pytest.mark.parametrize('indices', [[-1], [2]])
def test_invalid_indices(indices):
    with pytest.raises(ValueError, match='out of range'):
        read_mdi_pixel(io.BytesIO(pixel_file()), image_indices=indices)


@pytest.mark.parametrize('indices', [[True], [1.5], ['1'], None])
def test_noninteger_selection(indices):
    with pytest.raises(TypeError):
        read_mdi_pixel(io.BytesIO(pixel_file()), image_indices=indices)


def test_empty_files_and_selection():
    assert list_mdi_pixel(io.BytesIO(pixel_file(()))) == []
    assert read_mdi_pixel(io.BytesIO(pixel_file(())), image_indices=[]) == []
    assert read_mdi_pixel(io.BytesIO(pixel_file()), image_indices=[]) == []
    rows = [(100., 0., [])]
    image, = read_mdi_pixel(io.BytesIO(pixel_file(rows)), image_indices=[0])
    assert image['intensities'] == array('d')


@pytest.mark.parametrize('data', [
    b'', struct.pack('<i', -1), pixel_file()[:-1], pixel_file() + b'x',
    struct.pack('<iddi', 1, 100., 0., -1),
])
def test_corrupt_layout(data):
    with pytest.raises(ValueError):
        list_mdi_pixel(io.BytesIO(data))


def test_inconsistent_pixel_counts():
    data = bytearray(pixel_file())
    struct.pack_into('<i', data, 64, 2)
    with pytest.raises(ValueError, match='Inconsistent'):
        list_mdi_pixel(io.BytesIO(data))
    with pytest.raises(ValueError, match='Inconsistent'):
        read_mdi_pixel(io.BytesIO(data), image_indices=[1])


def test_offsets_beyond_four_gb_without_large_file():
    pixels = 600_000_000
    stride = 20 + pixels * 8

    class VirtualFile:
        position = 0

        def tell(self):
            return self.position

        def seek(self, offset, whence=0):
            assert whence in (0, 2)
            self.position = offset if whence == 0 else 4 + 2 * stride + offset

        def read(self, size):
            segments = {0: struct.pack('<i', 2),
                        4: struct.pack('<ddi', 100., 0., pixels),
                        4 + stride: struct.pack('<ddi', 200., 0., pixels)}
            result = segments[self.position]
            assert len(result) == size
            self.position += size
            return result

    assert list_mdi_pixel(VirtualFile())[1]['mz'] == 200.


def test_selected_intensity_beyond_four_gb():
    count = 200_000_000
    index = count - 1

    class VirtualFile:
        position = 0

        def tell(self):
            return self.position

        def seek(self, offset, whence=0):
            self.position = offset if whence == 0 else 4 + count * 28 + offset

        def read(self, size):
            segments = {
                0: struct.pack('<i', count),
                4: struct.pack('<ddi', 100., 0., 1),
                4 + index * 28: struct.pack('<ddi', 200., 0., 1),
                24 + index * 28: struct.pack('<d', 42.),
            }
            result = segments[self.position]
            assert len(result) == size
            self.position += size
            return result

    image, = read_mdi_pixel(VirtualFile(), image_indices=[index])
    assert image['intensities'] == array('d', [42.])


def test_reads_intensities_in_bounded_chunks():
    values = [1.25] * 150_000
    data = pixel_file([(100., 0., values)])

    class BoundedStream(io.BytesIO):
        def read(self, size=-1):
            assert 0 <= size <= 1024 * 1024
            return super().read(size)

    image, = read_mdi_pixel(BoundedStream(data), image_indices=[0])
    assert len(image['intensities']) == len(values)
    assert image['intensities'][0] == image['intensities'][-1] == 1.25
