"""Image and label arrays on disk. Needs neither cellprofiler_core nor the recorded fixtures."""
import os

import h5py
import numpy
import pytest

from scverse_export.raster import (DATASET, array_info, read_array, write_array, write_image,
                                   write_labels)


@pytest.mark.parametrize("dtype", ["uint8", "uint16", "int32", "int64", "float32", "float64"])
def test_round_trip_preserves_values_and_dtype(tmp_path, dtype):
    """Source dtype has to survive: a uint16 image silently widened to float64 would quadruple
    the store and lose the fact that it was integer intensity."""
    data = (numpy.arange(24).reshape(2, 3, 4)).astype(dtype)
    path = write_array(str(tmp_path / "a.h5"), data)
    back = read_array(path)
    numpy.testing.assert_array_equal(back, data)
    assert back.dtype == data.dtype


def test_dataset_is_named_data(tmp_path):
    """The importer looks up f["data"] by name, so this is part of the on-disk contract."""
    write_array(str(tmp_path / "a.h5"), numpy.zeros((2, 2)))
    with h5py.File(tmp_path / "a.h5", "r") as f:
        assert list(f.keys()) == [DATASET]


def test_creates_missing_parent_directories(tmp_path):
    path = str(tmp_path / "labels" / "P1_A01_1" / "Nuclei.h5")
    write_labels(path, numpy.zeros((4, 4), dtype="int32"))
    assert os.path.exists(path)


def test_no_temporary_files_are_left_behind(tmp_path):
    write_array(str(tmp_path / "a.h5"), numpy.ones((8, 8)))
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a.h5"]


def test_failed_write_leaves_the_previous_file_intact(tmp_path, monkeypatch):
    """The reason for temp-and-rename. A crash mid-write must not leave a truncated file that
    reads as valid, and must not destroy what was already there."""
    path = str(tmp_path / "a.h5")
    write_array(path, numpy.full((4, 4), 7, dtype="int32"))

    import scverse_export.raster as raster
    real_replace = os.replace

    def boom(src, dst):
        raise OSError("disk full")

    monkeypatch.setattr(raster.os, "replace", boom)
    with pytest.raises(OSError):
        write_array(path, numpy.full((4, 4), 9, dtype="int32"))

    monkeypatch.setattr(raster.os, "replace", real_replace)
    numpy.testing.assert_array_equal(read_array(path), numpy.full((4, 4), 7, dtype="int32"))
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a.h5"], "temp file was not cleaned up"


def test_overwrites_in_place(tmp_path):
    path = str(tmp_path / "a.h5")
    write_array(path, numpy.zeros((3, 3), dtype="int32"))
    write_array(path, numpy.ones((5, 5), dtype="int32"))
    assert read_array(path).shape == (5, 5)


def test_array_info_reports_shape_and_dtype_without_reading(tmp_path):
    path = write_image(str(tmp_path / "a.h5"), numpy.zeros((4, 16, 32), dtype="uint16"))
    assert array_info(path) == ((4, 16, 32), "uint16")


def test_compression_is_on_by_default_and_can_be_turned_off(tmp_path):
    """Label arrays are mostly zeros, so the default earns its CPU. Checked as a real size
    difference rather than by asserting the filter is set."""
    labels = numpy.zeros((512, 512), dtype="int32")
    labels[100:140, 100:140] = 1
    labels[200:260, 200:260] = 2
    small = write_labels(str(tmp_path / "z.h5"), labels)
    big = write_labels(str(tmp_path / "u.h5"), labels, compression=None)
    assert os.path.getsize(small) < os.path.getsize(big) / 10
    numpy.testing.assert_array_equal(read_array(small), read_array(big))


def test_empty_array_round_trips(tmp_path):
    """A field of view where an object type found nothing still has to produce a file, and h5py
    cannot chunk (so cannot compress) a zero-size dataset."""
    path = write_array(str(tmp_path / "a.h5"), numpy.zeros((0, 4), dtype="int32"))
    back = read_array(path)
    assert back.shape == (0, 4) and back.dtype == numpy.dtype("int32")


def test_all_zero_label_array_round_trips(tmp_path):
    """The commoner empty case: the array is full size, every pixel background."""
    path = write_labels(str(tmp_path / "n.h5"), numpy.zeros((64, 64), dtype="int32"))
    assert read_array(path).max() == 0


def test_write_image_requires_a_channel_axis(tmp_path):
    with pytest.raises(ValueError, match=r"\(C, Y, X\)"):
        write_image(str(tmp_path / "a.h5"), numpy.zeros((16, 16), dtype="uint16"))
    write_image(str(tmp_path / "b.h5"), numpy.zeros((1, 16, 16), dtype="uint16"))


def test_write_labels_rejects_wrong_shape_and_dtype(tmp_path):
    with pytest.raises(ValueError, match=r"\(Y, X\)"):
        write_labels(str(tmp_path / "a.h5"), numpy.zeros((2, 4, 4), dtype="int32"))
    with pytest.raises(ValueError, match="integers"):
        write_labels(str(tmp_path / "b.h5"), numpy.zeros((4, 4), dtype="float32"))


def test_labels_accept_the_dtypes_cellprofiler_produces(tmp_path):
    """objects.segmented comes back as one of these depending on object count and platform."""
    for dtype in ("int16", "int32", "int64", "uint16", "uint32", "bool"):
        data = numpy.zeros((4, 4), dtype=dtype)
        path = write_labels(str(tmp_path / f"{dtype}.h5"), data)
        assert read_array(path).dtype == numpy.dtype(dtype)


def test_large_label_values_survive(tmp_path):
    """A plate-scale run can exceed 32767 objects in one field of view, so the label dtype has to
    carry through rather than being narrowed on write."""
    data = numpy.zeros((4, 4), dtype="int32")
    data[0, 0] = 70000
    path = write_labels(str(tmp_path / "a.h5"), data)
    assert read_array(path)[0, 0] == 70000
