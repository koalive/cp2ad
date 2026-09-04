"""Image and label arrays on disk, one array per file, written with h5py alone.

CellProfiler ships numpy and h5py and nothing else this package can rely on, so the arrays that
will become a SpatialData store are plain HDF5 rather than zarr. Each file holds one dataset named
"data": an image stack shaped (C, Y, X), or a label array shaped (Y, X) whose integers are the
object labels CellProfiler assigned, 0 for background. The importer reads these and hands them to
Image2DModel or Labels2DModel. Nothing here depends on which library writes the bytes, so a future
CellProfiler that ships zarr only changes this module.

Writes land on a temporary name in the destination directory and are renamed into place. Rename
within one filesystem is atomic on POSIX, so a crash or a killed worker mid-write leaves either the
old file or the new one, never a truncated file that the importer would read as valid.
"""
from __future__ import annotations

import os
import tempfile
from typing import Optional, Tuple

import h5py
import numpy

DATASET = "data"

# Lossless, and level 1 is the measured sweet spot. On four real BBBC021 channels (10.5 MB raw)
# it saves 51% in 0.14 s, where level 4 saves 53% in 0.32 s. Label arrays are mostly runs of
# zeros and go from 16.8 MB to 0.2 MB in 0.02 s. Worth having on by default, because an
# uncompressed 4-channel 2048x2048 uint16 stack is 33.6 MB and a 384-well plate at 4 sites per
# well runs to about 50 GB of images before any label arrays.
DEFAULT_COMPRESSION = "gzip"
DEFAULT_COMPRESSION_OPTS = 1


def write_array(path: str, data: numpy.ndarray,
                compression: Optional[str] = DEFAULT_COMPRESSION,
                compression_opts: Optional[int] = DEFAULT_COMPRESSION_OPTS) -> str:
    """Write one array to `path` as dataset "data", creating parent directories. Returns `path`.

    Pass compression=None to trade disk for write speed. h5py chooses a chunk shape itself when
    compression is on, which it has to be for a chunked dataset.
    """
    data = numpy.asarray(data)
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    # Same directory as the destination, so the rename stays within one filesystem.
    fd, tmp = tempfile.mkstemp(prefix=".tmp-", suffix=".h5", dir=directory)
    os.close(fd)
    try:
        with h5py.File(tmp, "w") as f:
            if data.size and compression:
                f.create_dataset(DATASET, data=data, compression=compression,
                                 compression_opts=compression_opts)
            else:
                # h5py cannot chunk a zero-size dataset, and chunking is required to compress one.
                f.create_dataset(DATASET, data=data)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
    return path


def read_array(path: str) -> numpy.ndarray:
    """The array back, as numpy. Used by the tests and by anything checking an export by hand;
    the real importer reads these files itself."""
    with h5py.File(path, "r") as f:
        return f[DATASET][()]


def array_info(path: str) -> Tuple[Tuple[int, ...], str]:
    """(shape, dtype name) without reading the pixels, for the manifest rows that record them."""
    with h5py.File(path, "r") as f:
        dset = f[DATASET]
        return tuple(dset.shape), dset.dtype.name


def write_image(path: str, data: numpy.ndarray, **kwargs) -> str:
    """One field of view's channel stack, shaped (C, Y, X). A single-channel image still gets a
    length-1 channel axis, because Image2DModel.parse() expects the axis to exist and guessing
    later from ndim alone would be ambiguous against a 2D image."""
    data = numpy.asarray(data)
    if data.ndim != 3:
        raise ValueError(f"an image stack must be (C, Y, X), got shape {data.shape}. "
                         "Add a length-1 channel axis for a single-channel image.")
    return write_array(path, data, **kwargs)


def write_labels(path: str, data: numpy.ndarray, **kwargs) -> str:
    """One object type's label array for one field of view, shaped (Y, X), integer labels with 0
    for background. This is CellProfiler's `objects.segmented` unchanged."""
    data = numpy.asarray(data)
    if data.ndim != 2:
        raise ValueError(f"a label array must be (Y, X), got shape {data.shape}")
    if data.dtype.kind not in "iub":
        raise ValueError(f"label arrays must hold integers, got dtype {data.dtype}. "
                         "CellProfiler's objects.segmented is already integer-labelled.")
    return write_array(path, data, **kwargs)
