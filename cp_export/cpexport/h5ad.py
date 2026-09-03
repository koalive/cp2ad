"""Minimal AnnData .h5ad writer using only h5py (the CellProfiler app has no pandas/anndata).
Follows the anndata on-disk spec: https://anndata.readthedocs.io/en/latest/fileformat-prose.html
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List

import h5py
import numpy

from .assemble import Table

STR = h5py.string_dtype(encoding="utf-8")


@dataclass
class Frame:
    """A tabular value destined for uns, written in anndata's on-disk dataframe encoding so it
    reads back as a pandas DataFrame rather than a nested dict.

    Without this, _write_elem stores a list of row-dicts as {"0": {...}, "1": {...}}, which is
    technically lossless but useless for anything tabular: you cannot slice it, sort it, or hand it
    to pandas without rebuilding it by hand. `columns` is ordered, so column order survives the
    round trip; `index` defaults to the row positions as strings.
    """
    columns: Dict[str, List[Any]] = field(default_factory=dict)
    index: List[str] = field(default_factory=list)

    @classmethod
    def from_records(cls, records: List[Dict[str, Any]], fieldnames: List[str]) -> "Frame":
        """Column-orient a list of same-shaped row dicts. `fieldnames` fixes the column order, so
        an empty record list still produces a frame with the right (empty) columns."""
        columns = {name: [r[name] for r in records] for name in fieldnames}
        return cls(columns=columns, index=[str(i) for i in range(len(records))])


def _enc(g, etype: str, version: str):
    g.attrs["encoding-type"] = etype
    g.attrs["encoding-version"] = version


def _text(v) -> str:
    if isinstance(v, (bytes, numpy.bytes_)):
        return v.decode("utf-8")
    return "" if v is None else str(v)


def _write_string_array(g, key: str, values):
    arr = numpy.array([_text(v) for v in values], dtype=object)
    d = g.create_dataset(key, data=arr, dtype=STR)
    _enc(d, "string-array", "0.2.0")


def _write_pairs(g, key: str, pairs):
    """(text, value) pairs -- module setting_values -- as an n x 2 string array, so they read back
    as pairs instead of as the repr of a tuple."""
    arr = numpy.array([[_text(a), _text(b)] for a, b in pairs], dtype=object)
    d = g.create_dataset(key, data=arr, dtype=STR)
    _enc(d, "string-array", "0.2.0")


def _write_array(g, key: str, arr: numpy.ndarray):
    d = g.create_dataset(key, data=arr)
    _enc(d, "array", "0.2.0")


def _write_column(g, key: str, arr: numpy.ndarray):
    arr = numpy.asarray(arr)
    if arr.dtype.kind in ("O", "U", "S"):
        _write_string_array(g, key, arr)
    elif arr.dtype.kind == "b":
        _write_array(g, key, arr.astype(bool))
    else:
        _write_array(g, key, arr)


def _write_dataframe(parent, key: str, index, columns: Dict[str, numpy.ndarray]):
    g = parent.create_group(key)
    _enc(g, "dataframe", "0.2.0")
    g.attrs["_index"] = "_index"
    g.attrs["column-order"] = numpy.array(list(columns), dtype=STR)
    _write_string_array(g, "_index", index)
    for name, arr in columns.items():
        _write_column(g, name, arr)


def _write_elem(g, key: str, value: Any):
    if value is None:
        d = g.create_dataset(key, data="", dtype=STR)  # anndata has no null; "" is the least surprising stand-in
        _enc(d, "string", "0.2.0")
    elif isinstance(value, Frame):
        _write_dataframe(g, key, value.index, {k: numpy.asarray(v) for k, v in value.columns.items()})
    elif isinstance(value, dict):
        sub = g.create_group(key)
        _enc(sub, "dict", "0.1.0")
        for k, v in value.items():
            # '/' is the HDF5 path separator; keep setting texts as flat keys
            _write_elem(sub, str(k).replace("/", "|"), v)
    elif isinstance(value, str):
        d = g.create_dataset(key, data=value, dtype=STR)
        _enc(d, "string", "0.2.0")
    elif isinstance(value, (bool, numpy.bool_)):
        d = g.create_dataset(key, data=bool(value)); _enc(d, "numeric-scalar", "0.2.0")
    elif isinstance(value, (int, float, numpy.integer, numpy.floating)):
        d = g.create_dataset(key, data=value); _enc(d, "numeric-scalar", "0.2.0")
    elif isinstance(value, numpy.ndarray):
        _write_column(g, key, value)
    elif isinstance(value, (list, tuple)):
        if len(value) and all(isinstance(v, dict) for v in value):
            _write_elem(g, key, {str(i): v for i, v in enumerate(value)})
        elif len(value) and all(isinstance(v, (tuple, list)) and len(v) == 2 for v in value):
            _write_pairs(g, key, value)
        else:
            nums = [v for v in value if v is not None]
            if nums and all(isinstance(v, (int, float, numpy.integer, numpy.floating)) and not isinstance(v, bool) for v in nums):
                _write_array(g, key, numpy.array([numpy.nan if v is None else v for v in value], dtype=float)
                             if any(v is None or isinstance(v, float) for v in value)
                             else numpy.asarray(value, dtype=int))
            else:
                _write_string_array(g, key, value)
    else:
        _write_elem(g, key, str(value))


def write_h5ad(table: Table, path: str) -> None:
    with h5py.File(path, "w") as f:
        _enc(f, "anndata", "0.1.0")
        _write_array(f, "X", numpy.asarray(table.X, dtype=numpy.float32))
        _write_dataframe(f, "obs", table.obs_names, table.obs)
        _write_dataframe(f, "var", table.var_names, table.var)
        obsm = f.create_group("obsm"); _enc(obsm, "dict", "0.1.0")
        for k, v in table.obsm.items():
            _write_array(obsm, k, numpy.asarray(v))
        for empty in ("varm", "obsp", "varp", "layers"):
            _enc(f.create_group(empty), "dict", "0.1.0")
        uns = f.create_group("uns"); _enc(uns, "dict", "0.1.0")
        for k, v in table.uns.items():
            _write_elem(uns, k, v)
