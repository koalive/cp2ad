"""Folder layout and manifest for ExportForSpatialData.

One CellProfiler run writes one folder per plate, because one SpatialData object per plate is the
target and keeping the on-disk unit the same as the target unit is what keeps the importer simple:

    <prefix>_export/<plate>/
        images/<sample key>.h5              (C, Y, X), one stack per field of view
        labels/<sample key>/<object>.h5     (Y, X), integer labels, 0 for background
        tables/<prefix>.h5ad                one row per cell, plus the manifest in uns

The importer finds everything else from uns["cellprofiler_mapping"]["elements"], so it never walks
the filesystem or parses a file name. Manifest paths are relative to the plate folder, so the
folder can be moved or renamed without invalidating them.

Everything here is pure, with no CellProfiler imports, so it is all testable directly.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, fields
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy

from .assemble import Table
from .h5ad import Frame
from .samples import PLATE, SampleNaming, format_tag, sample_key

IMAGES, LABELS, TABLES = "images", "labels", "tables"
ELEMENT_IMAGE, ELEMENT_LABELS = "image", "labels"
STATUS_OK, STATUS_FAILED = "ok", "failed"

# One folder when the pipeline defines no plate tag: every image set is assumed to come from the
# same plate. That assumption is reported as a module warning and in the documentation, since it is
# wrong for a run spanning several plates and there is no way to tell from the metadata.
UNKNOWN_PLATE = "plate"

# Where run() records why a cycle failed, so post_run() can put the reason in the manifest instead
# of only reporting that a file is missing. Read back per image set from Image measurements, which
# is the one channel that carries values from a worker process to the main one.
ERROR_MEASUREMENT = "ExportForSpatialData_Error"

# obs columns SpatialData joins a table to its label arrays by: which element a row annotates, and
# which integer in that element is the row. label_id is already there for every exported object.
REGION_KEY, INSTANCE_KEY = "region_key", "label_id"

_UNSAFE = re.compile(r"[^A-Za-z0-9._-]+")


@dataclass(frozen=True)
class ElementRow:
    """One image stack or one label array on disk, as the importer needs to see it."""
    sample_key: str
    image_number: int
    element_type: str           # ELEMENT_IMAGE | ELEMENT_LABELS
    element_name: str           # the object name; empty on an image row
    path: str                   # relative to the plate folder
    shape: str                  # comma-joined, e.g. "4,2048,2048"
    # Not `dtype`: pandas resolves df.dtype to a column of that name, and anndata's writer reads
    # elem.dtype.kind to dispatch, so the column made the whole table unwritable to zarr.
    element_dtype: str
    region_key_value: str       # joins to obs["region_key"]; empty on an image row
    status: str                 # STATUS_OK | STATUS_FAILED
    error: str


@dataclass(frozen=True)
class ChannelAxisRow:
    """One position on an image stack's channel axis. Separate from the elements table because it
    is identical for every field of view in a run, and because a flat two-column table is less
    error-prone for the importer than a delimited string it would have to split."""
    channel: str
    stack_index: int


def safe_segment(text: str, fallback: str = "unnamed") -> str:
    """One path segment from a metadata value. Plate and well names come from acquisition software
    and can carry slashes or spaces, and a slash would silently create a nested folder."""
    cleaned = _UNSAFE.sub("_", str(text).strip()).strip("._-")
    return cleaned or fallback


def plate_of(md: Dict[str, Any], plate_tag: str = PLATE) -> str:
    """The plate folder name for one image set, or UNKNOWN_PLATE when the pipeline has no plate
    tag or the value is blank for this image set."""
    return safe_segment(format_tag(md.get(plate_tag)), UNKNOWN_PLATE)


def plates_by_image(values: Dict[int, Dict[str, Any]], plate_tag: str = PLATE) -> Dict[str, List[int]]:
    """{plate folder: image numbers}, in image-number order. A run covering two plates writes two
    folders, which is what makes one folder correspond to one SpatialData object."""
    out: Dict[str, List[int]] = {}
    for image_number in sorted(values):
        out.setdefault(plate_of(values[image_number], plate_tag), []).append(image_number)
    return out


def image_path(key: str) -> str:
    return f"{IMAGES}/{key}.h5"


def labels_path(key: str, obj: str) -> str:
    return f"{LABELS}/{key}/{safe_segment(obj)}.h5"


def table_path(prefix: str) -> str:
    return f"{TABLES}/{safe_segment(prefix, 'table')}.h5ad"


def region_key_value(key: str, obj: str) -> str:
    """The labels element's name, which is also what obs["region_key"] holds for its rows. One
    function so the element name and the column can never disagree."""
    return f"{key}__{obj}"


def spatialdata_attrs(regions: Iterable[str]) -> Dict[str, Any]:
    """uns["spatialdata_attrs"] for an exported table: which elements its rows annotate.

    The per-object table builder names the region after the object, `Cells`, because that is the
    only region an ExportToAnnData file has. Here the rows annotate one labels element per field of
    view, `<field>__Cells`, so the region is the list of those. Getting this wrong is not a detail
    SpatialData tolerates: TableModel.parse refuses a table whose region names no element it holds,
    which is [issue #414](https://github.com/scverse/spatialdata/issues/414).
    """
    return {"region": sorted(set(str(r) for r in regions)),
            "region_key": REGION_KEY, "instance_key": INSTANCE_KEY}


def channel_axis_rows(channels: Sequence[str]) -> List[ChannelAxisRow]:
    return [ChannelAxisRow(channel=c, stack_index=i) for i, c in enumerate(channels)]


def _row(sample_key_value: str, image_number: int, element_type: str, element_name: str,
         path: str, root: str, error: str) -> ElementRow:
    """One manifest row, with shape and dtype read back off the file that was actually written.

    Reading rather than remembering is deliberate: it confirms the file exists and is readable, so
    a cycle that failed or was never run shows up as failed without run() having to report
    anything in the ordinary case.
    """
    from .raster import array_info                      # local: keeps h5py off this module's import
    shape, element_dtype, status = "", "", STATUS_FAILED
    if not error:
        try:
            dims, element_dtype = array_info(os.path.join(root, path))
            shape, status = ",".join(str(d) for d in dims), STATUS_OK
        except Exception as exc:                        # missing, truncated, or not readable
            error = f"{type(exc).__name__}: {exc}"
    return ElementRow(
        sample_key=sample_key_value, image_number=image_number, element_type=element_type,
        element_name=element_name, path=path, shape=shape, element_dtype=element_dtype,
        region_key_value=region_key_value(sample_key_value, element_name) if element_name else "",
        status=status, error=error)


def element_rows(root: str, image_numbers: Iterable[int], values: Dict[int, Dict[str, Any]],
                 naming: SampleNaming, objects: Sequence[str],
                 errors: Optional[Dict[int, str]] = None) -> List[ElementRow]:
    """The manifest for one plate folder: one row per (field of view, element).

    `root` is the plate folder on disk, `values` the metadata per image set, `objects` the label
    arrays that were exported, and `errors` whatever run() recorded per image set. Every row's
    shape and dtype come from the file, so this doubles as a check that the export is complete.
    """
    errors = errors or {}
    rows: List[ElementRow] = []
    for image_number in image_numbers:
        md = values.get(image_number, {})
        key = sample_key(md, image_number, naming)
        error = errors.get(image_number, "")
        rows.append(_row(key, image_number, ELEMENT_IMAGE, "", image_path(key), root, error))
        for obj in objects:
            rows.append(_row(key, image_number, ELEMENT_LABELS, obj, labels_path(key, obj),
                             root, error))
    return rows


def manifest_to_uns(rows: Sequence[ElementRow],
                    channels: Sequence[ChannelAxisRow]) -> Dict[str, Frame]:
    """The two entries ExportForSpatialData adds to uns["cellprofiler_mapping"], in the same
    dataframe encoding as the channel, object and measurement tables already there, so they read
    back as pandas DataFrames the importer can filter."""
    return {
        "elements": Frame.from_records([_asdict(r) for r in rows],
                                       [f.name for f in fields(ElementRow)]),
        "image_channels": Frame.from_records([_asdict(r) for r in channels],
                                             [f.name for f in fields(ChannelAxisRow)]),
    }


def _asdict(row) -> Dict[str, Any]:
    return {f.name: getattr(row, f.name) for f in fields(type(row))}


def region_key_column(table: Table, values: Dict[int, Dict[str, Any]], naming: SampleNaming,
                      base_object: str) -> numpy.ndarray:
    """obs["region_key"] for a joined table: the labels element each row belongs to.

    `values` is the per-image-set metadata, the same dict element_rows() names files from, so the
    column and the manifest cannot disagree. Reading the tags out of table.obs instead would look
    equivalent and is not: build_object_table drops a Metadata column that is missing for every
    image set, so a key built from obs could come out a component short while the manifest, built
    from measurements, kept it.
    """
    keys = []
    for image_number in table.obs["ImageNumber"]:
        image_number = int(image_number)
        key = sample_key(values.get(image_number, {}), image_number, naming)
        keys.append(region_key_value(key, base_object))
    return numpy.array(keys, dtype=object)


def subset_table(table: Table, keep: numpy.ndarray) -> Table:
    """The rows `keep` selects, as a new Table. Splitting one run's joined table into one table
    per plate, so each plate folder is self-contained.

    var is shared unchanged. uns is shared too, except qc_summary, which is recounted: a run-level
    count copied into every plate would overstate each of them.
    """
    keep = numpy.asarray(keep, dtype=bool)
    uns = dict(table.uns)
    if "qc_flag" in table.obs and "qc_summary" in uns:
        flags = numpy.asarray(table.obs["qc_flag"])[keep]
        summary = {str(k): int(v) for k, v in zip(*numpy.unique(flags, return_counts=True))}
        # primaries_without_secondary is a run-level count that cannot be split per plate, so it
        # is carried over rather than guessed at.
        if "primaries_without_secondary" in uns["qc_summary"]:
            summary["primaries_without_secondary"] = int(
                uns["qc_summary"]["primaries_without_secondary"])
        uns["qc_summary"] = summary
    return Table(
        X=table.X[keep],
        obs_names=[name for name, take in zip(table.obs_names, keep) if take],
        var_names=list(table.var_names),
        obs={k: numpy.asarray(v)[keep] for k, v in table.obs.items()},
        var=dict(table.var),
        obsm={k: numpy.asarray(v)[keep] for k, v in table.obsm.items()},
        uns=uns)


def selected_channels(ctx, requested: Sequence[str]) -> List[str]:
    """The channels to export: exactly those requested, or every channel loaded from a file when
    nothing was requested. Derived images stay out of the default because a pipeline can produce
    many of them and a folder of intermediates is not what anyone asked for."""
    if requested:
        return [c for c in ctx.channels if c in set(requested)]
    return [c for c in ctx.channels
            if getattr(ctx.channel_info.get(c), "source", None) == "file"]


def selected_objects(ctx, requested: Sequence[str]) -> List[str]:
    """The label arrays to export: exactly those requested, or every object the pipeline made.

    Every object by default, unlike the per-cell table, which is scoped to the primary, secondary
    and tertiary roles because joining compartments into one row is what it is for. A label array
    has no such constraint, and a pipeline that segments spots outside the role chain still wants
    them visible in the viewer.
    """
    if requested:
        return [o for o in ctx.objects if o in set(requested)]
    return list(ctx.objects)
