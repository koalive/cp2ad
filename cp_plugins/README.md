# CellProfiler plugins for the scverse ecosystem

Two modules, sharing one export backend:

| Module | Writes | Read it with |
|---|---|---|
| **ExportToAnnData** | one `.h5ad` per run, one row per cell | anndata, scanpy, squidpy |
| **ExportForSpatialData** | one folder per plate: image stacks, label arrays, and the same table | a SpatialData importer (separate) |

Feature names, `obs`/`uns` keys and layout match
`squidpy.experimental.im.calculate_image_features` (cp_measure naming), so a CellProfiler run and a
squidpy run on the same masks are interchangeable per compartment. Both modules build their table
the same way, so [Object roles](#object-roles), [What lands where](#what-lands-where) and
[Row names](#row-names) describe both.

ExportForSpatialData's design notes, including what is deferred to a second phase:
[`../docs/2026-09-03-exportforspatialdata-design.md`](../docs/2026-09-03-exportforspatialdata-design.md).

## Install
Point CellProfiler at this folder: `Preferences -> CellProfiler plugins directory`, or headless
`--plugins-directory=/path/to/cp_plugins`. Both modules appear under *File Processing*.

```sh
/Applications/CellProfiler.app/Contents/MacOS/cp -c -r -p pipeline.cppipe -i images -o out --plugins-directory=cp_plugins
```

Add either module at the end of the pipeline. A module after it that makes features may not have
them exported, which raises a warning rather than an error: another exporter or SaveImages after it
is normal. Both modules in one pipeline is normal too.

## ExportToAnnData
Writes `<prefix>.h5ad` (joined, `var_names = <Object>__<feature>`), optionally also
`<prefix>_<Object>.h5ad` per compartment. Defaults work for IdentifyPrimary -> Secondary ->
Tertiary pipelines. Nothing happens during the run: the whole object is assembled and written in
`post_run`.

## ExportForSpatialData
Exports the pixels alongside the table, which is what ExportToAnnData cannot do: image stacks and
segmentation masks on disk, so a SpatialData object can be assembled with images, labels and the
per-cell table in one coordinate system.

One folder per plate, because one SpatialData object per plate is the target:

```
<prefix>_export/
    <plate>/
        images/<sample key>.h5              (C, Y, X) float32, one stack per field of view
        labels/<sample key>/<object>.h5     (Y, X) integer labels, 0 for background
        tables/<prefix>.h5ad                that plate's rows, manifest in uns
```

Image stacks and label arrays are written cycle by cycle, so a long run does not hold them in
memory. The tables are written in `post_run`, one per plate, each a subset of the run's rows with
`uns["qc_summary"]` recounted for that plate.

An importer reads `uns["cellprofiler_mapping"]["elements"]` and never walks the folder or parses a
file name. One row per (field of view, element), with `path` relative to the plate folder, so the
folder can be moved or renamed. Each row carries the `shape` and `dtype` read back off the file
that was written, and `status` is `failed` with the reason in `error` when a file is missing or
unreadable. `image_channels` gives the channel axis as (`channel`, `stack_index`) pairs.
`obs["region_key"]` on the table joins to the `region_key_value` of a labels row.

A cycle that fails to write does not stop the run: the error is logged, recorded as an image
measurement, and reported in the manifest.

*Images to export* and *Segmentations to export* are prefilled on first open, with every channel
loaded from a file and every object the pipeline made. Narrow them as you like; a narrowed
selection stays narrowed.

Read a plate folder back with `cell_painting_io.read_cellprofiler_export`, in
[cell-painting-io](https://github.com/scverse/cell-painting-io), which turns one folder into a
SpatialData object: an Image and a Labels element per field of view, and the table annotating them.
`cellprofiler_export_plates` lists the plate folders of a multi-plate run. Each field gets its own
coordinate system and the fields are not placed relative to each other, since stage coordinates are
phase 2.

**With no `Metadata_Plate` in the pipeline, every image set is assumed to come from one plate** and
lands in a single folder named `plate`. That is wrong for a run spanning several plates, and there
is no way to tell from the metadata, so the module reports it as a warning and logs it at the start
of the run. Define a plate tag in Metadata if the run covers more than one plate.

## Object roles
Objects carry one of three roles: `primary`, `secondary`, `tertiary`. In a standard pipeline
primary = nucleus, secondary = cell, tertiary = cytoplasm, but the names stay generic because
nothing in the exporter assumes that biology. The joined row is one secondary object (or, when
there is none, one primary object) with its primary and tertiary object joined onto it.

*How to pick primary / secondary / tertiary objects* = **Automatic** (default) reads them off the
pipeline: IdentifyTertiaryObjects names all three (smaller = primary, larger = secondary, output =
tertiary); else IdentifySecondaryObjects names primary (its input) and secondary (its output); else
the primary candidates.

Ambiguity at any of those three levels is resolved, never refused, by the same score: how many
*other* modules name the object in a setting, i.e. use it as an input (RelateObjects, FilterObjects,
MaskObjects, MeasureObject* selections, whose multi-object values are split on `", "`). Highest score
wins, ties go to the object produced first.

- several IdentifyTertiaryObjects -> the most-related tertiary output object; primary and secondary
  come from that same module's smaller/larger settings, so the chosen chain keeps all three
  compartments (`fallback: most_related_tertiary`).
- no tertiary but several IdentifySecondaryObjects -> the most-related secondary output object,
  primary = its input (`fallback: most_related_secondary`).
- several primary candidates and no chain at all -> the most-related one
  (`fallback: most_related_primary`).

The choice is logged as a warning, repeated as advice and recorded in
`uns["cellprofiler"]["role_detection"]` (`{mode: automatic|manual, fallback: most_related_tertiary |
most_related_secondary | most_related_primary | null, candidates: [...]}`). The run is only refused
when the pipeline produces no objects at all.

**Manual** shows *Primary objects (e.g. nuclei)*, *Secondary objects (e.g. cells)* and *Tertiary
objects (e.g. cytoplasm)*. Manual mode exports exactly the roles you set. A role left at `None` is
not exported, and nothing is auto-detected alongside your choice, so naming only the primary object
gives a one-compartment export of that object.

Roles scope the per-cell table only. ExportForSpatialData exports label arrays for every object by
default, roles included or not, since a pipeline that segments spots outside the role chain still
wants them in the viewer.

A `FilterObjects` module is followed as a *rename* only when it relabels several objects together
(the JUMP border-cleanup step that relabels nucleus and cell in one go). A single-pair
`FilterObjects` makes a filtered subset: its output inherits its input's role and is scored as an
ordinary candidate, but it never stands in for its input. So a pipeline that measures
`FilteredNuclei` exports `FilteredNuclei`, while one that filters `Nuclei` into `PH3PosNuclei` and
keeps measuring `Nuclei` exports `Nuclei` rather than silently exporting only the positive nuclei.

## What lands where
- `X`: float32, NaN kept, no constant-column drop. Only numeric CellProfiler columns (`float`,
  `integer`) become features, and only the ones that describe the object's own morphology,
  intensity or texture. `varchar`/`blob` measurements are left out, and so is anything extrinsic to
  the object's own biology (see below).
- `var`: `cp_name`, `module_num`, `module_name`, `category`, `measurement`, `channel`, `channel2`,
  `other_object`, `scale`, `coltype`, `parsed_by` (`api`|`fallback`), `region`.
- **Extrinsic measurements never reach `X`/`var`.** Two objects with identical biology should land
  at the same point in feature space no matter where they were imaged, how the sample was rotated,
  or what label CellProfiler happened to assign them or their neighbors. These go to `obs` instead:
  `cp_measure`-named columns (`Center_X`, `Orientation`, `Parent_Nuclei`, ...) on the per-object
  files, `<Object>__<name>` (e.g. `Nuclei__Center_X`, `Cells__Number_Object_Number`) on the joined
  file. `is_extrinsic()` in `scverse_export/names.py` is the single source of truth for the split; it
  currently covers:
  - **Position and orientation**: the whole `Location` category (an object's own center and, per
    channel, its intensity-weighted/max-intensity center, all absolute pixel coordinates) plus the
    `AreaShape` measurements that are themselves a coordinate or an angle in the image frame
    (`Center_X`/`Y`, the four `BoundingBox{Minimum,Maximum}_{X,Y}` corners, `Orientation`).
  - **Identity and linkage**: `Number_Object_Number` (the object's own arbitrary label, already
    carried as `obs["label_id"]`), `Parent_<Object>` (a label referencing another row, not a
    measurement of this one), `Children_<Object>_Count` (the joined file already surfaces this as
    `count_<child>` for non-role children, see below, so it stays out of `X` rather than duplicating
    that value under a second, inconsistent name), and `Neighbors_{First,Second}ClosestObjectNumber`
    (another object's label, the same problem as `Parent_*`).
  - Deliberately **kept in `X`**: `Neighbors_NumberOfNeighbors`, `Neighbors_PercentTouching`,
    `Neighbors_{First,Second}ClosestDistance` and `Neighbors_AngleBetweenNeighbors`. Local crowding,
    density, and the alignment between neighboring cells count as biology here, unlike an object's
    own absolute position.
  - Press **"See where each measurement will land"** on ExportToAnnData to see this split
    applied to your own pipeline, before running it. Three tables cover the per-cell columns
    (channels measured, objects, measurements), each naming the module that produced the entry and,
    for measurements, the exact resulting name. They list only what becomes a per-cell column, so
    an absence is informative: a channel loaded but never measured (a raw image only used to build
    a derived one, say) gets no row. A fourth table, **Also exported, in uns**, accounts for
    everything else the file carries, so an absence upstream never has to be read as data loss.
    Turn on **"Show advanced features?"** for **"Also write the mapping tables to .uns?"**, which
    saves those first three tables into `uns["cellprofiler_mapping"]` on every run as three
    DataFrames.
- `obs`: `region`, `label_id`, `ImageNumber`, `Metadata_*`, `n_missing_features`, the extrinsic
  columns above, plus `qc_flag` and `count_<child>` **on the joined file only** (the per-object
  files have neither). ExportForSpatialData adds `region_key`.
- `obs_names`: `<sample key>_<label>`, where the sample key identifies the field of view and the
  label is the integer CellProfiler gave the object there. See "Row names" below.
- `obsm["spatial"]`: cell centers in site pixel coordinates, the base object's `Location_Center_X/Y`.
  This duplicates `obs["Center_X"/"Center_Y"]` / `obs["<base object>__Center_X"/"Center_Y"]` on
  purpose: `obsm["spatial"]` is for spatial-analysis tooling that expects it there, e.g. squidpy.
- `uns["spatialdata_attrs"]`: region/instance keys (drops into SpatialData as a TableModel).
- `uns["qc_summary"]` counts are computed before any `Drop` filtering. `qc_flag` values: `ok`,
  `no_primary`, `multi_secondary_per_primary`, `no_tertiary`, `multi_tertiary`; plus the
  summary-only count `primaries_without_secondary`.
- `uns["cellprofiler_join"]` (joined file), keyed by role (`primary`, `tertiary`): how each
  compartment was matched to the base row. Either the `Parent_<Object>` column, or
  `shared_label_id` when one FilterObjects module relabelled both objects (the JUMP chain writes no
  `Parent_Nuclei`). If neither is available the run stops with a `JoinError`.
- setting texts stored as `uns` keys have `/` replaced by `|` (HDF5 path separator); a `None` value
  is written as a string scalar `""` (`.h5ad` has no null).
- `uns["cellprofiler"]`: CP version, pipeline text, modules (`setting_values` keeps repeated setting
  texts, e.g. one row per texture scale; the collapsed first-occurrence-only `settings` dict is
  internal to the exporter and is not exported), channels, objects (with `source`: `pipeline`|
  `file`), roles, role_detection, relationships, all image-level measurements (thresholds, counts,
  QC, timings), experiment measurements and exporter settings.
- `uns["cellprofiler"]["image"]` deserves calling out, because it is where everything that is not a
  per-cell measurement goes: one entry per image set, holding `FileName_*`, `PathName_*` and `URL_*`
  per channel, `Metadata_*`, `Count_*`, `Threshold_*`, `ExecutionTime_*` and the rest. File names and
  paths are therefore exported, just keyed by image set rather than repeated on every cell. It is
  column-oriented, so `pd.DataFrame(adata.uns["cellprofiler"]["image"])` gives the table directly.
  `uns["cellprofiler"]["channels"]` likewise names **every** channel, including ones no measurement
  reads.
- `uns["cellprofiler_mapping"]`: the preview's `channels`, `objects` and `measurements` tables,
  written in AnnData's dataframe encoding, so each reads back as a `pd.DataFrame` ready to filter
  and sort. Optional on ExportToAnnData (*Also write the mapping tables to .uns?*); always written
  by ExportForSpatialData, which adds `elements` and `image_channels` to it.

## Row names

Each row is named `<sample key>_<label>`: the field of view, then the integer CellProfiler gave the
object in it. The sample key is an advanced setting, called *How to build row names* on
ExportToAnnData and *How to name fields of view* on ExportForSpatialData, where the same key also
names the image stacks, the label arrays and the coordinate systems.

**Automatic** (default) reads the key off the pipeline's Metadata tags, taking at most one plate
part (`Plate`, `Barcode`, `PlateID`), one well part (`Well`, or `Row` plus `Column`), and one site
part (`Site`, `Field`, `FieldIndex`, `Position`). Plate/Well/Site gives `P1_A01_1_5`; Opera's `Well`
and `Field` give `A02_03_5`. `Frame`, `Series` and `Channel` are never used: they index a z plane, a
timepoint or a channel inside one field of view, not the field.

The run checks the chosen tags give every image set a different key, and appends `img<n>` when they
do not, so two fields of view never share a row name. With no usable tags the key is `img<n>`
alone, unique within a run but not across runs.

**Manual** uses the comma-separated tags you list, where the `Metadata_` prefix is optional and an
empty list means image number alone. A tag the pipeline lacks is reported, and the key falls back
to `img<n>` rather than leaving a blank mid-name.

Manual exists because Automatic depends on what a run contains: the same pipeline over a different
subset can name rows differently. Pin the tags before exporting anything you mean to concatenate.

**"Auto-configure from this pipeline"** (also advanced) fills in the object roles and the row-name
tags from the current pipeline and switches both to Manual. It cannot check uniqueness, which needs
real values, so the run still appends `img<n>` if the pinned tags turn out not to separate the
image sets.

The resolved scheme is logged and recorded in `uns["cellprofiler"]["sample_naming"]`.
ExportToAnnData also shows it at the top of "See where each measurement will land".

ExportForSpatialData **always** appends `img<n>`, giving `A02_03_img1_5`. It names files during the
run and table rows afterwards, possibly from different processes, and the two have to agree exactly
or a manifest row points at a file written under another name. The uniqueness check needs every
image set at once, which `run()` cannot see, so the image number makes keys unique by construction
instead. The cost is a longer key when the tags alone would have done.

## Limits
- **Image stacks are written as CellProfiler's `pixel_data`**, which is float32 rescaled to [0, 1],
  not the source `uint16`. Values are the ones the pipeline measured, so the export is faithful to
  the run, but it is not the raw file and it is 2x the size of `uint16`. Recovering the original
  range needs the source images.
- **Label dtype varies between fields of view.** CellProfiler narrows `objects.segmented` to fit
  the object count, so a field with 102 objects gives `int8` and a busier one `int16`. Every array
  is correct on its own, and the manifest reports each dtype, but an importer stacking them has to
  promote rather than assume one type.
- **Nothing is streamed on the table side.** See below.

## Scaling
Sized for per-well or per-site-batch runs. Both modules assemble the table in `post_run`, so every
cell of every image set in one CellProfiler run is in memory at once. ExportForSpatialData streams
its image stacks and label arrays per cycle, but not its table.

`build_object_table` fills `X` feature-major. One `get_measurement` call reads one feature's
float64 series across every image set, writes it into `X`'s float32 rows, then drops it before the
next feature. Peak memory is therefore `X` plus one feature's float64 series, not a float64 copy of
every feature at once.

A whole plate in a single run will still exhaust memory: on the order of 5M cells x 8k features is
~160 GB as float32 for `X` alone. Split the plate into per-well or per-site runs and concatenate
afterwards. `obs_names` stay unique across runs when the sample key includes a plate and a site tag,
since those identify a field of view independently of the run that processed it. A key that falls
back to `img<n>` does not, because image numbers restart at 1 every run. The log line before each
write reports the shape and size.

## Tests
Most of the suite needs only numpy and h5py:
```sh
cd cp_plugins && python -m pytest -q
```
The tests that import `cellprofiler_core` are skipped without it. To run those too:
```sh
cd plugin_sandbox && pixi install && pixi run pip install "numpy<2" "cython<3" setuptools wheel \
  && pixi run pip install --no-build-isolation python-javabridge==4.0.5 \
  && pixi run pip install "cellprofiler-core==4.2.8.1" pytest "anndata<0.11" pandas && cd ..
cd cp_plugins && ../plugin_sandbox/.pixi/envs/default/bin/python -m pytest -q
```
`tests/test_integration_app.py` runs the real app binary headless and is skipped when it is absent.
