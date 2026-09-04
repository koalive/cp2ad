# ExportToAnnData — CellProfiler plugin

Writes one `.h5ad` per run with one row per cell. Feature names, `obs`/`uns` keys and layout match
`squidpy.experimental.im.calculate_image_features` (cp_measure naming), so a CellProfiler run and a
squidpy run on the same masks are interchangeable per compartment. Design: `../docs/2026-08-29-exporttoanndata-design.md`.

## Use
1. Point CellProfiler at this folder: `Preferences -> CellProfiler plugins directory`, or headless
   `--plugins-directory=/path/to/cp_plugins`.
2. Add **ExportToAnnData** as the last module. Defaults work for IdentifyPrimary -> Secondary -> Tertiary pipelines.
3. Output: `<prefix>.h5ad` (joined, `var_names = <Object>__<feature>`), optionally `<prefix>_<Object>.h5ad`.

## Object roles
Objects carry one of three roles -- `primary`, `secondary`, `tertiary`. In a standard pipeline
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
objects (e.g. cytoplasm)*. Manual mode exports exactly the roles you set -- a role left at `None` is
not exported and nothing is auto-detected alongside your choice, so naming only the primary object
gives a one-compartment export of that object.

A `FilterObjects` module is followed as a *rename* only when it relabels several objects together
(the JUMP border-cleanup step that relabels nucleus and cell in one go). A single-pair
`FilterObjects` makes a filtered subset: its output inherits its input's role and is scored as an
ordinary candidate, but it never stands in for its input. So a pipeline that measures
`FilteredNuclei` exports `FilteredNuclei`, while one that filters `Nuclei` into `PH3PosNuclei` and
keeps measuring `Nuclei` exports `Nuclei` rather than silently exporting only the positive nuclei.

```sh
/Applications/CellProfiler.app/Contents/MacOS/cp -c -r -p pipeline.cppipe -i images -o out --plugins-directory=cp_plugins
```

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
  - Press **"See where each measurement will land"** on the module itself to see this split
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
  files have neither). `obs_names`: `<Plate>_<Well>_<Site>_<label>` when Metadata Plate/Well/Site
  exist, else `img<n>_<label>`.
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
- `uns["cellprofiler_mapping"]` (only with *Also write the mapping tables to .uns?*): the preview's
  `channels`, `objects` and `measurements` tables, written in AnnData's dataframe encoding, so each
  reads back as a `pd.DataFrame` ready to filter and sort.

## Scaling
Sized for **per-well or per-site-batch CellProfiler runs**. Nothing is streamed: `run()` is a no-op and the
whole object is assembled in `post_run`, so every cell of every image set in one CellProfiler run is held in
memory at once. `build_object_table` fills `X` feature-major -- one `get_measurement` call reads one
feature's float64 series across every image set, writes it into `X`'s float32 rows, then drops it before
the next feature -- so peak memory is `X` (float32, the whole object) plus one feature's float64 series,
not every feature's float64 copy at once. A whole plate in a single run -- of the order of 5M cells x 8k
features, ~160 GB as float32 for `X` alone -- will still exhaust memory. Split the plate into per-well or
per-site runs and concatenate the `.h5ad` files afterwards; `obs_names` stay unique across runs as long as
Metadata Plate/Well/Site are set. The log line before each write reports the shape and size.

## Tests
```sh
cd plugin_sandbox && pixi install && pixi run pip install "numpy<2" "cython<3" setuptools wheel \
  && pixi run pip install --no-build-isolation python-javabridge==4.0.5 \
  && pixi run pip install "cellprofiler-core==4.2.8.1" pytest "anndata<0.11" pandas && cd ..
cd cp_plugins && ../plugin_sandbox/.pixi/envs/default/bin/python -m pytest -q
```
`tests/test_integration_app.py` runs the real app binary headless and is skipped when it is absent.
