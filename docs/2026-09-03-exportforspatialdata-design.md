# ExportForSpatialData: design plan

Status: planning only. Nothing here is implemented or committed. Scope is the CellProfiler module. The SpatialData importer is a separate piece of work, covered here only as the contract the module has to satisfy.

## Goal

ExportToAnnData exports numbers. ExportForSpatialData exports the pixels too: images and segmentation labels alongside the per-cell table, so the result becomes a real SpatialData object rather than a feature matrix. Cell 47 stops being a row of numbers and becomes a labelled region in a specific image, which you can overlay, crop, and re-measure.

CellProfiler cannot write a `.zarr` store, so the module writes HDF5 arrays in a folder tree shaped like the store the importer will build.

One SpatialData object per plate is the requirement, so the module writes one folder per plate. The on-disk unit and the target object are the same unit, which is what keeps the importer simple.

## What the module writes

```
<prefix>_export/
└── P1/                            one plate, one SpatialData object
    ├── images/                    the pixels: one stack per field of view
    │   ├── P1_A01_1.h5            dataset "data", (C, Y, X), source dtype
    │   └── P1_A01_2.h5
    ├── labels/                    the segmentations: one array per object per FOV
    │   ├── P1_A01_1/
    │   │   ├── Nuclei.h5          dataset "data", (Y, X), integer labels, 0 = background
    │   │   ├── Cells.h5
    │   │   └── Spots.h5           non-role objects included by default
    │   └── P1_A01_2/
    │       └── ...
    └── tables/
        └── <prefix>.h5ad          one row per cell, every FOV this run covered on P1
```

`sample_key` is `<Plate>_<Well>_<Site>` from the Metadata tags, the same tags `assemble._obs_name` already uses for `obs_names`.

When the pipeline defines no `Metadata_Plate`, the module assumes every image set comes from one plate and writes a single folder. That is the right assumption for the common case, but it is an assumption, so it surfaces twice: as a module warning through the existing `advice.py` mechanism, which `validate_module_warnings` raises in the GUI and the run logs, and in the module documentation. Worth stating the consequence too, since element names then carry no plate: two such exports are still each a valid object, but concatenating them needs distinct prefixes or the importer's suffix dictionary to keep names apart. Adding a Metadata module that extracts Plate is the real fix, and the warning should say so.

How one plate folder becomes one object:

```
P1/ on disk                              P1.zarr, one SpatialData object
────────────────────────────────────────────────────────────────────────────────
images/P1_A01_1.h5        (C,Y,X)  ──▶  images["P1_A01_1"]           Image2DModel
labels/P1_A01_1/Nuclei.h5   (Y,X)  ──▶  labels["P1_A01_1__Nuclei"]   Labels2DModel
labels/P1_A01_1/Cells.h5    (Y,X)  ──▶  labels["P1_A01_1__Cells"]    Labels2DModel
tables/<prefix>.h5ad               ──▶  tables["table"]              TableModel
                                        annotates the labels having rows,
                                        region_key -> label_id
every element of one FOV           ──▶  coordinate system "P1_A01_1"
```

Images are one stacked `(C, Y, X)` array per FOV because `Image2DModel.parse()` wants that shape, and `spatialdata_io`'s generic reader uses the same `cyx` convention. One array per channel would only force the importer to re-stack. Label arrays are one per object type because `Labels2DModel.parse()` wants a single integer-labelled raster per call, which is exactly what `objects.segmented` already is. No `shapes/` or `points/`, matching your sketch.

A run that spans two plates writes two folders. A plate split across several runs gets several tables in its `tables/` folder, one per run, and the importer reads all of them. That is why there is no shared file to append to and no concurrent-write problem to solve: give each run a distinct prefix and nothing collides. The existing `wants_overwrite` setting guards the rest.

## The manifest

The importer opens the table first and finds everything else from `uns["cellprofiler_mapping"]`, which already holds `channels`, `objects`, and `measurements` as `Frame`s that read back as pandas DataFrames. This module adds two entries and writes the dict unconditionally, since the importer cannot work without it.

`elements`, one row per (FOV, element) on disk:

| column | example | purpose |
|---|---|---|
| `sample_key` | `P1_A01_1` | names the image element, groups its labels, names the coordinate system |
| `image_number` | `3` | ties the row to `obs["ImageNumber"]` |
| `element_type` | `image` or `labels` | picks `Image2DModel` or `Labels2DModel` |
| `element_name` | `Nuclei` | names the labels element; empty on an image row |
| `path` | `labels/P1_A01_1/Nuclei.h5` | array location, relative to the plate folder |
| `shape`, `dtype` | `2048,2048`, `int32` | checks before parsing |
| `region_key_value` | `P1_A01_1__Nuclei` | joins to `obs["region_key"]`; empty on an image row |
| `status`, `error` | `ok`, or `failed` and a message | skip a FOV whose write failed |
| `Metadata_*` | `P1`, `A01`, `1` | whatever tags the pipeline has |

`image_channels`, one row per exported channel, giving the `c` axis of every stack: `channel` and `stack_index`. A separate two-column table rather than a delimited string on every `elements` row, so the importer sorts by index and reads off names with no parsing.

Most `elements` fields already come out of `provenance()`, which walks image-level measurements per `ImageNumber`. Both should share one row-per-image-number helper rather than traversing twice.

## Coordinate systems

Placing an element in a shared coordinate system means naming that system and giving the element a transform into it. Names carry the plate so two plates never collide: `<Plate>_<Well>_<Site>` for a FOV's own pixel space, `<Plate>_<Well>` for a well space.

Phase 1 gives every element the identity transform into its own FOV system. An image and the labels from the same FOV share it, which is what makes overlay work, and it is free because they are the same pixel grid. Naming it explicitly rather than taking SpatialData's default name `global` matters: elements from different plates all called `global` get suffixed by `concatenate()`, which is the problem in [issue #541](https://github.com/scverse/spatialdata/issues/541).

Phase 2 records per-FOV stage offsets and pixel size in `elements` when they exist as image-level measurements, so the importer can build a `<Plate>_<Well>` system. When that metadata is absent there is no stitching. The values are not guessed, defaulted, or inferred, and the module says so rather than silently producing nothing. Phase 1's per-FOV systems stay correct either way.

Plate-level coordinates are out of scope. Well pitch is a property of the plate, 9 mm on a standard 96 and 4.5 mm on a 384, and essentially never in image metadata, so a plate grid would need the user to state the format explicitly.

## What gets exported

Raw images by default, meaning every channel where `ctx.channel_info[name].source == "file"`. That filter already exists, since it drives the preview's channels table. Derived images stay out by default because pipelines produce many of them and nobody wants a folder of intermediates.

Every segmentation by default, not only the primary, secondary, and tertiary role objects. The table is role-scoped because joining compartments into one row per cell is its job, but a label array has no such constraint, and a pipeline that segments spots outside the role chain still wants them visible. Both defaults are permissive, which is right for pixels you cannot recover without re-running the pipeline.

Multi-select overrides for each. Settings otherwise mirror ExportToAnnData: `directory`, `prefix`, the role settings, `autoconfig`, the preview button, `policy`, `wants_advanced`, `wants_overwrite`. Phase 2 adds an opt-in stitching toggle, off by default.

The preview dialog gains a fifth table listing the planned elements, so the user sees the file list before starting a run.

## The table stays cell-level

One table per run, one row per cell, which is what `join_tables()` already produces. Adding a `region_key` column holding `<sample_key>__<object>` lets it annotate every labels element; `instance_key` is the existing `label_id`, whose values are already the integers in the label array.

The module does not compute a well-level aggregate, which is a deliberate divergence from the harmonized spec. A well-level table is a groupby away from a cell-level one and `spatialdata.aggregate()` exists for it, but aggregating at export time destroys the single-cell data. The spec's memory concern is real and belongs to the cell table itself, not to having several tables; per-plate folders and per-run tables already keep any one file bounded. Whether the rollup ends up in the importer or a downstream function is open.

## Failure isolation

CellProfiler already survives a bad image set: the worker catches per-image-set exceptions and only aborts on the `ED_STOP` disposition, keeping measurements recorded before the failure. Which disposition comes back depends on how the run was launched, so the module should not depend on it.

Each cycle's writes go in a try/except. On failure the module logs, records `status="failed"` and the message in that FOV's `elements` rows, and returns normally, so CellProfiler never sees an exception and the remaining FOVs still export. Writes go to a temporary name and get renamed into place, which is atomic within a filesystem, so a crash leaves no half-written file that looks valid.

## Naming invariants the importer depends on

These are the module's contract. Four things, all load-bearing for concatenating two plates without renaming anything:

- Element names carry the plate: `<Plate>_<Well>_<Site>` and `<Plate>_<Well>_<Site>__<Object>`.
- Obs names carry the plate: `<Plate>_<Well>_<Site>_<label>`, which ExportToAnnData already writes.
- Coordinate system names carry the plate.
- `region_key` is always `"region_key"` and `instance_key` always `"label_id"`, so `concatenate()` preserves both without arguments.

## The importer, in outline

Separate work, so only the shape matters here. It reads a plate folder's table, takes `elements` and `image_channels` from `uns`, parses each `ok` row into an `Image2DModel` or `Labels2DModel` under the coordinate system named by `sample_key`, and writes one `.zarr` per plate folder. Several tables in one `tables/` folder are read together.

One trap worth recording, because the module's own choices cause it: `spatialdata_attrs["region"]` must come from `sorted(obs["region_key"].unique())`, not from the labels dict. Every object is exported but the table is role-scoped, so a non-role element like `P1_A01_1__Spots` has no rows, and listing it anyway raises `ValueError: Regions in the AnnData object and cells_region do not match`, which is [issue #414](https://github.com/scverse/spatialdata/issues/414). An element with no rows is a valid outcome: the segmentation is viewable, it just has no feature vector.

## Shared backend

Reused unchanged: `introspect.build_context`, `assemble.build_object_table`, `assemble.join_tables`, `h5ad.write_h5ad` and `h5ad.Frame`, `advice.advice`, and all of `preview`'s report functions with `mapping_to_uns`. The X/obs split matters as much for a table headed into SpatialData as for ExportToAnnData. `channel_report` is a direct win, since `introspect.py` already tracks each channel's module and file-versus-pipeline source, which is what the image-selection setting needs.

Extracted rather than copy-pasted: the "build context, join tables, gather provenance, write the `.h5ad`, log advice" sequence now inline in `ExportToAnnData.post_run()` moves to `scverse_export/export.py` and both modules call it. Same reasoning that made `is_extrinsic` shared.

New: `scverse_export/raster.py` for the h5py array writer and reader, under the same h5py-only contract as `h5ad.py`, including temp-and-rename; the `elements` and `image_channels` builders; a shared sample-key helper pulled out of `assemble._obs_name`'s tag formatting, so element naming and row naming cannot drift; and `exportforspatialdata.py`.

The module needs a real `run(self, workspace)`, unlike ExportToAnnData. Measurements accumulate for the whole run and can be read in `post_run()`, but pixels do not survive that long, so each cycle's images and label arrays have to be pulled and written as it goes. The exact `Workspace`, `ImageSet`, and `ObjectSet` calls need confirming against a real `cellprofiler_core`.

## Plugin folder layout

Two Python files at the top of the plugins folder, shared code in a subfolder. `plugin_list()` in `cellprofiler_core/utilities/core/plugins.py` globs `[!_]*.py` non-recursively, so subfolders are invisible to the loader, and `load_plugins` puts the plugins directory on `sys.path`, which is what makes the subfolder importable. `load_plugin` takes the first `Module` subclass per file and warns if it finds none, so a stray top-level helper would produce a startup warning.

```
cp_plugins/                    <- CellProfiler plugins directory
  exporttoanndata.py           <- ExportToAnnData
  exportforspatialdata.py      <- ExportForSpatialData
  scverse_export/              <- shared package, never scanned
    introspect.py  names.py  assemble.py  h5ad.py  preview.py  advice.py
    raster.py                  <- new
    export.py                  <- new
  tests/  pytest.ini  README.md
```

Renaming from `cp_export` and `cpexport`, which differ by one character. Tests stay put until these go upstream, where the official repo has its own layout. The rename touches every import across roughly 16 files, all mechanical, so it should be its own commit with the suite passing unchanged on both sides.

## Why HDF5 and not Zarr

CellProfiler 4.2.8, the version installed here and the one `cellprofiler-core==4.2.8.1` pins, ships no `zarr`. Its [setup.py at that tag](https://github.com/CellProfiler/core/blob/v4.2.8.1/setup.py) lists neither `zarr` nor `spatialdata`, `ome-zarr-py`, `dask`, `xarray`, `anndata`, or `pandas`, and no newer release exists. `zarr~=2.16.1` is a dependency on the unreleased `master` branch, so CellProfiler 5 may ship it.

`numpy` and `h5py` are always present, the same constraint `h5ad.py` already works under. Every leaf array is an `.h5` file rather than a `.zarray`, but nothing about the folder layout, naming, or manifest depends on which library writes the bytes, so swapping the writer later is contained.

## Phases and order

Phase 1, the export:

1. Rename `cp_export` to `cp_plugins` and `cpexport` to `scverse_export`, suite passing unchanged either side.
2. Pull `_obs_name`'s tag formatting into a shared sample-key helper. No behavior change.
3. `scverse_export/raster.py`, unit-tested against plain numpy arrays, no CellProfiler dependency needed.
4. `scverse_export/export.py`, extracted from `ExportToAnnData.post_run()`, with its tests still passing as proof nothing changed.
5. `exportforspatialdata.py`: settings and `post_run()` first, then `run()` once the Workspace calls are confirmed.

Phase 2, stitching metadata:

6. Opt-in toggle, the metadata probe, the offset and pixel-size fields in `elements`, and the refuse-rather-than-guess path.

## Open questions

1. Whether the CellProfiler build you will run carries real `zarr`, which would let the module skip HDF5 entirely. Worth checking before step 3.
2. Where the well-level rollup eventually lives, importer or downstream. Nothing is blocked on it.
