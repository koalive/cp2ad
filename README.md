# cp2ad

Two CellProfiler plugin modules that export a pipeline's results straight into the scverse
ecosystem, no intermediate CSV.

| Module | Writes | Read it with |
|---|---|---|
| **ExportToAnnData** | one `.h5ad` per run, one row per cell | anndata, scanpy, squidpy |
| **ExportForSpatialData** | one folder per plate: image stacks, segmentation masks, and the same table | a SpatialData importer (separate) |

Point CellProfiler's plugins directory at [`cp_plugins/`](cp_plugins) and both appear under *File
Processing*. They share one backend, so the table is identical either way; ExportForSpatialData adds
the pixels around it.

Feature names and the `obs`/`uns` layout follow
`squidpy.experimental.im.calculate_image_features` (cp_measure naming), so a CellProfiler run and a
squidpy run on the same masks are interchangeable per compartment.

## What they do
- Join a pipeline's primary/secondary/tertiary objects (nuclei/cells/cytoplasm in a standard
  pipeline) into one row per cell, detecting those roles from the pipeline by default.
- Write measurements to `X` using cp_measure/squidpy-style feature names, but only the ones that
  describe an object's own morphology, intensity or texture. Position, orientation, and
  CellProfiler's own object identifiers and linkage (arbitrary object numbers, `Parent_*`/
  `Children_*` references) stay out of `X` and land in `obs` instead. `X`-derived similarity then
  reflects biology, not where a cell happened to sit in the image or which label CellProfiler
  assigned it.
- Record full pipeline provenance in `uns`: modules and their settings, object roles, channels,
  the pipeline text, and every image-level measurement in `uns["cellprofiler"]["image"]`, one entry
  per image set. That last one is where the file name, path and URL of each channel end up, so the
  export keeps track of which raw images it came from even though those are not per-cell values.
- Name rows from the pipeline's own Metadata tags (`A02_03_5` for well A02, field 3, object 5),
  either detected or pinned by hand so they stay stable across runs you mean to concatenate.

**ExportToAnnData** also shows where a name comes from before you run anything: "See where each
measurement will land" opens three tables for the per-cell columns (channels measured, objects,
measurements), naming which module produced each entry and a measurement's exact destination, plus a
fourth accounting for what the file carries outside those columns. An advanced-features checkbox
saves the first three into `uns["cellprofiler_mapping"]` as DataFrames on every run.

**ExportForSpatialData** also exports the pixels: an image stack and a label array per field of
view, written as the run goes, plus a manifest in `uns["cellprofiler_mapping"]["elements"]` giving
each file's path, shape, dtype and status. An importer builds a SpatialData object per plate from
that manifest alone, without walking the folder or parsing a file name. A field of view that fails
to write is recorded there rather than stopping the run.

See [`cp_plugins/README.md`](cp_plugins/README.md) for installation, both modules' settings, the
folder layout, the full export schema and the known limits.

## Notebook
[`test-adapter-to-viewer.ipynb`](test-adapter-to-viewer.ipynb) compares the CSV
(`ExportToSpreadsheet`) and AnnData (`ExportToAnnData`) exports of the same pipeline, and overlays
cell and nucleus centers on the source image via napari.

## Origin
A team effort from the scverse x Cell Painting Hackathon, Berlin 2026. Code base started by Tim
Treis, developed further by Loan Vulliard.
