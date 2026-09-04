# ExportToAnnData

A CellProfiler plugin that exports per-cell measurements straight from a CellProfiler pipeline to
an AnnData `.h5ad` file. Feature names and the `obs`/`uns` layout follow
`squidpy.experimental.im.calculate_image_features` (cp_measure naming), so a CellProfiler run and a
squidpy run on the same masks are interchangeable per compartment.

## What it does
- Adds an **ExportToAnnData** module to CellProfiler (point CellProfiler's plugins directory at
  [`cp_plugins/`](cp_plugins)).
- Joins a pipeline's primary/secondary/tertiary objects (nuclei/cells/cytoplasm in a standard
  pipeline) into one row per cell.
- Writes measurements to `X` using cp_measure/squidpy-style feature names, but only the ones that
  describe an object's own morphology, intensity or texture. Position, orientation, and
  CellProfiler's own object identifiers and linkage (arbitrary object numbers, `Parent_*`/
  `Children_*` references) stay out of `X` and land in `obs` instead. `X`-derived similarity then
  reflects biology, not where a cell happened to sit in the image or which label CellProfiler
  assigned it.
- Records full pipeline provenance in `uns`: modules and their settings, object roles, channels,
  the pipeline text, and every image-level measurement in `uns["cellprofiler"]["image"]`, one entry
  per image set. That last one is where the file name, path and URL of each channel end up, so the
  export keeps track of which raw images it came from even though those are not per-cell values.
- Shows exactly where a name comes from before you run anything: "See where each measurement will
  land" opens three tables for the per-cell columns (channels measured, objects, measurements),
  naming which module produced each entry and a measurement's exact destination, plus a fourth
  accounting for what the file carries outside those columns. An advanced-features checkbox saves
  the first three into `uns["cellprofiler_mapping"]` as DataFrames on every run.

See [`cp_plugins/README.md`](cp_plugins/README.md) for installation, usage, and the full export
schema.

## Notebook
[`test-adapter-to-viewer.ipynb`](test-adapter-to-viewer.ipynb) compares the CSV
(`ExportToSpreadsheet`) and AnnData (`ExportToAnnData`) exports of the same pipeline, and overlays
cell and nucleus centers on the source image via napari.

## Origin
A team effort from the scverse x Cell Painting Hackathon, Berlin 2026. Code base started by Tim
Treis, developed further by Loan Vulliard.
