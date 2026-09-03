# ExportToAnnData

A CellProfiler plugin that exports per-cell measurements straight from a CellProfiler pipeline to
an AnnData `.h5ad` file. Feature names and the `obs`/`uns` layout follow
`squidpy.experimental.im.calculate_image_features` (cp_measure naming), so a CellProfiler run and a
squidpy run on the same masks are interchangeable per compartment.

## What it does
- Adds an **ExportToAnnData** module to CellProfiler (point CellProfiler's plugins directory at
  [`cp_export/`](cp_export)).
- Joins a pipeline's primary/secondary/tertiary objects (nuclei/cells/cytoplasm in a standard
  pipeline) into one row per cell.
- Writes measurements to `X` using cp_measure/squidpy-style feature names, but only the ones that
  describe an object's own morphology, intensity or texture. Position, orientation, and
  CellProfiler's own object identifiers and linkage (arbitrary object numbers, `Parent_*`/
  `Children_*` references) stay out of `X` and land in `obs` instead. `X`-derived similarity then
  reflects biology, not where a cell happened to sit in the image or which label CellProfiler
  assigned it.
- Records full pipeline provenance (modules, object roles, channels, image/experiment metadata) in
  `uns`.
- Shows exactly where a name comes from before you run anything: "Press button to see where each
  object and measurement will land" opens a table of every object and measurement, which module
  produced it, and its exact destination in the export, with an optional checkbox to write that
  same information as companion CSV files on every run.

See [`cp_export/README.md`](cp_export/README.md) for installation, usage, and the full export
schema.

## Notebook
[`test-adapter-to-viewer.ipynb`](test-adapter-to-viewer.ipynb) compares the CSV
(`ExportToSpreadsheet`) and AnnData (`ExportToAnnData`) exports of the same pipeline, and overlays
cell and nucleus centers on the source image via napari.

## Origin
A team effort from the scverse x Cell Painting Hackathon, Berlin 2026. Code base started by Tim
Treis, developed further by Loan Vulliard.
