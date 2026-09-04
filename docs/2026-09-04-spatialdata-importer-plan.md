# Where the SpatialData importer belongs

> **Status, 2026-09-04.** Built and working. `read_cellprofiler_export` and `cellprofiler_export_plates` are in `cell-painting-io` at `src/cell_painting_io/cellprofiler.py`, with 17 tests in `tests/test_cellprofiler.py`. Validated on the real four-field export: 4 Images, 8 Labels, a 151 x 126 table, 4 coordinate systems, the label-to-row join exact for every field, and a zarr round-trip preserving `X` and the provenance. Two exporter fixes were needed and are made but not yet committed, since they need a GUI run. Naming alignment with `read_plate` stays open by decision, and is documented in both repos rather than worked around.



Both candidate repos are checked out at `~/Documents/Tests/`. I read them against what ExportForSpatialData actually writes.

## Recommendation: `cell-painting-io`, as a second reader

`cell_painting_io.read_plate` already builds the object we are aiming at. One plate of the Cell Painting Gallery becomes fields of view as Images, the CellProfiler Nuclei, Cells and Cytoplasm as Labels, wells as Shapes, and well- and cell-level measurements as Tables, with every element in a field, well and plate coordinate system and the plate barcode in its name so two plates concatenate. Our exporter produces the same kind of object from a different source, so the importer is a sibling of `read_plate`, not a change to it.

What it reuses from that repo: `parse_well`, `PLATE_FORMATS`, `fov_offsets` and `_well_shapes` for phase 2, the `{plate}_{well}_s{site}` / `{plate}_{well}` / `{plate}` coordinate-system convention, plate-prefixed element names, and the `TableModel.parse` pattern. What it does not need: `labels_from_outlines`, because we export real masks rather than the gallery's one-pixel outlines, and `read_profiles` / `annotate_features`, because our table arrives as an AnnData with `var` already annotated by the exporter.

`spatialdata-io` is the wrong home for now, on its own stated terms. Its [contributing guide](https://github.com/scverse/spatialdata-io/blob/main/docs/contributing.md) asks a new reader for a public specification of the raw layout, small public example data with a permissive licence, and tests. We have none of the three: the plugin is unreleased, the `.h5` layout is our own invention and undocumented outside this repo, and no public dataset exists in it. Their guide routes exactly this case to the `experimental` module, where `iss` sits. That is where to propose it later, once the plugin is merged into CellProfiler-plugins and the layout is written down. Proving the format in `cell-painting-io` first costs nothing and makes that proposal stronger.

## Form

One function per plate folder, in a new `src/cell_painting_io/cellprofiler.py`:

```python
def read_cellprofiler_export(path: Path | str, *, lazy: bool = True) -> SpatialData
```

Exposed from `__init__.py` through the same lazy `__getattr__` that defers `read_plate`, so importing the package still does not pull in `spatialdata`. A second helper lists the plate folders under one export root, so a multi-plate run reads as a list of objects that `spatialdata.concatenate()` accepts.

It is manifest-driven and nothing else. Read `tables/<prefix>.h5ad`, take `uns["cellprofiler_mapping"]["elements"]`, and build only what those rows name, resolving each `path` relative to the plate folder. No directory walking, no filename parsing. Rows with `status == "failed"` are skipped and reported, which is what the status column is for. `image_channels` gives `c_coords` in stack order. `obs["region_key"]` joins rows to Labels elements.

Lazy by default is worth having: every leaf is one h5py dataset, so `dask.array.from_array(h5py.File(p)[ "data"])` defers the read and matches what the spatialdata-io guide asks of raster loading.

## What the exporter had to change

Both of these came out of trying to build a real SpatialData object rather than reading the code, and neither is visible from inside the exporter's own tests.

**1. `uns["spatialdata_attrs"]` is wrong on ExportForSpatialData's table.** Verified against the file from the last run:

```
spatialdata_attrs: {'region': 'Cells', 'region_key': 'region', 'instance_key': 'label_id'}
obs['region']      : all 'Cells'
obs['region_key']  : 'A02_03_img1__Cells', 'B04_01_img2__Cells', ...
labels elements    : 'A02_03_img1__Cells', 'A02_03_img1__Nuclei', ...
```

No element is called `Cells`, so the attrs point at a region that does not exist, which is [spatialdata issue #414](https://github.com/scverse/spatialdata/issues/414) waiting to happen. The join column is correct and complete: every `region_key` value is a real Labels element. Only the attrs were wrong.

It is stricter than a wrong label. `TableModel.parse` refuses outright when the key is already set, `region`, `region_key` and / or `instance_key` is/has been passed as argument(s), so the exported table could not be parsed at all without deleting the attrs first. Fixed by `spatial.spatialdata_attrs()`, which names the label elements the rows annotate, called per plate in `post_run`.

**2. A manifest column named `dtype` made the object unwritable to zarr.** This one only appears at `sdata.write()`:

```
AttributeError: 'Series' object has no attribute 'kind'
Error raised while writing key 'elements' of ... to /tables/cells/uns/cellprofiler_mapping
```

Pandas resolves `frame.dtype` to a column of that name, and anndata's writer dispatches on `elem.dtype.kind`, so it got a Series where it expected a numpy dtype. The manifest travels in the table's `uns`, so one unwritable column made the whole SpatialData object unwritable. A column named `shape` is safe, because the real `DataFrame.shape` attribute shadows it; I checked the other four manifest tables and only `elements` was affected. The column is `element_dtype` now, with a test asserting no manifest column is named `dtype`, since nothing upstream of the zarr write catches it.

## Two divergences from `read_plate`, left in place on purpose

**Element naming.** Ours is `{sample_key}__{Object}`, e.g. `A02_03_img1__Cells`; theirs is `{plate}_{well}_s{site}_{object}` with a lowercase object, e.g. `BR00000001_A01_s1_cells`. Images are `{field}_image` in both, which is the one thing that already agrees. Exact string agreement is not reachable anyway, because our key adapts to whichever Metadata tags a pipeline carries. Structural agreement is: element name = field-of-view key plus role, one separator, lowercase role. Deferred by decision, and stated in both READMEs and in the reader's docstring so nobody assumes an object from the gallery and an object from a run name the same thing the same way.

**Coordinate systems.** `read_plate` gives every element three, field, well and plate, built from stage coordinates. We give one, the per-field identity, because the module exports no stage coordinates or pixel size. That is phase 2 of the ExportForSpatialData plan. Until then a plate reads as unstitched fields, which the reader says plainly rather than faking with a nominal grid.

## What is left

1. A GUI run to confirm the two exporter fixes, then commit them here.
2. Naming alignment, when it is decided.
3. A round-trip test running both exports over one pipeline and asserting the per-cell tables agree. It needs a real CellProfiler run, so it belongs here rather than in `cell-painting-io`.
4. Phase 2, once the exporter records stage offsets: well and plate systems through `fov_offsets`, and well Shapes through `_well_shapes`.
5. `instance_key` is `label_id` here and `ObjectNumber` in `read_plate`. Both hold the CellProfiler object number. Part of the naming decision.
