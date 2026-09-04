"""Backend shared by ExportToAnnData and ExportForSpatialData.

The two plugin files next to this package hold the CellProfiler Module subclasses, because
CellProfiler's plugin loader globs the top level of the plugins folder only. Everything they can
share lives here:

    introspect  read roles, channels and objects off a pipeline, without measurements
    assemble    Measurements -> squidpy-shaped tables
    names       CellProfiler feature identity -> column name, and the X/obs split
    samples     sample keys and row names, built from Metadata tags
    export      the table-building step both modules call
    h5ad        .h5ad writer on h5py alone, since the CellProfiler app has no anndata
    raster      image and label arrays on disk, one array per file
    spatial     ExportForSpatialData's folder layout and manifest
    advice      what to change in the pipeline, as warnings and log lines
    preview     ExportToAnnData's "where will each measurement land" tables

No module here imports CellProfiler. They take pipelines, measurements and workspaces as plain
arguments, so the whole package is testable on numpy and h5py.
"""
