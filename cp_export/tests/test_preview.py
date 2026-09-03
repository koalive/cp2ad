from cpexport.assemble import build_object_table, join_tables
from cpexport.introspect import build_context
from cpexport.preview import (DEST_MERGED, DEST_OBS, DEST_X, channel_report, mapping_to_uns,
                              measurement_report, object_report, render_channel_table,
                              render_object_table, render_report_html, report_summary,
                              uns_report)
from test_assemble import make_measurements


def test_measurement_report_matches_the_actual_export(fake_pipeline, meas_arrays):
    """measurement_report() must never drift from what build_object_table/join_tables actually
    produce: every row must name a real var_name (DEST_X) or obs key (DEST_OBS/DEST_MERGED) on the
    joined table. A non-numeric CellProfiler column gets no row at all, never a row with an empty
    name."""
    ctx = build_context(fake_pipeline)
    m = make_measurements(meas_arrays)
    tables = {obj: build_object_table(ctx, m, obj) for obj in ctx.roles.values()}
    joined = join_tables(ctx, m, tables)

    rows = measurement_report(ctx)
    assert rows  # the fixture pipeline has real measurements; an empty report would hide a bug
    for r in rows:
        assert r.anndata_name
        if r.destination == DEST_X:
            assert r.anndata_name in joined.var_names, r
        else:
            assert r.destination in (DEST_OBS, DEST_MERGED), r
            assert r.anndata_name in joined.obs, r


def test_measurement_report_scoped_to_role_objects(fake_pipeline, meas_arrays):
    """Only the objects ctx.roles names appear. An object the pipeline made but that isn't part of
    the export has nothing to report."""
    ctx = build_context(fake_pipeline)
    rows = measurement_report(ctx)
    assert {r.object for r in rows} <= set(ctx.roles.values())


def test_report_summary_counts_add_up(fake_pipeline):
    ctx = build_context(fake_pipeline)
    rows = measurement_report(ctx)
    summary = report_summary(rows)
    assert sum(summary.values()) == len(rows)
    assert set(summary) == {DEST_X, DEST_OBS, DEST_MERGED}


def test_object_report_names_the_producing_module(fake_pipeline):
    """object_report() answers "where does this object come from": one row per role, naming the
    module that produced it."""
    ctx = build_context(fake_pipeline)
    rows = object_report(ctx)
    assert {r.role for r in rows} <= {"primary", "secondary", "tertiary"}
    assert {r.object for r in rows} == set(ctx.roles.values())
    for r in rows:
        info = ctx.objects[r.object]
        assert r.module_name == info.module_name
        assert r.module_num == info.module_num
        assert r.source in ("pipeline", "file")


def test_channel_report_names_the_producing_module(fake_pipeline):
    """channel_report() answers "where does this channel come from": one row per channel, naming
    the module that loaded or computed it."""
    ctx = build_context(fake_pipeline)
    rows = channel_report(ctx)
    assert rows  # the fixture pipeline measures at least one channel
    for r in rows:
        info = ctx.channel_info[r.channel]
        assert r.module_name == info.module_name
        assert r.module_num == info.module_num
        assert r.source == info.source
        assert r.source in ("pipeline", "file")


def test_channel_report_excludes_channels_no_export_measurement_reads(fake_pipeline):
    """A channel only feeding a derived image (never itself measured) gets no row: channel_report
    lists what's exported, not every named image the pipeline has."""
    ctx = build_context(fake_pipeline)
    measured = {r.channel for r in measurement_report(ctx) if r.channel}
    measured |= {r.channel2 for r in measurement_report(ctx) if r.channel2}
    rows = channel_report(ctx)
    assert {r.channel for r in rows} == measured
    assert {r.channel for r in rows} <= set(ctx.channels)


def test_render_channel_table_reports_when_nothing_is_measured():
    assert "No exported measurement reads a channel" in render_channel_table([])


def test_render_report_html_has_all_four_tables_and_filters(fake_pipeline):
    ctx = build_context(fake_pipeline)
    channel_rows = channel_report(ctx)
    object_rows = object_report(ctx)
    measurement_rows = measurement_report(ctx)
    uns_rows = uns_report(ctx)
    html = render_report_html(channel_rows, object_rows, measurement_rows, uns_rows)
    assert html.startswith("<html>") and html.endswith("</html>")
    assert (html.index("<h3>Channels measured</h3>") < html.index("<h3>Objects</h3>")
            < html.index("<h3>Measurements</h3>") < html.index("<h3>Also exported, in uns</h3>"))
    # one header row per table, plus one row per channel/object/measurement/uns key
    assert html.count("<tr") == (len(channel_rows) + len(object_rows) + len(measurement_rows)
                                 + len(uns_rows) + 4)

    some_object = measurement_rows[0].object
    filtered = render_report_html(channel_rows, object_rows, measurement_rows, uns_rows,
                                  filter_text=some_object)
    assert some_object in filtered
    assert filtered.count("<tr") <= html.count("<tr")

    # a filter matching nothing still renders valid, header-only tables, not an error
    empty = render_report_html(channel_rows, object_rows, measurement_rows, uns_rows,
                               filter_text="no such thing exists 12345")
    assert empty.count("<tr") == 4


def test_uns_report_accounts_for_the_image_level_columns(fake_pipeline):
    """The per-cell tables omit image-level columns, so the uns table has to account for them.
    FileName_*/PathName_* are exported, to uns["cellprofiler"]["image"], and saying otherwise
    (as an earlier version of the preview did) is wrong."""
    ctx = build_context(fake_pipeline)
    rows = uns_report(ctx)
    image_row = next(r for r in rows if r.key == 'cellprofiler["image"]')
    n_image_feats = len([f for f in ctx.features if f.object == "Image"])
    assert image_row.detail == f"{n_image_feats} columns"
    if any(f.cp_name.startswith("FileName_") for f in ctx.features if f.object == "Image"):
        assert "FileName_*" in image_row.holds

    # every channel stays named in uns, including the ones "Channels measured" filters out
    channels_row = next(r for r in rows if r.key == 'cellprofiler["channels"]')
    assert channels_row.detail == f"{len(ctx.channels)} channels"
    assert len(channel_report(ctx)) <= len(ctx.channels)


def test_uns_report_lists_the_mapping_tables_only_when_enabled(fake_pipeline):
    ctx = build_context(fake_pipeline)
    assert not any(r.key == "cellprofiler_mapping" for r in uns_report(ctx))
    assert any(r.key == "cellprofiler_mapping" for r in uns_report(ctx, wants_mapping_uns=True))


def test_render_object_table_reports_when_nothing_is_configured():
    assert "No object is currently assigned" in render_object_table([])


def test_mapping_to_uns_round_trips_as_dataframes(fake_pipeline, tmp_path):
    """These tables are tabular, so they must come back as DataFrames, not as a dict keyed "0",
    "1", "2" (what _write_elem does with a plain list of row dicts). A dict cannot be filtered or
    sorted without being rebuilt first."""
    import numpy
    import pandas
    anndata_lib = __import__("anndata")
    from cpexport.assemble import Table
    from cpexport.h5ad import write_h5ad

    ctx = build_context(fake_pipeline)
    mapping = mapping_to_uns(channel_report(ctx), object_report(ctx), measurement_report(ctx))

    t = Table(X=numpy.zeros((1, 1), dtype=numpy.float32), obs_names=["a"], var_names=["x"],
             obsm={"spatial": numpy.zeros((1, 2))}, uns={"cellprofiler_mapping": mapping})
    path = tmp_path / "mapping.h5ad"
    write_h5ad(t, str(path))

    back = anndata_lib.read_h5ad(str(path)).uns["cellprofiler_mapping"]
    assert set(back.keys()) == {"channels", "objects", "measurements"}
    for name, frame in mapping.items():
        table = back[name]
        assert isinstance(table, pandas.DataFrame), (name, type(table))
        assert list(table.columns) == list(frame.columns)          # declared column order survives
        assert len(table) == len(next(iter(frame.columns.values())))

    measurements = back["measurements"]
    assert measurements["destination"].isin({DEST_X, DEST_OBS, DEST_MERGED}).all()
    # the point of a DataFrame: this is a one-liner rather than a rebuild
    assert len(measurements[measurements.destination == DEST_X]) >= 1
