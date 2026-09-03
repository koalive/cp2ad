from cpexport.assemble import build_object_table, join_tables
from cpexport.introspect import build_context
from cpexport.preview import (DEST_DROPPED, DEST_MERGED, DEST_OBS, DEST_X, feature_report,
                              render_html, report_summary)
from test_assemble import make_measurements


def test_feature_report_matches_the_actual_export(fake_pipeline, meas_arrays):
    """feature_report() must never drift from what build_object_table/join_tables actually
    produce: every row it calls DEST_X/DEST_OBS must name a real var_name/obs key on the joined
    table, and every DEST_DROPPED row must carry no name."""
    ctx = build_context(fake_pipeline)
    m = make_measurements(meas_arrays)
    tables = {obj: build_object_table(ctx, m, obj) for obj in ctx.roles.values()}
    joined = join_tables(ctx, m, tables)

    rows = feature_report(ctx)
    assert rows  # the fixture pipeline has real features; an empty report would hide a bug
    for r in rows:
        if r.destination == DEST_X:
            assert r.anndata_name in joined.var_names, r
        elif r.destination in (DEST_OBS, DEST_MERGED):
            assert r.anndata_name in joined.obs, r
        else:
            assert r.destination == DEST_DROPPED and r.anndata_name == "", r


def test_feature_report_scoped_to_role_objects(fake_pipeline, meas_arrays):
    """Only the objects ctx.roles names should appear -- an object the pipeline made but that is
    not part of the export (not a role) has nothing to preview."""
    ctx = build_context(fake_pipeline)
    rows = feature_report(ctx)
    assert {r.object for r in rows} <= set(ctx.roles.values())


def test_report_summary_counts_add_up(fake_pipeline):
    ctx = build_context(fake_pipeline)
    rows = feature_report(ctx)
    summary = report_summary(rows)
    assert sum(summary.values()) == len(rows)
    assert set(summary) == {DEST_X, DEST_OBS, DEST_MERGED, DEST_DROPPED}


def test_render_html_is_well_formed_and_filters(fake_pipeline):
    ctx = build_context(fake_pipeline)
    rows = feature_report(ctx)
    html = render_html(rows)
    assert html.startswith("<html>") and html.endswith("</html>")
    assert html.count("<tr") == len(rows) + 1  # +1 header row

    some_object = rows[0].object
    filtered = render_html(rows, filter_text=some_object)
    expected = sum(1 for r in rows if some_object.lower() in
                   " ".join([r.object, r.cp_name, r.module_name, r.category, r.destination,
                             r.anndata_name, r.reason]).lower())
    assert filtered.count("<tr") == expected + 1

    # a filter matching nothing still renders a valid (empty-body) table, not an error
    empty = render_html(rows, filter_text="no such feature exists 12345")
    assert empty.count("<tr") == 1
