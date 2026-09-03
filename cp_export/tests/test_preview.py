import csv

from cpexport.assemble import build_object_table, join_tables
from cpexport.introspect import build_context
from cpexport.preview import (DEST_DROPPED, DEST_MERGED, DEST_OBS, DEST_X, measurement_report,
                              object_report, render_object_table, render_report_html,
                              report_summary, write_csv)
from test_assemble import make_measurements


def test_measurement_report_matches_the_actual_export(fake_pipeline, meas_arrays):
    """measurement_report() must never drift from what build_object_table/join_tables actually
    produce: every DEST_X/DEST_OBS row must name a real var_name/obs key on the joined table, and
    every DEST_DROPPED row must carry no name."""
    ctx = build_context(fake_pipeline)
    m = make_measurements(meas_arrays)
    tables = {obj: build_object_table(ctx, m, obj) for obj in ctx.roles.values()}
    joined = join_tables(ctx, m, tables)

    rows = measurement_report(ctx)
    assert rows  # the fixture pipeline has real measurements; an empty report would hide a bug
    for r in rows:
        if r.destination == DEST_X:
            assert r.anndata_name in joined.var_names, r
        elif r.destination in (DEST_OBS, DEST_MERGED):
            assert r.anndata_name in joined.obs, r
        else:
            assert r.destination == DEST_DROPPED and r.anndata_name == "", r


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
    assert set(summary) == {DEST_X, DEST_OBS, DEST_MERGED, DEST_DROPPED}


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


def test_render_report_html_has_both_tables_and_filters(fake_pipeline):
    ctx = build_context(fake_pipeline)
    object_rows = object_report(ctx)
    measurement_rows = measurement_report(ctx)
    html = render_report_html(object_rows, measurement_rows)
    assert html.startswith("<html>") and html.endswith("</html>")
    assert "<h3>Objects</h3>" in html and "<h3>Measurements</h3>" in html
    # one header row per table, plus one row per object/measurement
    assert html.count("<tr") == len(object_rows) + len(measurement_rows) + 2

    some_object = measurement_rows[0].object
    filtered = render_report_html(object_rows, measurement_rows, filter_text=some_object)
    assert some_object in filtered
    assert filtered.count("<tr") <= html.count("<tr")

    # a filter matching nothing still renders valid, header-only tables, not an error
    empty = render_report_html(object_rows, measurement_rows, filter_text="no such thing exists 12345")
    assert empty.count("<tr") == 2


def test_render_object_table_reports_when_nothing_is_configured():
    assert "No object is currently assigned" in render_object_table([])


def test_write_csv_round_trips(fake_pipeline, tmp_path):
    ctx = build_context(fake_pipeline)
    rows = object_report(ctx)
    path = tmp_path / "objects.csv"
    write_csv(rows, str(path))
    with open(path, newline="") as fh:
        read_back = list(csv.DictReader(fh))
    assert len(read_back) == len(rows)
    assert read_back[0]["object"] == rows[0].object
    assert read_back[0]["module_name"] == rows[0].module_name
