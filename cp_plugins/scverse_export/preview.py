"""Report of where each channel, object, and measurement lands in the AnnData export: which module
made it, and for measurements, whether it ends up in X or obs (merged into an existing column
counts as obs), and under what exact name.

The channel, object and measurement tables cover the per-cell part of the export, and list only
what reaches it: a channel no measurement reads, an object with no role, a per-object column that
is not numeric. None of those get a row. Anything that is exported but is not a per-cell column
still shows up, in the fourth table uns_report builds. Per-image-set facts in particular, which
file each channel was read from included, are exported to uns["cellprofiler"]["image"] rather than
being dropped.

Powers ExportToAnnData's "See where each measurement will land" dialog. Pure Python, no GUI or
CellProfiler-runtime dependency, so it stays testable outside the module.

The measurement report replays assemble.py's own is_numeric/is_extrinsic/to_cpm_names rules and its
per-object dedup, instead of re-deriving them, so it can't drift from what a run actually produces.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, fields
from html import escape
from typing import Dict, List

from .assemble import _is_numeric
from .h5ad import Frame
from .introspect import IMAGE, Context, file_loaded_objects
from .names import is_extrinsic, to_cpm_names

DEST_X, DEST_OBS, DEST_MERGED = "X", "obs", "merged"

REASON_EXTRINSIC = "position, orientation, or identity/linkage, not this object's own biology"


@dataclass(frozen=True)
class ChannelReportRow:
    channel: str
    module_name: str
    module_num: int
    source: str                 # "file" (loaded by NamesAndTypes/LoadData/LoadImages) or "pipeline"


@dataclass(frozen=True)
class ObjectReportRow:
    role: str                   # primary | secondary | tertiary
    object: str
    module_name: str
    module_num: int
    source: str                 # "pipeline" (an Identify* module) or "file" (a label image)


@dataclass(frozen=True)
class MeasurementReportRow:
    object: str
    cp_name: str                # CellProfiler's own measurement name, e.g. "AreaShape_Area"
    module_name: str
    module_num: int
    category: str
    channel: str                 # the image this measurement reads, or "" if it isn't per-channel
    channel2: str                # the second image, for a Correlation pair; "" otherwise
    destination: str            # DEST_X | DEST_OBS | DEST_MERGED
    anndata_name: str           # "<object>__<name>" on the joined file
    reason: str                 # empty for a plain X measurement; explains obs/merged otherwise


def channel_report(ctx: Context) -> List[ChannelReportRow]:
    """One row per channel that at least one exported measurement actually reads, in pipeline
    order. Answers "where does LogDNA come from": which module produced it, and whether that
    module loaded it from a file or computed it from other images. A channel nothing measures
    (e.g. a raw image only used to build a derived one) gets no row: it was never exported."""
    used = set()
    for row in measurement_report(ctx):
        if row.channel:
            used.add(row.channel)
        if row.channel2:
            used.add(row.channel2)
    rows = []
    for name in ctx.channels:
        if name not in used:
            continue
        info = ctx.channel_info.get(name)
        if info is None:
            rows.append(ChannelReportRow(channel=name, module_name="?", module_num=-1, source="?"))
            continue
        rows.append(ChannelReportRow(channel=name, module_name=info.module_name,
                                     module_num=info.module_num, source=info.source))
    return rows


def object_report(ctx: Context) -> List[ObjectReportRow]:
    """One row per role CellProfiler has assigned an object to. Answers "where does Nuclei come
    from" directly: which module produced it, and whether that module identified it from an image
    or loaded it from a label file."""
    from_files = set(file_loaded_objects(ctx))
    rows = []
    for role in ("primary", "secondary", "tertiary"):
        obj = ctx.roles.get(role)
        if obj is None:
            continue
        info = ctx.objects[obj]
        rows.append(ObjectReportRow(
            role=role, object=obj, module_name=info.module_name or "?",
            module_num=-1 if info.module_num is None else info.module_num,
            source="file" if obj in from_files else "pipeline"))
    return rows


def _measurement_row(f, destination: str, anndata_name: str, reason: str) -> MeasurementReportRow:
    return MeasurementReportRow(
        object=f.object, cp_name=f.cp_name,
        module_name=f.module_name or "?", module_num=-1 if f.module_num is None else f.module_num,
        category=f.category, channel=f.image or "", channel2=f.image2 or "",
        destination=destination, anndata_name=anndata_name, reason=reason)


def measurement_report(ctx: Context) -> List[MeasurementReportRow]:
    """One row per (object, measurement, cpm name) the export machinery actually touches, scoped to
    the objects ctx.roles names, i.e. exactly what reaches the joined .h5ad as a per-cell column.
    A per-object column that is not numeric (a varchar or blob measurement) cannot enter a float
    matrix, so it gets no row here.

    Image-level columns are a different matter and are not in scope here at all: FileName_*,
    PathName_*, Metadata_*, Count_*, Threshold_* and the rest belong to the Image object, not to a
    role object, and provenance() exports them to uns["cellprofiler"]["image"]. See uns_report.

    anndata_name is the joined file's name (<object>__<name>); a per-object file, if the module
    also writes one, uses the same name without the object prefix."""
    rows: List[MeasurementReportRow] = []
    role_objects = list(dict.fromkeys(ctx.roles.values()))  # role order, de-duplicated
    for obj in role_objects:
        claimed_by: Dict[tuple, str] = {}  # (is_extrinsic, cpm name) -> claiming cp_name
        for f in ctx.features:
            if f.object != obj or not _is_numeric(f.coltype):
                continue
            extrinsic = is_extrinsic(f)
            for name, backend in to_cpm_names(f, ctx.channels):
                key = (extrinsic, name)
                anndata_name = f"{obj}__{name}"
                if key in claimed_by:
                    rows.append(_measurement_row(f, DEST_MERGED, anndata_name,
                                                 f"same value as {claimed_by[key]} under this name"))
                    continue
                claimed_by[key] = f.cp_name
                rows.append(_measurement_row(f, DEST_OBS if extrinsic else DEST_X, anndata_name,
                                             REASON_EXTRINSIC if extrinsic else ""))
    return rows


def report_summary(rows: List[MeasurementReportRow]) -> Dict[str, int]:
    """Counts per destination, for the one-line summary above the measurement table."""
    out = {DEST_X: 0, DEST_OBS: 0, DEST_MERGED: 0}
    for r in rows:
        out[r.destination] += 1
    return out


def _frame(rows, row_type) -> Frame:
    return Frame.from_records([asdict(r) for r in rows], [f.name for f in fields(row_type)])


def mapping_to_uns(channel_rows: List[ChannelReportRow], object_rows: List[ObjectReportRow],
                   measurement_rows: List[MeasurementReportRow]) -> Dict[str, Frame]:
    """The same three tables the preview dialog shows, for uns["cellprofiler_mapping"].

    Each one is a Frame, so write_h5ad stores it in anndata's dataframe encoding and it reads back
    as a pandas DataFrame. These tables are tabular, and a dict of row-dicts keyed "0", "1", "2"
    (what _write_elem does with a plain list of dicts) cannot be sliced, sorted, or filtered
    without rebuilding it first."""
    return {
        "channels": _frame(channel_rows, ChannelReportRow),
        "objects": _frame(object_rows, ObjectReportRow),
        "measurements": _frame(measurement_rows, MeasurementReportRow),
    }


# ---- what else the export carries in uns -------------------------------------------------------

@dataclass(frozen=True)
class UnsReportRow:
    key: str                    # the uns path, e.g. 'cellprofiler["image"]'
    holds: str                  # what lives there, in the reader's terms
    detail: str                 # how much of it, counted from this pipeline where possible


def uns_report(ctx: Context, wants_mapping_uns: bool = False) -> List[UnsReportRow]:
    """What the export puts in uns, beyond X/obs/var. Answers the question the other three tables
    cannot: measurements that never become a per-cell column are not lost, most of them are
    per-image-set facts (which file each channel was read from, plate/well metadata, thresholds)
    and they land in uns["cellprofiler"]["image"] instead, one row per image set."""
    image_feats = [f.cp_name for f in ctx.features if f.object == IMAGE]
    located = sorted({n.split("_", 1)[0] for n in image_feats
                      if n.startswith(("FileName_", "PathName_", "URL_"))})
    rows = [
        UnsReportRow(
            key='cellprofiler["image"]',
            holds="One row per image set: where each channel was read from (" +
                  (", ".join(f"{p}_*" for p in located) if located else "no file columns found") +
                  "), plate/well metadata, object counts, thresholds, per-module timings",
            detail=f"{len(image_feats)} columns"),
        UnsReportRow(
            key='cellprofiler["channels"]',
            holds="Every channel name the pipeline defines, in order, including any no measurement "
                  "reads (so a channel absent from the Channels table above still appears here)",
            detail=f"{len(ctx.channels)} channels"),
        UnsReportRow(
            key='cellprofiler["objects"], ["roles"], ["role_detection"]',
            holds="Every object the pipeline makes, which module made it, and how the "
                  "primary/secondary/tertiary roles were chosen",
            detail=f"{len(ctx.objects)} object{'' if len(ctx.objects) == 1 else 's'}, "
                   f"{len(ctx.roles)} with a role"),
        UnsReportRow(
            key='cellprofiler["modules"]',
            holds="Every module in the pipeline with all of its settings, in order",
            detail=f"{len(ctx.modules)} modules"),
        UnsReportRow(
            key='cellprofiler["pipeline_text"], ["version"], ["run_timestamp"], ["experiment"]',
            holds="The full pipeline as text, the CellProfiler version, and when the run happened",
            detail="filled in at run time"),
        UnsReportRow(
            key='cellprofiler["relationships"]',
            holds="Parent/child object relationships CellProfiler recorded (RelateObjects and "
                  "similar)",
            detail="filled in at run time"),
        UnsReportRow(
            key="cellprofiler_join, qc_summary",
            holds="How each compartment was matched onto the base row, and how many rows were not "
                  "exactly one primary plus one secondary plus one tertiary",
            detail="joined file only"),
        UnsReportRow(
            key="spatialdata_attrs",
            holds="Region and instance keys, so the table loads into SpatialData as a TableModel",
            detail="3 keys"),
    ]
    if wants_mapping_uns:
        rows.append(UnsReportRow(
            key="cellprofiler_mapping",
            holds="The Channels, Objects and Measurements tables above, as three DataFrames",
            detail="this setting is on"))
    return rows


# ---- HTML rendering (for the preview dialog's wx.html.HtmlWindow) ----
#
# wx.html.HtmlWindow is a lightweight HTML 3.2-ish renderer with only partial CSS support, and no
# external browser engine (unlike wx.html2.WebView, which needs one). Colour-coding here uses the
# classic <td bgcolor="..."> / <font color="..."> attributes it understands natively; a <style>
# block or inline CSS renders inconsistently across wx versions.
#
# Font size gets two independent pushes, because a wx.html.HtmlWindow's default text still reads
# small even after the dialog raises the widget's base size (HtmlWindow.SetStandardFonts): the
# whole body is also wrapped in an explicit <font size="4">, one HTML level above the "normal"
# level 3, so the page is legible even if the widget-level setting has no visible effect on a given
# platform.
_BODY_FONT_SIZE = "4"

_DEST_STYLE = {
    DEST_X: ("X (morphology)", "#e6f4ea", "#1e7e34"),
    DEST_OBS: ("obs (extrinsic)", "#e8f0fe", "#1a56c4"),
    DEST_MERGED: ("obs (merged)", "#fff6e0", "#8a6d1a"),
}


def _measurement_matches(row: MeasurementReportRow, needle: str) -> bool:
    haystack = " ".join([row.object, row.cp_name, row.module_name, row.category,
                         row.destination, row.anndata_name, row.reason]).lower()
    return needle in haystack


def _object_matches(row: ObjectReportRow, needle: str) -> bool:
    haystack = " ".join([row.role, row.object, row.module_name, row.source]).lower()
    return needle in haystack


def _channel_matches(row: ChannelReportRow, needle: str) -> bool:
    haystack = " ".join([row.channel, row.module_name, row.source]).lower()
    return needle in haystack


def render_channel_table(rows: List[ChannelReportRow], filter_text: str = "") -> str:
    """An HTML table (no <html>/<body> wrapper) mapping each channel some exported measurement
    reads to the module that loaded or computed it. A channel with no such measurement (a raw
    image only used to build a derived one, say) is not here; it is still named in
    uns["cellprofiler"]["channels"], and where its file came from is in ["image"]."""
    needle = filter_text.strip().lower()
    shown = [r for r in rows if not needle or _channel_matches(r, needle)]
    parts = ["<h3>Channels measured</h3>"]
    if not rows:
        parts.append("<p>No exported measurement reads a channel directly.</p>")
        return "".join(parts)
    parts.append('<table border="1" cellpadding="5" cellspacing="0">'
                 '<tr bgcolor="#dddddd"><th>Channel</th><th>Defined in module</th>'
                 "<th>Source</th></tr>")
    for r in shown:
        module = f"{escape(r.module_name)} (#{r.module_num})" if r.module_num >= 0 else escape(r.module_name)
        parts.append(f"<tr><td><b>{escape(r.channel)}</b></td><td>{module}</td>"
                     f"<td>{escape(r.source)}</td></tr>")
    parts.append("</table>")
    if needle and len(shown) < len(rows):
        parts.append(f"<p>Showing {len(shown)} of {len(rows)} channels.</p>")
    return "".join(parts)


def render_object_table(rows: List[ObjectReportRow], filter_text: str = "") -> str:
    """An HTML table (no <html>/<body> wrapper, meant to be embedded) mapping each role to the
    object filling it, the module that made it, and whether that module identified it from an
    image or loaded it from a label file."""
    needle = filter_text.strip().lower()
    shown = [r for r in rows if not needle or _object_matches(r, needle)]
    parts = ["<h3>Objects</h3>"]
    if not rows:
        parts.append("<p>No object is currently assigned to a role, so nothing will be exported.</p>")
        return "".join(parts)
    parts.append('<table border="1" cellpadding="5" cellspacing="0">'
                 '<tr bgcolor="#dddddd"><th>Role</th><th>Object</th><th>Defined in module</th>'
                 "<th>Source</th></tr>")
    for r in shown:
        module = f"{escape(r.module_name)} (#{r.module_num})" if r.module_num >= 0 else escape(r.module_name)
        parts.append(f"<tr><td>{escape(r.role)}</td><td><b>{escape(r.object)}</b></td><td>{module}</td>"
                     f"<td>{escape(r.source)}</td></tr>")
    parts.append("</table>")
    if needle and len(shown) < len(rows):
        parts.append(f"<p>Showing {len(shown)} of {len(rows)} objects.</p>")
    return "".join(parts)


def render_measurement_table(rows: List[MeasurementReportRow], filter_text: str = "") -> str:
    """An HTML table (no <html>/<body> wrapper) of every measurement that becomes a per-cell
    column, colour-coded by destination, with the exact name it carries in the joined .h5ad.
    Per-image-set columns are not here; the uns table covers those."""
    needle = filter_text.strip().lower()
    shown = [r for r in rows if not needle or _measurement_matches(r, needle)]
    summary = report_summary(rows)
    parts = ["<h3>Measurements</h3>"]
    labels = {DEST_X: "&rarr; X", DEST_OBS: "&rarr; obs", DEST_MERGED: "merged"}
    parts.append("<p>" + " &nbsp; ".join(
        f'<font color="{_DEST_STYLE[dest][2]}"><b>{summary[dest]} {labels[dest]}</b></font>'
        for dest in (DEST_X, DEST_OBS, DEST_MERGED)) + "</p>")
    if needle:
        parts.append(f"<p>Showing {len(shown)} of {len(rows)} measurements.</p>")
    parts.append('<table border="1" cellpadding="5" cellspacing="0">'
                 '<tr bgcolor="#dddddd"><th>Object</th><th>CellProfiler name</th>'
                 "<th>Defined in module</th><th>Category</th><th>Destination</th>"
                 "<th>Name in the joined .h5ad</th><th>Why</th></tr>")
    for r in shown:
        label, bg, fg = _DEST_STYLE[r.destination]
        module = f"{escape(r.module_name)} (#{r.module_num})" if r.module_num >= 0 else escape(r.module_name)
        name_cell = escape(r.anndata_name)
        reason_cell = f'<font color="#555555"><i>{escape(r.reason)}</i></font>' if r.reason else ""
        parts.append(
            f"<tr><td>{escape(r.object)}</td><td>{escape(r.cp_name)}</td><td>{module}</td>"
            f"<td>{escape(r.category)}</td>"
            f'<td bgcolor="{bg}"><font color="{fg}"><b>{label}</b></font></td>'
            f"<td>{name_cell}</td><td>{reason_cell}</td></tr>")
    parts.append("</table>")
    return "".join(parts)


def _uns_matches(row: UnsReportRow, needle: str) -> bool:
    return needle in " ".join([row.key, row.holds, row.detail]).lower()


def render_uns_table(rows: List[UnsReportRow], filter_text: str = "") -> str:
    """An HTML table (no <html>/<body> wrapper) of what the export carries outside X/obs/var, so a
    measurement missing from the tables above can be tracked down instead of assumed lost."""
    needle = filter_text.strip().lower()
    shown = [r for r in rows if not needle or _uns_matches(r, needle)]
    parts = ["<h3>Also exported, in uns</h3>",
             "<p>Not per-cell columns, so not in the tables above. Nothing here is lost, it is "
             "keyed by image set or by run instead of by cell.</p>"]
    if not rows:
        return "".join(parts[:1]) + "<p>Nothing.</p>"
    parts.append('<table border="1" cellpadding="5" cellspacing="0">'
                 '<tr bgcolor="#dddddd"><th>Key under uns</th><th>What it holds</th>'
                 "<th>Size</th></tr>")
    for r in shown:
        parts.append(f"<tr><td><tt>{escape(r.key)}</tt></td><td>{escape(r.holds)}</td>"
                     f"<td>{escape(r.detail)}</td></tr>")
    parts.append("</table>")
    if needle and len(shown) < len(rows):
        parts.append(f"<p>Showing {len(shown)} of {len(rows)} keys.</p>")
    return "".join(parts)


def render_report_html(channel_rows: List[ChannelReportRow], object_rows: List[ObjectReportRow],
                       measurement_rows: List[MeasurementReportRow],
                       uns_rows: List[UnsReportRow] = (), filter_text: str = "") -> str:
    """The full preview page: the three per-cell tables (channels, objects, measurements), then
    what else the export carries in uns. One filter box covers all four."""
    body = (render_channel_table(channel_rows, filter_text) + render_object_table(object_rows, filter_text) +
           render_measurement_table(measurement_rows, filter_text) +
           render_uns_table(list(uns_rows), filter_text))
    return f'<html><body><font size="{_BODY_FONT_SIZE}">{body}</font></body></html>'
