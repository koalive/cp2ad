"""Report of where each object and each measurement lands in the AnnData export: which module
made it, and whether it ends up in X, obs, merged into an existing column, or dropped, and under
what exact name.

Powers ExportToAnnData's "Press button to see where each feature will land" dialog. Pure Python,
no GUI or CellProfiler-runtime dependency, so it stays testable outside the module.

The measurement report replays assemble.py's own is_numeric/is_extrinsic/to_cpm_names rules and its
per-object dedup, instead of re-deriving them, so it can't drift from what a run actually produces.
"""
from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, fields
from html import escape
from typing import Dict, List

from .assemble import _is_numeric
from .introspect import Context, file_loaded_objects
from .names import is_extrinsic, to_cpm_names

DEST_X, DEST_OBS, DEST_MERGED, DEST_DROPPED = "X", "obs", "merged", "dropped"

REASON_EXTRINSIC = "position, orientation, or identity/linkage, not this object's own biology"


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
    destination: str            # DEST_X | DEST_OBS | DEST_MERGED | DEST_DROPPED
    anndata_name: str           # "<object>__<name>" on the joined file; "" when dropped
    reason: str                 # empty for a plain X measurement; explains the others


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
        category=f.category, destination=destination, anndata_name=anndata_name, reason=reason)


def measurement_report(ctx: Context) -> List[MeasurementReportRow]:
    """One row per (object, measurement, cpm name) the export machinery actually touches, scoped to
    the objects ctx.roles names, i.e. exactly what reaches the joined .h5ad. anndata_name is that
    joined file's name (<object>__<name>); a per-object file, if the module also writes one, uses
    the same name without the object prefix."""
    rows: List[MeasurementReportRow] = []
    role_objects = list(dict.fromkeys(ctx.roles.values()))  # role order, de-duplicated
    for obj in role_objects:
        claimed_by: Dict[tuple, str] = {}  # (is_extrinsic, cpm name) -> claiming cp_name
        for f in ctx.features:
            if f.object != obj:
                continue
            if not _is_numeric(f.coltype):
                rows.append(_measurement_row(f, DEST_DROPPED, "", f"column type is {f.coltype!r}, not numeric"))
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
    out = {DEST_X: 0, DEST_OBS: 0, DEST_MERGED: 0, DEST_DROPPED: 0}
    for r in rows:
        out[r.destination] += 1
    return out


def write_csv(rows, path: str) -> None:
    """Write a list of ObjectReportRow or MeasurementReportRow to a CSV file, one row per line,
    columns in field order. An empty list still writes a header-only file."""
    fieldnames = [f.name for f in fields(rows[0])] if rows else []
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(asdict(r))


# ---- HTML rendering (for the preview dialog's wx.html.HtmlWindow) ----
#
# wx.html.HtmlWindow is a lightweight HTML 3.2-ish renderer with only partial CSS support, and no
# external browser engine (unlike wx.html2.WebView, which needs one). Colour-coding here uses the
# classic <td bgcolor="..."> / <font color="..."> attributes it understands natively; a <style>
# block or inline CSS renders inconsistently across wx versions. Font size is set by the dialog
# itself, via HtmlWindow.SetStandardFonts, so it follows the user's system font size instead of a
# fixed size baked into this HTML.

_DEST_STYLE = {
    DEST_X: ("X (morphology)", "#e6f4ea", "#1e7e34"),
    DEST_OBS: ("obs (extrinsic)", "#e8f0fe", "#1a56c4"),
    DEST_MERGED: ("obs (merged)", "#fff6e0", "#8a6d1a"),
    DEST_DROPPED: ("dropped", "#f1f1f1", "#777777"),
}


def _measurement_matches(row: MeasurementReportRow, needle: str) -> bool:
    haystack = " ".join([row.object, row.cp_name, row.module_name, row.category,
                         row.destination, row.anndata_name, row.reason]).lower()
    return needle in haystack


def _object_matches(row: ObjectReportRow, needle: str) -> bool:
    haystack = " ".join([row.role, row.object, row.module_name, row.source]).lower()
    return needle in haystack


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
    parts.append('<table border="1" cellpadding="4" cellspacing="0">'
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
    """An HTML table (no <html>/<body> wrapper) of every measurement, colour-coded by destination,
    with the exact name it carries in the joined .h5ad or the reason it gets none."""
    needle = filter_text.strip().lower()
    shown = [r for r in rows if not needle or _measurement_matches(r, needle)]
    summary = report_summary(rows)
    parts = ["<h3>Measurements</h3>"]
    labels = {DEST_X: "&rarr; X", DEST_OBS: "&rarr; obs", DEST_MERGED: "merged", DEST_DROPPED: "dropped"}
    parts.append("<p>" + " &nbsp; ".join(
        f'<font color="{_DEST_STYLE[dest][2]}"><b>{summary[dest]} {labels[dest]}</b></font>'
        for dest in (DEST_X, DEST_OBS, DEST_MERGED, DEST_DROPPED)) + "</p>")
    if needle:
        parts.append(f"<p>Showing {len(shown)} of {len(rows)} measurements.</p>")
    parts.append('<table border="1" cellpadding="3" cellspacing="0">'
                 '<tr bgcolor="#dddddd"><th>Object</th><th>CellProfiler name</th>'
                 "<th>Defined in module</th><th>Category</th><th>Destination</th>"
                 "<th>Name in the joined .h5ad</th><th>Why</th></tr>")
    for r in shown:
        label, bg, fg = _DEST_STYLE[r.destination]
        module = f"{escape(r.module_name)} (#{r.module_num})" if r.module_num >= 0 else escape(r.module_name)
        name_cell = escape(r.anndata_name) if r.anndata_name else "<i>not exported</i>"
        reason_cell = f'<font color="#555555"><i>{escape(r.reason)}</i></font>' if r.reason else ""
        parts.append(
            f"<tr><td>{escape(r.object)}</td><td>{escape(r.cp_name)}</td><td>{module}</td>"
            f"<td>{escape(r.category)}</td>"
            f'<td bgcolor="{bg}"><font color="{fg}"><b>{label}</b></font></td>'
            f"<td>{name_cell}</td><td>{reason_cell}</td></tr>")
    parts.append("</table>")
    return "".join(parts)


def render_report_html(object_rows: List[ObjectReportRow], measurement_rows: List[MeasurementReportRow],
                       filter_text: str = "") -> str:
    """The full preview page: the objects table above the measurements table, one filter box
    covering both."""
    body = render_object_table(object_rows, filter_text) + render_measurement_table(measurement_rows, filter_text)
    return f"<html><body>{body}</body></html>"
