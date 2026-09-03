"""Feature-by-feature report of where each CellProfiler measurement lands in the AnnData export --
X, obs, merged into an existing obs/X column, or dropped -- and under what exact name.

Powers ExportToAnnData's "Show features and where they will land" preview dialog. Pure Python, no
GUI or CellProfiler-runtime dependency (like the rest of cpexport), so it is unit-testable and
reusable outside the module (e.g. a --dry-run report).

This mirrors assemble.py's own decisions by construction, not by re-deriving them: it calls the
exact same is_numeric/is_extrinsic/to_cpm_names rules assemble._var_columns and
assemble._extrinsic_columns use, and replays their per-object, per-branch "seen name" dedup, so the
report can never drift from what a run actually produces.
"""
from __future__ import annotations

from dataclasses import dataclass
from html import escape
from typing import Dict, List, Tuple

from .assemble import _is_numeric
from .introspect import Context
from .names import Feature, is_extrinsic, to_cpm_names

DEST_X, DEST_OBS, DEST_MERGED, DEST_DROPPED = "X", "obs", "merged", "dropped"

REASON_EXTRINSIC = "position, orientation, or identity/linkage -- not this object's own biology"


@dataclass(frozen=True)
class FeatureReportRow:
    object: str
    cp_name: str
    module_name: str
    module_num: int
    category: str
    destination: str            # DEST_X | DEST_OBS | DEST_MERGED | DEST_DROPPED
    anndata_name: str           # "<object>__<name>" on the joined file; "" when dropped
    reason: str                 # empty for a plain X feature; explains obs/merged/dropped otherwise


def _row(f: Feature, destination: str, anndata_name: str, reason: str) -> FeatureReportRow:
    return FeatureReportRow(
        object=f.object, cp_name=f.cp_name,
        module_name=f.module_name or "?", module_num=-1 if f.module_num is None else f.module_num,
        category=f.category, destination=destination, anndata_name=anndata_name, reason=reason)


def feature_report(ctx: Context) -> List[FeatureReportRow]:
    """One row per (object, feature, cpm-name) the export machinery actually touches, restricted to
    the objects ctx.roles names -- exactly the set that reaches the joined .h5ad. `anndata_name` is
    that joined file's name (`<object>__<name>`); a per-object file, if the module also writes one,
    uses the same name without the object prefix."""
    rows: List[FeatureReportRow] = []
    role_objects = list(dict.fromkeys(ctx.roles.values()))  # role order, de-duplicated
    for obj in role_objects:
        claimed_by: Dict[Tuple[bool, str], str] = {}  # (is_extrinsic, cpm name) -> claiming cp_name
        for f in ctx.features:
            if f.object != obj:
                continue
            if not _is_numeric(f.coltype):
                rows.append(_row(f, DEST_DROPPED, "", f"column type is {f.coltype!r}, not numeric"))
                continue
            extrinsic = is_extrinsic(f)
            for name, backend in to_cpm_names(f, ctx.channels):
                key = (extrinsic, name)
                anndata_name = f"{obj}__{name}"
                if key in claimed_by:
                    rows.append(_row(f, DEST_MERGED, anndata_name,
                                     f"same value as {claimed_by[key]} under this name"))
                    continue
                claimed_by[key] = f.cp_name
                rows.append(_row(f, DEST_OBS if extrinsic else DEST_X, anndata_name,
                                 REASON_EXTRINSIC if extrinsic else ""))
    return rows


def report_summary(rows: List[FeatureReportRow]) -> Dict[str, int]:
    """Counts per destination, for a one-line summary above the detail table."""
    out = {DEST_X: 0, DEST_OBS: 0, DEST_MERGED: 0, DEST_DROPPED: 0}
    for r in rows:
        out[r.destination] += 1
    return out


# ---- HTML rendering (for the preview dialog's wx.html.HtmlWindow) ----
#
# wx.html.HtmlWindow is a lightweight, dependency-free HTML 3.2-ish renderer (no external browser
# engine required, unlike wx.html2.WebView) with only partial CSS support. Colour-coding therefore
# uses the classic <td bgcolor="..."> / <font color="..."> attributes it is built to understand,
# rather than relying on a <style> block or inline "style=" CSS, which render inconsistently
# across wx versions.

_DEST_STYLE = {
    DEST_X: ("X (morphology)", "#e6f4ea", "#1e7e34"),
    DEST_OBS: ("obs (extrinsic)", "#e8f0fe", "#1a56c4"),
    DEST_MERGED: ("obs (merged)", "#fff6e0", "#8a6d1a"),
    DEST_DROPPED: ("dropped", "#f1f1f1", "#777777"),
}


def _matches(row: FeatureReportRow, needle: str) -> bool:
    haystack = " ".join([row.object, row.cp_name, row.module_name, row.category,
                         row.destination, row.anndata_name, row.reason]).lower()
    return needle in haystack


def _summary_line(summary: Dict[str, int], n_objects: int, n_total: int) -> str:
    parts = [f"<b>{n_total} features</b> across {n_objects} object(s): "]
    labels = {DEST_X: "&rarr; X", DEST_OBS: "&rarr; obs", DEST_MERGED: "merged", DEST_DROPPED: "dropped"}
    parts.append(" &nbsp; ".join(
        f'<font color="{_DEST_STYLE[dest][2]}"><b>{summary[dest]} {labels[dest]}</b></font>'
        for dest in (DEST_X, DEST_OBS, DEST_MERGED, DEST_DROPPED)))
    return "".join(parts)


def render_html(rows: List[FeatureReportRow], filter_text: str = "") -> str:
    """A self-contained HTML table for wx.html.HtmlWindow: one row per feature, colour-coded by
    destination, with the exact name it will carry in the joined .h5ad (or the reason it does not
    get one). Filters case-insensitively across every column when `filter_text` is non-empty."""
    needle = filter_text.strip().lower()
    shown = [r for r in rows if not needle or _matches(r, needle)]
    summary = report_summary(rows)
    parts = ['<html><body><font face="sans-serif" size="2">']
    parts.append(f"<p>{_summary_line(summary, len({r.object for r in rows}), len(rows))}</p>")
    if needle:
        parts.append(f"<p>Showing {len(shown)} of {len(rows)} matching &quot;{escape(filter_text)}&quot;.</p>")
    parts.append('<table border="1" cellpadding="3" cellspacing="0">'
                 "<tr bgcolor=\"#dddddd\"><th>Object</th><th>CellProfiler name</th>"
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
    parts.append("</table></font></body></html>")
    return "".join(parts)
