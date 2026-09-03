"""
ExportToAnnData
===============

**ExportToAnnData** exports per-cell measurements as an AnnData ``.h5ad`` file whose feature
names, ``obs`` and ``uns`` keys follow ``squidpy.experimental.im.calculate_image_features``
(cp_measure naming). Each row is one cell: the primary, secondary and tertiary objects joined
together, nuclei, cells and cytoplasm in a standard pipeline, with the full pipeline provenance in
``uns["cellprofiler"]``.

Extrinsic measurements never enter ``X``. An object's own position and orientation (its center, its
per-channel intensity-weighted center, the corners of its bounding box, its angle relative to the
image axes) and its identity and linkage (its own arbitrary object number, ``Parent_*``/
``Neighbors_*ClosestObjectNumber`` references to another object's label, ``Children_*_Count``)
describe where or how an object was imaged, or which arbitrary label CellProfiler assigned it or
its neighbors, not its biology. Left in ``X`` they would bias downstream similarity and clustering
on those instead. They are exported to ``obs`` instead, per object on the compartment's own file
and as ``<Object>__<name>`` on the joined file; see ``cpexport.names.is_extrinsic`` for the exact
rules, or press "Press button to see where each object and measurement will land" in the module's
own settings to inspect them for your pipeline directly.

Place this module **last** in the pipeline. It reads the pipeline structure itself: channels,
objects, which module produced them, and how measurements are named. Nothing has to be
configured for a standard IdentifyPrimaryObjects -> IdentifySecondaryObjects ->
IdentifyTertiaryObjects pipeline.

|
============ ============ ===============
Supports 2D? Supports 3D? Respects masks?
============ ============ ===============
YES          NO           YES
============ ============ ===============
"""
import logging
import os
import sys

# Dev loop: CellProfiler's "Test > Reload Modules' Source" re-imports this file but leaves the
# cpexport package cached in sys.modules. With CPEXPORT_DEV=1 we purge it here so a reload picks
# up library edits too (settings changes still need a CellProfiler restart).
if os.environ.get("CPEXPORT_DEV"):
    for _name in [m for m in sys.modules if m == "cpexport" or m.startswith("cpexport.")]:
        del sys.modules[_name]

from cellprofiler_core.constants.module import IO_FOLDER_CHOICE_HELP_TEXT
from cellprofiler_core.module import Module
from cellprofiler_core.preferences import (ABSOLUTE_FOLDER_NAME, DEFAULT_INPUT_FOLDER_NAME,
                                           DEFAULT_INPUT_SUBFOLDER_NAME, DEFAULT_OUTPUT_FOLDER_NAME,
                                           DEFAULT_OUTPUT_SUBFOLDER_NAME, get_headless)
from cellprofiler_core.setting import Binary, ValidationError
from cellprofiler_core.setting.choice import Choice
from cellprofiler_core.setting.do_something import DoSomething
from cellprofiler_core.setting.subscriber import LabelSubscriber
from cellprofiler_core.setting.text import Directory, Text

from cpexport.advice import advice as _advice
from cpexport.assemble import POLICIES, build_object_table, join_tables, provenance
from cpexport.h5ad import write_h5ad
from cpexport.introspect import RoleError, build_context
from cpexport.preview import measurement_report, object_report, render_report_html, write_csv

LOGGER = logging.getLogger(__name__)
POLICY_CHOICES = tuple(p.capitalize() for p in POLICIES)  # single source of truth: cpexport.assemble.POLICIES
ROLE_AUTO, ROLE_MANUAL = "Automatic", "Manual"


class ExportToAnnData(Module):
    module_name = "ExportToAnnData"
    category = ["File Processing", "Data Tools"]
    variable_revision_number = 2

    def create_settings(self):
        self.directory = Directory(
            "Output file location",
            dir_choices=[DEFAULT_OUTPUT_FOLDER_NAME, DEFAULT_OUTPUT_SUBFOLDER_NAME, ABSOLUTE_FOLDER_NAME,
                         DEFAULT_INPUT_FOLDER_NAME, DEFAULT_INPUT_SUBFOLDER_NAME],
            doc=IO_FOLDER_CHOICE_HELP_TEXT)
        self.directory.dir_choice = DEFAULT_OUTPUT_FOLDER_NAME
        self.prefix = Text("File name prefix", "cellprofiler",
                           doc="Files are written as <prefix>.h5ad and, optionally, <prefix>_<Object>.h5ad.")
        self.wants_per_object = Binary(
            "Also write one file per object?", False,
            doc="Write squidpy-identical per-object tables (unprefixed feature names) next to the joined file.")
        self.policy = Choice(
            "When an object is not exactly one primary + one secondary + one tertiary", POLICY_CHOICES,
            doc="Flag: keep the row and mark obs['qc_flag']. Drop: remove it. Error: abort the run.")
        self.role_mode = Choice(
            "How to pick primary / secondary / tertiary objects", (ROLE_AUTO, ROLE_MANUAL),
            doc="""Automatic follows IdentifyPrimary/Secondary/TertiaryObjects and FilterObjects; with several
candidate primary objects and no secondary/tertiary chain it picks the one the most other modules use as an
input. Manual exports exactly the roles you set: a role left as "None" is simply not exported.""")
        self.autoconfig = DoSomething(
            "", "Auto-configure from this pipeline", self.do_autoconfig,
            doc="""Detect the primary/secondary/tertiary objects from the modules above, write them into the
three role settings (switching to Manual so later pipeline edits cannot silently change the export), and
explain every choice. Warnings (missing plate metadata, ambiguous objects, ...) are shown in the same dialog.""")
        self.preview = DoSomething(
            "", "Press button to see where each object and measurement will land", self.do_preview,
            doc="""Opens two tables for the objects currently configured above (Automatic or Manual):
which module made each object, and where each of its measurements ends up in the exported
``.h5ad``. Use it before running the pipeline to check a name before you go looking for it.

The objects table has one row per role (primary, secondary, tertiary): the object filling that
role, the module that produced it, and whether that module identified it from an image (source
"pipeline") or loaded it from a label file (source "file").

The measurements table has one row per (object, measurement). **CellProfiler name** is the
measurement as CellProfiler itself names it, e.g. ``AreaShape_Area``; **defined in module** names
the module (and module number) that computed it. **Destination** is where it lands in the AnnData
file:

-  *X (morphology)*: a real measurement of the object's own shape, intensity or texture. It becomes
   a column of ``X``/``var``, named ``<Object>__<name>`` on the joined file, so it contributes to
   any similarity or clustering computed from ``X``.
-  *obs (extrinsic)*: position, orientation, or an identifier or link to another row: an object's
   own arbitrary number, a ``Parent_``/``Children_`` reference, a neighbor's object number. These
   describe where or how an object was imaged, or an arbitrary CellProfiler label, not its biology,
   so they stay out of ``X`` and land in ``obs`` instead, also named ``<Object>__<name>``. The *Why*
   column gives the specific reason.
-  *obs (merged)*: a second CellProfiler measurement with the exact same value as one that already
   claimed this name (``AreaShape_Center_X`` and ``Location_Center_X`` are the same coordinate, for
   example). Its value lives under the name shown, not a separate column.
-  *dropped*: not a numeric column, e.g. a file name or path, so it cannot enter ``X`` at all and
   carries no information into the AnnData file.

**Name in the joined .h5ad** is the exact column (for *X*) or ``obs`` key (for *obs*) the value
will have. A per-object file, if "Also write one file per object?" is enabled, uses this same name
without the ``<Object>__`` prefix.

Type in the filter box to narrow both tables by object, measurement name, module, category,
destination, or exported name. The tables reflect the role settings as they are right now; change
Automatic vs Manual or the role objects and press the button again to refresh them. Check "Also
write the object/measurement tables as CSV?" below to save this same information next to the
``.h5ad`` on every run, not just when you press this button.""")
        self.primary_object = LabelSubscriber("Primary objects (e.g. nuclei)", "None")
        self.secondary_object = LabelSubscriber("Secondary objects (e.g. cells)", "None")
        self.tertiary_object = LabelSubscriber("Tertiary objects (e.g. cytoplasm)", "None")
        self.wants_feature_tables = Binary(
            "Also write the object/measurement tables as CSV?", False,
            doc="""Select "Yes" to write ``<prefix>_objects.csv`` and ``<prefix>_measurements.csv``
next to the ``.h5ad`` on every run: the same two tables "Press button to see where each object and
measurement will land" shows, as plain CSV, so provenance travels with the export without opening
CellProfiler again.""")
        self.wants_overwrite = Binary("Overwrite existing files without warning?", True)

    def settings(self):
        return [self.directory, self.prefix, self.wants_per_object, self.policy, self.role_mode,
                self.primary_object, self.secondary_object, self.tertiary_object,
                self.wants_feature_tables, self.wants_overwrite]

    def visible_settings(self):
        result = [self.directory, self.prefix, self.wants_per_object, self.policy, self.role_mode,
                  self.autoconfig]
        if self.role_mode == ROLE_MANUAL:
            result += [self.primary_object, self.secondary_object, self.tertiary_object]
        return result + [self.preview, self.wants_feature_tables, self.wants_overwrite]

    # ---- auto-configuration button -------------------------------------------------------
    _gui_pipeline = None

    def on_activated(self, workspace):
        # the GUI calls this when the module is selected; it is the only place a module can grab
        # the pipeline for use in a DoSomething callback (same pattern as the Metadata module)
        self._gui_pipeline = workspace.pipeline

    def on_deactivated(self):
        self._gui_pipeline = None

    def apply_autoconfig(self, pipeline):
        """Detect roles afresh (ignoring any Manual values), pin them as Manual, and return the
        explanation text. Raises RoleError when the pipeline identifies no objects at all."""
        ctx = build_context(pipeline, roles=None)
        lines = []
        for role, setting in (("primary", self.primary_object), ("secondary", self.secondary_object),
                              ("tertiary", self.tertiary_object)):
            name = ctx.roles.get(role)
            setting.value = name if name else "None"
            if name:
                info = ctx.objects[name]
                lines.append('%s = "%s", made by %s (module #%d)' % (role, name, info.module_name, info.module_num))
            else:
                lines.append("%s = None. No %s object exists in this pipeline, so nothing is exported "
                             "for that role." % (role, role))
        self.role_mode.value = ROLE_MANUAL
        text = ("These settings were filled in from the pipeline:\n\n- " + "\n- ".join(lines) +
                "\n\nThe mode was switched to Manual so that later pipeline edits cannot silently "
                "change which objects are exported.")
        fallback = (ctx.role_note or {}).get("fallback")
        if fallback:
            text += ("\n\nNote: the choice used the '%s' fallback between the candidates %s." %
                     (fallback, ", ".join(ctx.role_note.get("candidates", []))))
        warnings = _advice(ctx)
        if warnings:
            text += "\n\nWarnings:\n\n- " + "\n- ".join(warnings)
        return text, bool(warnings)

    def do_autoconfig(self):
        import wx
        if self._gui_pipeline is None:
            wx.MessageBox("The pipeline is not available yet - click another module and then this one "
                          "again, then retry.", caption="ExportToAnnData", style=wx.OK | wx.ICON_WARNING)
            return
        try:
            text, warned = self.apply_autoconfig(self._gui_pipeline)
        except RoleError as e:
            wx.MessageBox("Nothing was changed: %s\n\nAdd an IdentifyPrimaryObjects (or similar) module "
                          "before this one, or set the roles by hand in Manual mode." % e,
                          caption="ExportToAnnData auto-configuration", style=wx.OK | wx.ICON_ERROR)
            return
        wx.MessageBox(text, caption="ExportToAnnData auto-configuration",
                      style=wx.OK | (wx.ICON_WARNING if warned else wx.ICON_INFORMATION))

    # ---- object/measurement preview button --------------------------------------------------
    def do_preview(self):
        import wx
        if self._gui_pipeline is None:
            wx.MessageBox("The pipeline is not available yet - click another module and then this one "
                          "again, then retry.", caption="ExportToAnnData", style=wx.OK | wx.ICON_WARNING)
            return
        try:
            ctx = build_context(self._gui_pipeline, roles=self._roles())
        except RoleError as e:
            wx.MessageBox("Nothing to preview: %s\n\nAdd an IdentifyPrimaryObjects (or similar) module "
                          "before this one, or set the roles by hand in Manual mode." % e,
                          caption="ExportToAnnData feature preview", style=wx.OK | wx.ICON_ERROR)
            return
        _show_feature_preview(object_report(ctx), measurement_report(ctx))

    # ---- validation / advice -------------------------------------------------------------
    def _roles(self):
        if self.role_mode != ROLE_MANUAL:
            return None
        roles = {}
        for role, setting in (("primary", self.primary_object), ("secondary", self.secondary_object),
                              ("tertiary", self.tertiary_object)):
            if setting.value not in (None, "", "None"):
                roles[role] = setting.value
        return roles or None

    def validate_module(self, pipeline):
        if pipeline.modules()[-1] is not self:
            raise ValidationError("ExportToAnnData must be the last module; measurements made after it "
                                  "are not exported.", self.prefix)
        try:
            build_context(pipeline, roles=self._roles())
        except RoleError as e:
            raise ValidationError(str(e), self.role_mode)

    def validate_module_warnings(self, pipeline):
        try:
            ctx = build_context(pipeline, roles=self._roles())
        except RoleError:
            return  # a hard failure; validate_module already reports it as a ValidationError
        messages = _advice(ctx)
        if messages:
            raise ValidationError("\n\n".join(messages), self.prefix)

    # ---- run ---------------------------------------------------------------------------
    def _paths(self, ctx):
        base = self.directory.get_absolute_path()
        paths = {"__joined__": os.path.join(base, f"{self.prefix.value}.h5ad")}
        if self.wants_per_object.value:
            for obj in ctx.roles.values():
                paths[obj] = os.path.join(base, f"{self.prefix.value}_{obj}.h5ad")
        if self.wants_feature_tables.value:
            paths["__objects_csv__"] = os.path.join(base, f"{self.prefix.value}_objects.csv")
            paths["__measurements_csv__"] = os.path.join(base, f"{self.prefix.value}_measurements.csv")
        return paths

    def prepare_run(self, workspace):
        try:
            ctx = build_context(workspace.pipeline, roles=self._roles())
        except RoleError as e:  # a clean refusal: CellProfiler aborts the run without a traceback
            LOGGER.error("ExportToAnnData: %s", e)
            return False
        existing = [p for p in self._paths(ctx).values() if os.path.exists(p)]
        if existing and not self.wants_overwrite.value:
            if get_headless():
                raise ValueError("ExportToAnnData: output exists and overwrite is off: " + ", ".join(existing))
            import wx
            if wx.MessageBox("Overwrite %s?" % ", ".join(existing), caption="ExportToAnnData",
                             style=wx.YES_NO | wx.ICON_QUESTION) != wx.YES:
                return False
        return True

    def run(self, workspace):
        pass  # everything happens in post_run, like ExportToSpreadsheet

    def display(self, workspace, figure):
        figure.set_subplots((1, 1))
        figure.subplot_table(0, 0, [["Output", "written in post_run"]])

    def post_run(self, workspace):
        if workspace.pipeline.test_mode:
            return
        m = workspace.measurements
        try:  # the analysis worker may not have run prepare_run in this process
            ctx = build_context(workspace.pipeline, roles=self._roles())
        except RoleError as e:
            LOGGER.error("ExportToAnnData: %s", e)
            return
        paths = self._paths(ctx)
        os.makedirs(os.path.dirname(paths["__joined__"]), exist_ok=True)
        settings = {s.text: s.value_text for s in self.settings()}
        prov = provenance(ctx, m, settings)
        tables = {obj: build_object_table(ctx, m, obj) for obj in ctx.roles.values()}
        joined = join_tables(ctx, m, tables, policy=self.policy.value.lower())
        joined.uns["cellprofiler"] = prov
        for message in _advice(ctx):
            LOGGER.warning("ExportToAnnData: %s", message)
        LOGGER.info("ExportToAnnData: writing %d cells x %d features (~%.1f MB of float32) to %s",
                    joined.X.shape[0], joined.X.shape[1], joined.X.nbytes / 1e6, paths["__joined__"])
        write_h5ad(joined, paths["__joined__"])
        LOGGER.info("ExportToAnnData: wrote %s (%d cells x %d features)", paths["__joined__"], *joined.X.shape)
        for obj, path in paths.items():
            if obj in ("__joined__", "__objects_csv__", "__measurements_csv__"):
                continue
            t = tables[obj]
            t.uns["cellprofiler"] = prov
            write_h5ad(t, path)
            LOGGER.info("ExportToAnnData: wrote %s", path)
        if "__objects_csv__" in paths:
            write_csv(object_report(ctx), paths["__objects_csv__"])
            LOGGER.info("ExportToAnnData: wrote %s", paths["__objects_csv__"])
        if "__measurements_csv__" in paths:
            write_csv(measurement_report(ctx), paths["__measurements_csv__"])
            LOGGER.info("ExportToAnnData: wrote %s", paths["__measurements_csv__"])

    def upgrade_settings(self, setting_values, variable_revision_number, module_name):
        if variable_revision_number == 1:
            # v1 has no "write CSV provenance tables" setting. Default it to No, insert it right
            # before "Overwrite existing files without warning?" (the position settings() puts it
            # at now), and leave every other value untouched, so pipelines saved before this
            # feature existed keep loading.
            setting_values = setting_values[:8] + ["No"] + setting_values[8:]
            variable_revision_number = 2
        return setting_values, variable_revision_number

    def volumetric(self):
        return False


def _show_feature_preview(object_rows, measurement_rows):
    """The "Press button to see where each object and measurement will land" dialog: a filterable
    pair of tables (cpexport.preview.render_report_html), one for object provenance and one for
    measurement destinations, for whatever objects are currently configured.

    wx is imported here, not at module scope, so the module stays importable headless when wx
    isn't installed, the same reason do_autoconfig/prepare_run import it locally instead of at the
    top of the file.
    """
    import wx
    import wx.html

    class _FeaturePreviewDialog(wx.Dialog):
        def __init__(self, object_rows, measurement_rows):
            super().__init__(None, title="ExportToAnnData feature preview", size=(1100, 700),
                             style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER | wx.CLOSE_BOX)
            self.object_rows = object_rows
            self.measurement_rows = measurement_rows
            sizer = wx.BoxSizer(wx.VERTICAL)

            filter_row = wx.BoxSizer(wx.HORIZONTAL)
            filter_row.Add(wx.StaticText(self, label="Filter (object, measurement, module, destination):"),
                           0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
            self.filter_ctrl = wx.TextCtrl(self)
            self.filter_ctrl.Bind(wx.EVT_TEXT, self.on_filter)
            filter_row.Add(self.filter_ctrl, 1, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
            sizer.Add(filter_row, 0, wx.EXPAND)

            self.html_window = wx.html.HtmlWindow(self, style=wx.html.HW_SCROLLBAR_AUTO)
            # HtmlWindow's own default is legible but small for a dense table; follow the OS GUI
            # font size instead of a size baked into the HTML, so it scales with the user's own
            # font/accessibility settings rather than sitting fixed regardless of them.
            base_size = wx.SystemSettings.GetFont(wx.SYS_DEFAULT_GUI_FONT).GetPointSize()
            self.html_window.SetStandardFonts(size=base_size)
            self.html_window.SetPage(render_report_html(self.object_rows, self.measurement_rows))
            sizer.Add(self.html_window, 1, wx.EXPAND | wx.ALL, 5)

            close_button = wx.Button(self, id=wx.ID_OK, label="Close")
            close_button.Bind(wx.EVT_BUTTON, lambda evt: self.EndModal(wx.ID_OK))
            close_button.SetDefault()
            sizer.Add(close_button, 0, wx.ALIGN_RIGHT | wx.ALL, 8)

            self.SetSizer(sizer)

        def on_filter(self, evt):
            self.html_window.SetPage(
                render_report_html(self.object_rows, self.measurement_rows, self.filter_ctrl.GetValue()))

    dlg = _FeaturePreviewDialog(object_rows, measurement_rows)
    try:
        dlg.ShowModal()
    finally:
        dlg.Destroy()
