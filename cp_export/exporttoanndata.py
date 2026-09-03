"""
ExportToAnnData
===============

**ExportToAnnData** exports per-cell measurements as an AnnData ``.h5ad`` file whose feature
names, ``obs`` and ``uns`` keys follow ``squidpy.experimental.im.calculate_image_features``
(cp_measure naming), with one row per cell (the primary, secondary and tertiary objects joined --
in a standard pipeline nuclei, cells and cytoplasm) and the full pipeline provenance in
``uns["cellprofiler"]``.

Extrinsic measurements never enter ``X``: an object's own position/orientation (its center, its
per-channel intensity-weighted center, the corners of its bounding box, its angle relative to the
image axes) and its identity/linkage (its own arbitrary object number, ``Parent_*``/
``Neighbors_*ClosestObjectNumber`` references to another object's label, ``Children_*_Count``)
describe where or how an object was imaged, or which arbitrary label CellProfiler assigned it or
its neighbors -- not its biology -- and would bias downstream similarity/clustering accordingly.
They are exported to ``obs`` instead (per object on the compartment's own file,
``<Object>__<name>`` on the joined file); see ``cpexport.names.is_extrinsic`` for the exact rules.

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

LOGGER = logging.getLogger(__name__)
POLICY_CHOICES = tuple(p.capitalize() for p in POLICIES)  # single source of truth: cpexport.assemble.POLICIES
ROLE_AUTO, ROLE_MANUAL = "Automatic", "Manual"


class ExportToAnnData(Module):
    module_name = "ExportToAnnData"
    category = ["File Processing", "Data Tools"]
    variable_revision_number = 1

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
        self.primary_object = LabelSubscriber("Primary objects (e.g. nuclei)", "None")
        self.secondary_object = LabelSubscriber("Secondary objects (e.g. cells)", "None")
        self.tertiary_object = LabelSubscriber("Tertiary objects (e.g. cytoplasm)", "None")
        self.wants_overwrite = Binary("Overwrite existing files without warning?", True)

    def settings(self):
        return [self.directory, self.prefix, self.wants_per_object, self.policy, self.role_mode,
                self.primary_object, self.secondary_object, self.tertiary_object, self.wants_overwrite]

    def visible_settings(self):
        result = [self.directory, self.prefix, self.wants_per_object, self.policy, self.role_mode,
                  self.autoconfig]
        if self.role_mode == ROLE_MANUAL:
            result += [self.primary_object, self.secondary_object, self.tertiary_object]
        return result + [self.wants_overwrite]

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
                lines.append('%s = "%s" -- made by %s (module #%d)' % (role, name, info.module_name, info.module_num))
            else:
                lines.append("%s = None -- no %s object in this pipeline, so none is exported" % (role, role))
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
            if obj == "__joined__":
                continue
            t = tables[obj]
            t.uns["cellprofiler"] = prov
            write_h5ad(t, path)
            LOGGER.info("ExportToAnnData: wrote %s", path)

    def upgrade_settings(self, setting_values, variable_revision_number, module_name):
        return setting_values, variable_revision_number

    def volumetric(self):
        return False
