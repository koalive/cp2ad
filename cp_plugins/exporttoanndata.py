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
and as ``<Object>__<name>`` on the joined file; see ``scverse_export.names.is_extrinsic`` for the exact
rules, or press "See where each measurement will land" in the module's own settings to inspect them
for your pipeline directly.

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
# scverse_export package cached in sys.modules. With SCVERSE_EXPORT_DEV=1 we purge it here so a reload picks
# up library edits too (settings changes still need a CellProfiler restart).
if os.environ.get("SCVERSE_EXPORT_DEV"):
    for _name in [m for m in sys.modules if m == "scverse_export" or m.startswith("scverse_export.")]:
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

from scverse_export.advice import advice as _advice
from scverse_export.assemble import POLICIES
from scverse_export.export import build_export
from scverse_export.h5ad import write_h5ad
from scverse_export.introspect import RoleError, build_context
from scverse_export.preview import (channel_report, measurement_report, object_report,
                                    render_report_html, row_name_preview, uns_report)
from scverse_export.samples import PREFIX as TAG_PREFIX, detect_sample_tags, parse_tags, qualify

LOGGER = logging.getLogger(__name__)
POLICY_CHOICES = tuple(p.capitalize() for p in POLICIES)  # single source of truth: scverse_export.assemble.POLICIES
ROLE_AUTO, ROLE_MANUAL = "Automatic", "Manual"


class ExportToAnnData(Module):
    module_name = "ExportToAnnData"
    category = ["File Processing", "Data Tools"]
    variable_revision_number = 4

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
        self.role_mode = Choice(
            "How to pick primary / secondary / tertiary objects", (ROLE_AUTO, ROLE_MANUAL),
            doc="""Automatic follows IdentifyPrimary/Secondary/TertiaryObjects and FilterObjects; with several
candidate primary objects and no secondary/tertiary chain it picks the one the most other modules use as an
input. Manual exports exactly the roles you set: a role left as "None" is simply not exported.""")
        self.sample_key_mode = Choice(
            "How to build row names", (ROLE_AUTO, ROLE_MANUAL),
            doc="""Every row is named `<sample key>_<label>`, where the sample key identifies the field
of view the object came from and the label is the integer CellProfiler gave it there.

**Automatic** reads the sample key off the pipeline's Metadata tags, taking at most one plate part
(Plate, Barcode, PlateID), one well part (Well, or Row plus Column together), and one site part
(Site, Field, FieldIndex, Position). A Plate/Well/Site pipeline gets ``P1_A01_1_5``; an Opera or
Harmony pipeline with Well and Field gets ``A02_03_5``. Frame, Series and Channel are never used,
because they index a z plane, a timepoint, or a channel inside one field of view rather than the
field itself.

Whichever tags are used, they have to give every image set a different key, or rows from two
fields of view would collide. That is checked against the real values during the run, and
``img<n>`` is appended when it fails, so a key is never ambiguous. With no usable tags at all the
key is ``img<n>`` alone, which is unique within one run but not across runs.

**Manual** uses exactly the tags you list below, which makes the scheme stable. Automatic depends
on what a run contains: if a tag turns out not to separate the image sets, that run appends
``img<n>`` and the one before it may not have, so the same pipeline over different subsets can
produce differently shaped names. Naming the tags yourself fixes them, and is worth doing before
exporting anything you intend to concatenate.

"Auto-configure from this pipeline" fills these in and switches to Manual, which is the quickest
way to get a stable scheme that matches the pipeline. The resolved scheme is logged at the start
of the export, recorded in ``uns["cellprofiler"]["sample_naming"]``, and shown at the top of "See
where each measurement will land".""")
        self.sample_key_tags = Text(
            "Metadata tags for row names", "Plate,Well,Site",
            doc="""*(Used only when building row names Manually)*

Comma-separated Metadata tag names, in the order they should appear in the key. The ``Metadata_``
prefix is optional, so ``Well,Field`` and ``Metadata_Well,Metadata_Field`` mean the same thing.

A tag this pipeline does not have is reported and the key falls back to ``img<n>``, rather than
producing a key with a blank in it. Leave this empty to name rows by image number alone.""")
        self.primary_object = LabelSubscriber("Primary objects (e.g. nuclei)", "None")
        self.secondary_object = LabelSubscriber("Secondary objects (e.g. cells)", "None")
        self.tertiary_object = LabelSubscriber("Tertiary objects (e.g. cytoplasm)", "None")
        self.preview = DoSomething(
            "", "See where each measurement will land", self.do_preview,
            doc="""Opens three tables for the channels and objects currently configured above
(Automatic or Manual): which module made each channel and each object, and where each measurement
ends up in the exported ``.h5ad``. Use it before running the pipeline to check a name before you go
looking for it.

The first three tables cover the per-cell part of the export, and show only what actually reaches
it. A channel loaded or computed but never fed to a measurement (a raw image only used to build a
derived one, say) has no row in "Channels measured". An object CellProfiler produces but that isn't
assigned to a role has no row in the objects table. A per-object column that isn't numeric (a
varchar or blob measurement) has no row in the measurements table, since it cannot enter a float
matrix.

Being absent from those three tables does not mean being absent from the file. The fourth table,
"Also exported, in uns", lists what the export carries that isn't a per-cell column: which file
each channel was read from, plate and well metadata, object counts, thresholds, per-module timings,
the pipeline text itself. Those are per-image-set or per-run facts, so they go to
``uns["cellprofiler"]["image"]`` and its neighbours, one row per image set, instead of being
repeated on every cell. Read it together with the first three: if a channel is missing from
"Channels measured", ``uns["cellprofiler"]["channels"]`` still names it and ``["image"]`` still
records the file it came from.

The channels table has one row per channel that some exported measurement reads: the module that
loaded it (source "file") or computed it from other images (source "pipeline"). The objects table
has one row per role (primary, secondary, tertiary): the object filling that role, and the module
that produced it, from an image (source "pipeline") or a label file (source "file").

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

**Name in the joined .h5ad** is the exact column (for *X*) or ``obs`` key (for *obs*) the value
will have. A per-object file, if "Also write one file per object?" is enabled, uses this same name
without the ``<Object>__`` prefix.

Type in the filter box to narrow all three tables by name, module, category, destination, or
exported name. The tables reflect the role settings as they are right now; change Automatic vs
Manual or the role objects and press the button again to refresh them. Turn on "Show advanced
features?" below for "Also write the mapping tables to .uns?", which saves this same information
inside the ``.h5ad`` on every run, not just when you press this button.""")
        self.wants_advanced = Binary(
            "Show advanced features?", False,
            doc="""Select "Yes" to also show:

-  **How to build row names**, and the tag list it reveals when set to Manual: which Metadata tags
   name each row.
-  **Auto-configure from this pipeline**: read both the objects to track and the row-name tags off
   the modules above, fill them in, and pin them.
-  **When an object is not exactly one primary + one secondary + one tertiary**: how to handle a
   row that doesn't join 1:1:1.
-  **Also write the mapping tables to .uns?**: save the channel/object/measurement tables "See
   where each measurement will land" shows into ``uns["cellprofiler_mapping"]``.

These matter for troubleshooting an unusual pipeline, for pinning an export you intend to
concatenate, or for automating around the export, not for a first run.""")
        self.autoconfig = DoSomething(
            "", "Auto-configure from this pipeline", self.do_autoconfig,
            doc="""Read this pipeline and fill in both of the things the export otherwise decides for
itself, then explain every choice in one dialog.

-  The **primary/secondary/tertiary objects**, from IdentifyPrimary/Secondary/TertiaryObjects and
   FilterObjects.
-  The **Metadata tags that name each row**, from the tags the pipeline extracts.

Both modes switch to Manual, so a later pipeline edit cannot silently change which objects are
exported or how rows are named. That is the point of the button: it turns whatever the pipeline
currently implies into settings you can see and keep.

Warnings (missing plate metadata, ambiguous objects, and so on) appear in the same dialog. Whether
the chosen tags give every image set a different key can only be checked against real values, so
the run still appends ``img<n>`` if they turn out not to.""")
        self.policy = Choice(
            "When an object is not exactly one primary + one secondary + one tertiary", POLICY_CHOICES,
            doc="Flag: keep the row and mark obs['qc_flag']. Drop: remove it. Error: abort the run.")
        self.wants_mapping_uns = Binary(
            "Also write the mapping tables to .uns?", False,
            doc="""Select "Yes" to save the same channel/object/measurement tables "See where each
measurement will land" shows into ``uns["cellprofiler_mapping"]`` on every run, so the mapping
travels with the ``.h5ad`` rather than staying something you have to reopen CellProfiler to see.

They are written in AnnData's dataframe encoding, so ``adata.uns["cellprofiler_mapping"]["measurements"]``
reads back as a pandas DataFrame you can filter and sort directly, under the keys ``"channels"``,
``"objects"`` and ``"measurements"``.

This is separate from the provenance the export always writes (the pipeline text, per-image-set
columns, module settings). The "Also exported, in uns" table in the preview lists all of it.""")
        self.wants_overwrite = Binary("Overwrite existing files without warning?", True)

    def settings(self):
        return [self.directory, self.prefix, self.wants_per_object, self.policy, self.role_mode,
                self.primary_object, self.secondary_object, self.tertiary_object,
                self.wants_mapping_uns, self.wants_advanced, self.wants_overwrite,
                self.sample_key_mode, self.sample_key_tags]

    def visible_settings(self):
        result = [self.directory, self.prefix, self.wants_per_object, self.role_mode]
        if self.role_mode == ROLE_MANUAL:
            result += [self.primary_object, self.secondary_object, self.tertiary_object]
        result += [self.preview, self.wants_advanced]
        if self.wants_advanced.value:
            result += [self.sample_key_mode]
            if self.sample_key_mode == ROLE_MANUAL:
                result += [self.sample_key_tags]
            result += [self.autoconfig, self.policy, self.wants_mapping_uns]
        return result + [self.wants_overwrite]

    def _sample_tags(self):
        """The Manual tag list, or None to detect. Mirrors _roles()."""
        if self.sample_key_mode != ROLE_MANUAL:
            return None
        return parse_tags(self.sample_key_tags.value)

    # ---- auto-configuration button -------------------------------------------------------
    _gui_pipeline = None

    def on_activated(self, workspace):
        # the GUI calls this when the module is selected; it is the only place a module can grab
        # the pipeline for use in a DoSomething callback (same pattern as the Metadata module)
        self._gui_pipeline = workspace.pipeline

    def on_deactivated(self):
        self._gui_pipeline = None

    def apply_autoconfig(self, pipeline):
        """Detect the objects to track and the row-name tags afresh (ignoring any Manual values),
        pin both as Manual, and return the explanation text. Raises RoleError when the pipeline
        identifies no objects at all.

        Pinning both is the point: after this the export depends on what the pipeline looked like
        when the button was pressed, not on what it looks like at run time, so a later edit cannot
        quietly change which objects are exported or how rows are named."""
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

        detected, note = detect_sample_tags(qualify(t) for t in ctx.metadata_tags)
        self.sample_key_tags.value = ",".join(t[len(TAG_PREFIX):] for t in detected)
        self.sample_key_mode.value = ROLE_MANUAL
        if detected:
            lines.append("row names = <%s>_<label>, %s"
                         % (">_<".join(t[len(TAG_PREFIX):] for t in detected), note))
        else:
            lines.append("row names = img<n>_<label>. " + note[0].upper() + note[1:] +
                         ". Add a Metadata module extracting Plate, Well and Site for readable names.")

        text = ("These settings were filled in from the pipeline:\n\n- " + "\n- ".join(lines) +
                "\n\nBoth modes were switched to Manual so that later pipeline edits cannot "
                "silently change which objects are exported or how rows are named.")
        if detected:
            text += ("\n\nWhether those tags give every image set a different key can only be "
                     "checked against real values, so the run still appends img<n> if they do not.")
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
        _show_feature_preview(channel_report(ctx), object_report(ctx), measurement_report(ctx),
                              uns_report(ctx, self.wants_mapping_uns.value),
                              row_name_preview(ctx, self._sample_tags()))

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
        export = build_export(
            ctx, m, policy=self.policy.value.lower(),
            exporter_settings={s.text: s.value_text for s in self.settings()},
            wants_mapping_uns=self.wants_mapping_uns.value,
            sample_tags=self._sample_tags())
        joined = export.joined
        LOGGER.info("ExportToAnnData: row names are %s_<label>, %s",
                    "_".join(export.naming.parts), export.naming.note)
        for message in export.advice:
            LOGGER.warning("ExportToAnnData: %s", message)
        LOGGER.info("ExportToAnnData: writing %d cells x %d features (~%.1f MB of float32) to %s",
                    joined.X.shape[0], joined.X.shape[1], joined.X.nbytes / 1e6, paths["__joined__"])
        write_h5ad(joined, paths["__joined__"])
        LOGGER.info("ExportToAnnData: wrote %s (%d cells x %d features)", paths["__joined__"], *joined.X.shape)
        for obj, path in paths.items():
            if obj == "__joined__":
                continue
            write_h5ad(export.per_object[obj], path)
            LOGGER.info("ExportToAnnData: wrote %s", path)

    def upgrade_settings(self, setting_values, variable_revision_number, module_name):
        if variable_revision_number == 1:
            # v1 has no "write mapping tables" setting. Default it to No, insert it right before
            # "Overwrite existing files without warning?" (the position settings() puts it at now),
            # and leave every other value untouched, so pipelines saved before this feature existed
            # keep loading.
            setting_values = setting_values[:8] + ["No"] + setting_values[8:]
            variable_revision_number = 2
        if variable_revision_number == 2:
            # v2 has no "Show advanced features?" setting. Same treatment: default No, insert right
            # before "Overwrite existing files without warning?", which the new setting shifts one
            # slot further down.
            setting_values = setting_values[:9] + ["No"] + setting_values[9:]
            variable_revision_number = 3
        if variable_revision_number == 3:
            # v3 has no row-name settings. These two go on the end, where settings() puts them, so
            # nothing already saved shifts position. Automatic is the default, which means a
            # pipeline saved before this existed may get different row names than it did then:
            # names now use whatever Plate/Well/Site-like tags it has instead of falling back to
            # img<n> whenever one was missing. That is the point of the change, and the resolved
            # scheme is logged and recorded in uns; set Manual with Plate,Well,Site for the old
            # all-or-nothing behavior.
            setting_values = setting_values + [ROLE_AUTO, "Plate,Well,Site"]
            variable_revision_number = 4
        return setting_values, variable_revision_number

    def volumetric(self):
        return False


def _show_feature_preview(channel_rows, object_rows, measurement_rows, uns_rows, row_names=""):
    """The "See where each measurement will land" dialog: a filterable set of tables
    (scverse_export.preview.render_report_html) for channel, object, and measurement provenance, plus
    what the export carries in uns, for whatever objects are currently configured.

    wx is imported here, not at module scope, so the module stays importable headless when wx
    isn't installed, the same reason do_autoconfig/prepare_run import it locally instead of at the
    top of the file.
    """
    import wx
    import wx.html

    class _FeaturePreviewDialog(wx.Dialog):
        def __init__(self, channel_rows, object_rows, measurement_rows, uns_rows, row_names):
            super().__init__(None, title="ExportToAnnData feature preview", size=(1100, 700),
                             style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER | wx.CLOSE_BOX)
            self.channel_rows = channel_rows
            self.object_rows = object_rows
            self.measurement_rows = measurement_rows
            self.uns_rows = uns_rows
            self.row_names = row_names
            sizer = wx.BoxSizer(wx.VERTICAL)

            filter_row = wx.BoxSizer(wx.HORIZONTAL)
            filter_row.Add(wx.StaticText(self, label="Filter (channel, object, measurement, module, destination):"),
                           0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
            self.filter_ctrl = wx.TextCtrl(self)
            self.filter_ctrl.Bind(wx.EVT_TEXT, self.on_filter)
            filter_row.Add(self.filter_ctrl, 1, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
            sizer.Add(filter_row, 0, wx.EXPAND)

            self.html_window = wx.html.HtmlWindow(self, style=wx.html.HW_SCROLLBAR_AUTO)
            # A dense HTML table reads small in HtmlWindow at its own default size, so this pushes
            # the base size up twice: SetStandardFonts raises what the widget treats as "normal"
            # text, on top of which render_report_html wraps everything in an explicit <font
            # size="4">, one HTML level above normal. Either alone left the table looking small on
            # some platforms; both together give a size increase that actually shows up.
            system_size = wx.SystemSettings.GetFont(wx.SYS_DEFAULT_GUI_FONT).GetPointSize()
            self.html_window.SetStandardFonts(size=max(system_size + 2, 12))
            self.html_window.SetPage(render_report_html(
                self.channel_rows, self.object_rows, self.measurement_rows, self.uns_rows,
                row_names=self.row_names))
            sizer.Add(self.html_window, 1, wx.EXPAND | wx.ALL, 5)

            close_button = wx.Button(self, id=wx.ID_OK, label="Close")
            close_button.Bind(wx.EVT_BUTTON, lambda evt: self.EndModal(wx.ID_OK))
            close_button.SetDefault()
            sizer.Add(close_button, 0, wx.ALIGN_RIGHT | wx.ALL, 8)

            self.SetSizer(sizer)

        def on_filter(self, evt):
            self.html_window.SetPage(render_report_html(
                self.channel_rows, self.object_rows, self.measurement_rows, self.uns_rows,
                self.filter_ctrl.GetValue(), row_names=self.row_names))

    dlg = _FeaturePreviewDialog(channel_rows, object_rows, measurement_rows, uns_rows, row_names)
    try:
        dlg.ShowModal()
    finally:
        dlg.Destroy()
