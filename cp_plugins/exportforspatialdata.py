"""
ExportForSpatialData
====================

**ExportForSpatialData** exports the pixels as well as the numbers: the raw images and the
segmentation label arrays alongside the per-cell table, so the result can be assembled into a real
SpatialData object rather than a feature matrix on its own. A cell stops being a row of numbers and
becomes a labelled region in a specific image, which you can overlay, crop, and re-measure.

CellProfiler cannot write a ``.zarr`` store, so this module writes plain HDF5 arrays in a folder
tree shaped like the store a SpatialData importer will build from it. One folder per plate, because
one SpatialData object per plate is the target:

|

    <prefix>_export/<plate>/
        images/<sample key>.h5            (C, Y, X), one channel stack per field of view
        labels/<sample key>/<object>.h5   (Y, X), integer labels, 0 for background
        tables/<prefix>.h5ad              one row per cell, with the manifest in uns

|

The importer opens the table first and reads ``uns["cellprofiler_mapping"]["elements"]`` to find
and name everything else, so it never has to walk the folder or parse a file name. It is a separate
piece of work and is not shipped here.

The per-cell table is exactly what **ExportToAnnData** writes, including the rule that keeps
position, orientation and identity out of ``X``, plus one added ``obs["region_key"]`` column naming
the label array each row annotates.

Place this module **last** in the pipeline.

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
# scverse_export package cached in sys.modules. With SCVERSE_EXPORT_DEV=1 we purge it here so a
# reload picks up library edits too (settings changes still need a CellProfiler restart).
if os.environ.get("SCVERSE_EXPORT_DEV"):
    for _name in [m for m in sys.modules if m == "scverse_export" or m.startswith("scverse_export.")]:
        del sys.modules[_name]

import numpy

from cellprofiler_core.constants.module import IO_FOLDER_CHOICE_HELP_TEXT
from cellprofiler_core.module import Module
from cellprofiler_core.preferences import (ABSOLUTE_FOLDER_NAME, DEFAULT_INPUT_FOLDER_NAME,
                                           DEFAULT_INPUT_SUBFOLDER_NAME, DEFAULT_OUTPUT_FOLDER_NAME,
                                           DEFAULT_OUTPUT_SUBFOLDER_NAME, get_headless)
from cellprofiler_core.setting import Binary, ValidationError
from cellprofiler_core.setting.choice import Choice
from cellprofiler_core.setting.do_something import DoSomething
from cellprofiler_core.setting.subscriber import (ImageListSubscriber, LabelListSubscriber,
                                                  LabelSubscriber)
from cellprofiler_core.setting.text import Directory, Text

from scverse_export.advice import advice as _advice
from scverse_export.assemble import POLICIES, metadata_by_image, metadata_features
from scverse_export.export import build_export
from scverse_export.h5ad import write_h5ad
from scverse_export.introspect import RoleError, build_context
from scverse_export.raster import write_image, write_labels
from scverse_export.samples import PREFIX as TAG_PREFIX
from scverse_export.samples import detect_sample_tags, parse_tags, sample_key, stable_sample_naming
from scverse_export.spatial import (ERROR_MEASUREMENT, UNKNOWN_PLATE, channel_axis_rows,
                                    element_rows, image_path, labels_path, manifest_to_uns,
                                    plate_of, plates_by_image, region_key_column, selected_channels,
                                    selected_objects, subset_table, table_path)

LOGGER = logging.getLogger(__name__)
POLICY_CHOICES = tuple(p.capitalize() for p in POLICIES)
ROLE_AUTO, ROLE_MANUAL = "Automatic", "Manual"

NO_PLATE_ADVICE = (
    "This pipeline defines no Metadata_Plate, so every image set is assumed to come from one "
    "plate and everything is written to a single folder named '%s'. That is wrong for a run "
    "covering several plates, and nothing in the metadata can tell. Add a Metadata module "
    "extracting Plate to split them, which also makes element names unique across plates so two "
    "exports concatenate without renaming." % UNKNOWN_PLATE)


class ExportForSpatialData(Module):
    module_name = "ExportForSpatialData"
    category = ["File Processing", "Data Tools"]
    variable_revision_number = 1

    # ---- settings ------------------------------------------------------------------------
    def create_settings(self):
        self.directory = Directory(
            "Output file location",
            dir_choices=[DEFAULT_OUTPUT_FOLDER_NAME, DEFAULT_OUTPUT_SUBFOLDER_NAME,
                         ABSOLUTE_FOLDER_NAME, DEFAULT_INPUT_FOLDER_NAME,
                         DEFAULT_INPUT_SUBFOLDER_NAME],
            doc=IO_FOLDER_CHOICE_HELP_TEXT)
        self.directory.dir_choice = DEFAULT_OUTPUT_FOLDER_NAME
        self.prefix = Text(
            "Folder name prefix", "cellprofiler",
            doc="""One folder ``<prefix>_export`` is written here, holding one subfolder per plate.
Each plate folder has ``images/``, ``labels/`` and ``tables/``, and is what a SpatialData importer
turns into one ``.zarr`` store.""")
        self.channels = ImageListSubscriber(
            "Images to export", [],
            doc="""Leave this empty to export every image the pipeline loaded from a file, which is
usually what you want.

Derived images (illumination-corrected, filtered, masked, anything a module computed) are excluded
by default, because a pipeline can produce many of them and a folder full of intermediates is not
worth the disk. Select images here to export exactly those instead, which is how you keep a
derived image, an illumination-corrected channel for instance, in place of or alongside the raw
one.""")
        self.objects = LabelListSubscriber(
            "Segmentations to export", [],
            doc="""Leave this empty to export a label array for every object the pipeline made.

This is deliberately more permissive than the per-cell table, which only covers the primary,
secondary and tertiary objects because joining compartments into one row per cell is what it is
for. A label array has no such constraint, so a pipeline that segments spots or organelles outside
that chain still gets them as viewable elements. Such an element has no rows in the table, which is
a valid outcome: the segmentation is there to look at, it just has no feature vector.

Select objects here to export exactly those instead. Pixels cannot be recovered after the run
without re-running the pipeline, so the default errs on the side of keeping them.""")
        self.role_mode = Choice(
            "How to pick primary / secondary / tertiary objects", (ROLE_AUTO, ROLE_MANUAL),
            doc="""Which objects the per-cell table joins into one row. Automatic follows
IdentifyPrimary/Secondary/TertiaryObjects and FilterObjects. Manual uses exactly the roles you set.

This is separate from "Segmentations to export", which decides what gets a label array on disk.
An object can be exported as a label array without being part of the joined table.""")
        self.primary_object = LabelSubscriber("Primary objects (e.g. nuclei)", "None")
        self.secondary_object = LabelSubscriber("Secondary objects (e.g. cells)", "None")
        self.tertiary_object = LabelSubscriber("Tertiary objects (e.g. cytoplasm)", "None")
        self.wants_advanced = Binary(
            "Show advanced features?", False,
            doc="""Select "Yes" to also show:

-  **How to name fields of view**, and the tag list it reveals when set to Manual: which Metadata
   tags name each image, label array and coordinate system.
-  **Auto-configure from this pipeline**: read both the objects to track and the naming tags off
   the modules above, fill them in, and pin them.
-  **When an object is not exactly one primary + one secondary + one tertiary**: how to handle a
   row that doesn't join 1:1:1.

These matter for pinning an export you intend to concatenate with another, and for troubleshooting
an unusual pipeline, not for a first run.""")
        self.sample_key_mode = Choice(
            "How to name fields of view", (ROLE_AUTO, ROLE_MANUAL),
            doc="""Every image, label array and coordinate system is named after the field of view it
belongs to, and table rows are named ``<sample key>_<label>``, so this one choice decides all of
them at once.

**Automatic** reads the key off the pipeline's Metadata tags, taking at most one plate part (Plate,
Barcode, PlateID), one well part (Well, or Row plus Column together), and one site part (Site,
Field, FieldIndex, Position). Frame, Series and Channel are never used, because they index a z
plane, a timepoint or a channel inside one field of view rather than the field itself.

The image number is always part of the key here, unlike in ExportToAnnData. Files are named one
cycle at a time, possibly in a worker process, while table rows are named at the end in the main
process, and the two have to agree exactly or a manifest row would point at a file written under
another name. Including the image number makes the key unique by construction rather than by a
check that needs every image set at once.

**Manual** uses exactly the tags you list, which pins the scheme so it does not change if you run
the pipeline over a different subset. Worth doing before exporting anything you intend to
concatenate.""")
        self.sample_key_tags = Text(
            "Metadata tags for field-of-view names", "Plate,Well,Site",
            doc="""*(Used only when naming fields of view Manually)*

Comma-separated Metadata tag names, in the order they should appear in the key. The ``Metadata_``
prefix is optional, so ``Well,Field`` and ``Metadata_Well,Metadata_Field`` mean the same thing.
Leave it empty to name fields of view by image number alone.

A tag this pipeline does not have is reported and the key falls back to the image number, rather
than leaving a blank in the middle of a name.""")
        self.autoconfig = DoSomething(
            "", "Auto-configure from this pipeline", self.do_autoconfig,
            doc="""Read this pipeline and fill in both of the things the export otherwise decides for
itself, then explain every choice in one dialog.

-  The **primary/secondary/tertiary objects**, from IdentifyPrimary/Secondary/TertiaryObjects and
   FilterObjects.
-  The **Metadata tags that name each field of view**, from the tags the pipeline extracts.

Both modes switch to Manual, so a later pipeline edit cannot silently change which objects are
exported or how elements are named.""")
        self.policy = Choice(
            "When an object is not exactly one primary + one secondary + one tertiary",
            POLICY_CHOICES,
            doc="Flag: keep the row and mark obs['qc_flag']. Drop: remove it. Error: abort the run.")
        self.wants_overwrite = Binary(
            "Overwrite existing files without warning?", True,
            doc="""Images and label arrays are always overwritten as they are written, one cycle at a
time. This governs the tables, which are checked before the run starts.""")

    def settings(self):
        return [self.directory, self.prefix, self.channels, self.objects, self.role_mode,
                self.primary_object, self.secondary_object, self.tertiary_object,
                self.sample_key_mode, self.sample_key_tags, self.policy, self.wants_advanced,
                self.wants_overwrite]

    def visible_settings(self):
        result = [self.directory, self.prefix, self.channels, self.objects, self.role_mode]
        if self.role_mode == ROLE_MANUAL:
            result += [self.primary_object, self.secondary_object, self.tertiary_object]
        result += [self.wants_advanced]
        if self.wants_advanced.value:
            result += [self.sample_key_mode]
            if self.sample_key_mode == ROLE_MANUAL:
                result += [self.sample_key_tags]
            result += [self.autoconfig, self.policy]
        return result + [self.wants_overwrite]

    # ---- resolved settings ---------------------------------------------------------------
    def _roles(self):
        if self.role_mode != ROLE_MANUAL:
            return None
        roles = {}
        for role, setting in (("primary", self.primary_object), ("secondary", self.secondary_object),
                              ("tertiary", self.tertiary_object)):
            if setting.value not in (None, "", "None"):
                roles[role] = setting.value
        return roles or None

    def _sample_tags(self):
        if self.sample_key_mode != ROLE_MANUAL:
            return None
        return parse_tags(self.sample_key_tags.value)

    def _naming(self, ctx, m):
        """The one naming scheme this run uses for files, elements and rows. Deterministic given
        the pipeline, so run() and post_run() agree without passing state between processes."""
        return stable_sample_naming(metadata_features(ctx, m), self._sample_tags())

    def _root(self):
        return os.path.join(self.directory.get_absolute_path(), f"{self.prefix.value}_export")

    # ---- GUI -----------------------------------------------------------------------------
    _gui_pipeline = None

    def on_activated(self, workspace):
        self._gui_pipeline = workspace.pipeline

    def on_deactivated(self):
        self._gui_pipeline = None

    def apply_autoconfig(self, pipeline):
        """Detect the objects to track and the naming tags afresh, pin both as Manual, and return
        the explanation text."""
        ctx = build_context(pipeline, roles=None)
        lines = []
        for role, setting in (("primary", self.primary_object), ("secondary", self.secondary_object),
                              ("tertiary", self.tertiary_object)):
            name = ctx.roles.get(role)
            setting.value = name if name else "None"
            if name:
                info = ctx.objects[name]
                lines.append('%s = "%s", made by %s (module #%d)'
                             % (role, name, info.module_name, info.module_num))
            else:
                lines.append("%s = None. No %s object exists in this pipeline, so nothing is "
                             "exported for that role." % (role, role))
        self.role_mode.value = ROLE_MANUAL

        detected, note = detect_sample_tags(TAG_PREFIX + t for t in ctx.metadata_tags)
        self.sample_key_tags.value = ",".join(t[len(TAG_PREFIX):] for t in detected)
        self.sample_key_mode.value = ROLE_MANUAL
        if detected:
            lines.append("field-of-view names = <%s>_img<n>, %s"
                         % (">_<".join(t[len(TAG_PREFIX):] for t in detected), note))
        else:
            lines.append("field-of-view names = img<n>. " + note[0].upper() + note[1:])

        text = ("These settings were filled in from the pipeline:\n\n- " + "\n- ".join(lines) +
                "\n\nBoth modes were switched to Manual so that later pipeline edits cannot "
                "silently change which objects are exported or how elements are named.")
        warnings = _advice(ctx) + self._plate_advice(ctx)
        if warnings:
            text += "\n\nWarnings:\n\n- " + "\n- ".join(warnings)
        return text, bool(warnings)

    def do_autoconfig(self):
        import wx
        if self._gui_pipeline is None:
            wx.MessageBox("The pipeline is not available yet - click another module and then this "
                          "one again, then retry.", caption=self.module_name,
                          style=wx.OK | wx.ICON_WARNING)
            return
        try:
            text, warned = self.apply_autoconfig(self._gui_pipeline)
        except RoleError as e:
            wx.MessageBox("Nothing was changed: %s\n\nAdd an IdentifyPrimaryObjects (or similar) "
                          "module before this one, or set the roles by hand in Manual mode." % e,
                          caption=self.module_name + " auto-configuration",
                          style=wx.OK | wx.ICON_ERROR)
            return
        wx.MessageBox(text, caption=self.module_name + " auto-configuration",
                      style=wx.OK | (wx.ICON_WARNING if warned else wx.ICON_INFORMATION))

    # ---- validation ----------------------------------------------------------------------
    def _plate_advice(self, ctx):
        return [] if "Plate" in ctx.metadata_tags else [NO_PLATE_ADVICE]

    def validate_module(self, pipeline):
        if pipeline.modules()[-1] is not self:
            raise ValidationError("ExportForSpatialData must be the last module; measurements made "
                                  "after it are not exported.", self.prefix)
        try:
            build_context(pipeline, roles=self._roles())
        except RoleError as e:
            raise ValidationError(str(e), self.role_mode)

    def validate_module_warnings(self, pipeline):
        try:
            ctx = build_context(pipeline, roles=self._roles())
        except RoleError:
            return  # validate_module already reports this as an error
        messages = _advice(ctx) + self._plate_advice(ctx)
        if messages:
            raise ValidationError("\n\n".join(messages), self.prefix)

    # ---- run -----------------------------------------------------------------------------
    def prepare_run(self, workspace):
        try:
            ctx = build_context(workspace.pipeline, roles=self._roles())
        except RoleError as e:
            LOGGER.error("ExportForSpatialData: %s", e)
            return False
        for message in self._plate_advice(ctx):
            LOGGER.warning("ExportForSpatialData: %s", message)
        root = self._root()
        existing = [os.path.join(root, plate, table_path(self.prefix.value))
                    for plate in os.listdir(root)] if os.path.isdir(root) else []
        existing = [p for p in existing if os.path.exists(p)]
        if existing and not self.wants_overwrite.value:
            if get_headless():
                raise ValueError("ExportForSpatialData: output exists and overwrite is off: "
                                 + ", ".join(existing))
            import wx
            if wx.MessageBox("Overwrite %s?" % ", ".join(existing), caption=self.module_name,
                             style=wx.YES_NO | wx.ICON_QUESTION) != wx.YES:
                return False
        return True

    def run(self, workspace):
        """Write this image set's channel stack and label arrays, then let CellProfiler free them.

        This is where ExportForSpatialData differs from ExportToAnnData structurally. Measurements
        accumulate for the whole run so they can all be read at the end, but pixels do not survive
        that long: CellProfiler processes one image set at a time and does not keep every cycle's
        images alive, which is the right call, since otherwise memory would scale with plate size
        rather than with one field of view.

        A failure here is recorded and swallowed rather than raised. The run continues, the other
        fields of view still export, and post_run puts the reason in the manifest so the importer
        skips that field instead of crashing on a file that was never written.
        """
        m = workspace.measurements
        image_number = int(m.image_set_number)
        try:
            ctx = build_context(workspace.pipeline, roles=self._roles())
            naming = self._naming(ctx, m)
            md = {tag: _image_value(m, tag, image_number)
                  for tag in metadata_features(ctx, m)}
            key = sample_key(md, image_number, naming)
            root = os.path.join(self._root(), plate_of(md))

            channels = selected_channels(ctx, self.channels.value)
            stack = numpy.stack([
                numpy.asarray(workspace.image_set.get_image(name).pixel_data)
                for name in channels]) if channels else numpy.zeros((0, 0, 0))
            if stack.ndim != 3:
                raise ValueError(
                    "expected one 2D plane per channel, got a stack of shape %s. This module does "
                    "not support 3D or colour images yet." % (stack.shape,))
            write_image(os.path.join(root, image_path(key)), stack)

            for obj in selected_objects(ctx, self.objects.value):
                labels = workspace.object_set.get_objects(obj).segmented
                write_labels(os.path.join(root, labels_path(key, obj)), numpy.asarray(labels))
        except Exception as exc:
            LOGGER.error("ExportForSpatialData: image set %d was not exported: %s",
                         image_number, exc, exc_info=True)
            m.add_image_measurement(ERROR_MEASUREMENT, "%s: %s" % (type(exc).__name__, exc))

    def display(self, workspace, figure):
        figure.set_subplots((1, 1))
        figure.subplot_table(0, 0, [["Output", "one folder per plate, written as the run goes"]])

    def post_run(self, workspace):
        """Write one table per plate, each carrying the manifest for its own folder."""
        if workspace.pipeline.test_mode:
            return
        m = workspace.measurements
        try:
            ctx = build_context(workspace.pipeline, roles=self._roles())
        except RoleError as e:
            LOGGER.error("ExportForSpatialData: %s", e)
            return

        naming = self._naming(ctx, m)
        export = build_export(ctx, m, policy=self.policy.value.lower(),
                              exporter_settings={s.text: s.value_text for s in self.settings()},
                              wants_mapping_uns=True, naming=naming)
        for message in export.advice + self._plate_advice(ctx):
            LOGGER.warning("ExportForSpatialData: %s", message)
        LOGGER.info("ExportForSpatialData: naming fields of view %s, %s",
                    "_".join(naming.parts), naming.note)

        md_feats = metadata_features(ctx, m)
        image_numbers = [int(n) for n in m.get_image_numbers()]
        values = metadata_by_image(m, md_feats, image_numbers)

        base_object = ctx.roles.get("secondary") or ctx.roles.get("primary")
        joined = export.joined
        joined.obs["region_key"] = region_key_column(joined, values, naming, base_object)

        errors = {n: _image_value(m, ERROR_MEASUREMENT, n) for n in image_numbers}
        errors = {n: str(v) for n, v in errors.items() if v}
        channels = selected_channels(ctx, self.channels.value)
        objects = selected_objects(ctx, self.objects.value)

        root = self._root()
        for plate, plate_images in plates_by_image(values).items():
            plate_root = os.path.join(root, plate)
            keep = numpy.isin(numpy.asarray(joined.obs["ImageNumber"]), plate_images)
            table = subset_table(joined, keep)
            rows = element_rows(plate_root, plate_images, values, naming, objects, errors)
            table.uns.setdefault("cellprofiler_mapping", {}).update(
                manifest_to_uns(rows, channel_axis_rows(channels)))
            path = os.path.join(plate_root, table_path(self.prefix.value))
            os.makedirs(os.path.dirname(path), exist_ok=True)
            write_h5ad(table, path)
            failed = sum(1 for r in rows if r.status == "failed")
            LOGGER.info("ExportForSpatialData: wrote %s (%d cells x %d features, %d elements%s)",
                        path, table.X.shape[0], table.X.shape[1], len(rows),
                        ", %d failed" % failed if failed else "")

    def upgrade_settings(self, setting_values, variable_revision_number, module_name):
        return setting_values, variable_revision_number

    def volumetric(self):
        return False


def _image_value(m, feature, image_number):
    """One Image measurement, or None when it was never recorded for this image set. Errors are
    only recorded for cycles that failed, so absence is the ordinary case."""
    try:
        if not m.has_feature("Image", feature):
            return None
        return m.get_measurement("Image", feature, image_number)
    except Exception:
        return None
