"""Folder layout and manifest for ExportForSpatialData. Needs neither cellprofiler_core nor the
recorded fixtures."""
import os

import numpy
import pytest

from scverse_export.assemble import Table
from scverse_export.raster import array_info, write_image, write_labels
from scverse_export.samples import SampleNaming, sample_key
from scverse_export.spatial import (ELEMENT_IMAGE, ELEMENT_LABELS, STATUS_FAILED, STATUS_OK,
                                    UNKNOWN_PLATE, ChannelAxisRow, channel_axis_rows,
                                    element_rows, image_path, labels_path, manifest_to_uns,
                                    plate_of, plates_by_image, region_key_column,
                                    region_key_value, safe_segment, selected_channels,
                                    selected_objects, subset_table, table_path)

NAMING = SampleNaming(tags=("Metadata_Well", "Metadata_Field"), with_image_number=True)


# ---- path segments ---------------------------------------------------------------------------

@pytest.mark.parametrize("given,expected", [
    ("P1", "P1"), ("Plate 1", "Plate_1"),
    ("a/b", "a_b"),                     # a slash would silently create a nested folder
    ("../etc", "etc"), ("  P1  ", "P1"), ("", "unnamed"), ("///", "unnamed"),
    ("A01_s1", "A01_s1"), ("plate#2", "plate_2"),
])
def test_safe_segment(given, expected):
    assert safe_segment(given) == expected


def test_paths_are_relative_to_the_plate_folder():
    """The manifest stays valid when the folder is moved or renamed, which it would not if paths
    were absolute."""
    assert image_path("A02_03_img1") == "images/A02_03_img1.h5"
    assert labels_path("A02_03_img1", "Nuclei") == "labels/A02_03_img1/Nuclei.h5"
    assert table_path("cellprofiler") == "tables/cellprofiler.h5ad"
    for path in (image_path("k"), labels_path("k", "N"), table_path("p")):
        assert not path.startswith("/")


def test_labels_path_sanitises_the_object_name():
    assert labels_path("k", "My Objects") == "labels/k/My_Objects.h5"


# ---- plates ----------------------------------------------------------------------------------

def test_plate_of():
    assert plate_of({"Metadata_Plate": "P1"}) == "P1"
    assert plate_of({"Metadata_Plate": 1.0}) == "1"          # metadata is coerced to float in obs
    assert plate_of({}) == UNKNOWN_PLATE
    assert plate_of({"Metadata_Plate": None}) == UNKNOWN_PLATE
    assert plate_of({"Metadata_Plate": ""}) == UNKNOWN_PLATE


def test_plates_by_image_splits_a_run_covering_two_plates():
    values = {1: {"Metadata_Plate": "P1"}, 2: {"Metadata_Plate": "P2"},
              3: {"Metadata_Plate": "P1"}}
    assert plates_by_image(values) == {"P1": [1, 3], "P2": [2]}


def test_plates_by_image_puts_everything_in_one_folder_without_a_plate_tag():
    """The assumption the module warns about: no plate tag means one plate."""
    values = {1: {"Metadata_Well": "A01"}, 2: {"Metadata_Well": "A02"}}
    assert plates_by_image(values) == {UNKNOWN_PLATE: [1, 2]}


# ---- manifest --------------------------------------------------------------------------------

def test_element_rows_read_shape_and_dtype_off_the_written_files(tmp_path):
    """Reading rather than remembering: it confirms the file is there and readable, so a cycle
    that failed shows up without run() having to report anything in the ordinary case."""
    root = str(tmp_path)
    values = {1: {"Metadata_Well": "A02", "Metadata_Field": "03"}}
    key = "A02_03_img1"
    write_image(str(tmp_path / image_path(key)), numpy.zeros((3, 8, 16), dtype="uint16"))
    write_labels(str(tmp_path / labels_path(key, "Nuclei")), numpy.zeros((8, 16), dtype="int32"))

    rows = element_rows(root, [1], values, NAMING, ["Nuclei"])
    assert [r.element_type for r in rows] == [ELEMENT_IMAGE, ELEMENT_LABELS]
    image, labels = rows
    assert image.sample_key == key and image.element_name == "" and image.region_key_value == ""
    assert image.shape == "3,8,16" and image.dtype == "uint16" and image.status == STATUS_OK
    assert labels.element_name == "Nuclei" and labels.region_key_value == f"{key}__Nuclei"
    assert labels.shape == "8,16" and labels.dtype == "int32" and labels.status == STATUS_OK
    assert all(r.error == "" for r in rows)


def test_a_missing_file_is_reported_as_failed_not_as_a_crash(tmp_path):
    """The importer skips failed rows. Without this it would crash opening a file that a failed
    cycle never wrote."""
    values = {1: {"Metadata_Well": "A02", "Metadata_Field": "03"}}
    rows = element_rows(str(tmp_path), [1], values, NAMING, ["Nuclei"])
    assert all(r.status == STATUS_FAILED for r in rows)
    assert all(r.error for r in rows)
    assert all(r.shape == "" and r.dtype == "" for r in rows)


def test_a_recorded_error_is_carried_into_every_row_for_that_image_set(tmp_path):
    values = {1: {"Metadata_Well": "A02", "Metadata_Field": "03"},
              2: {"Metadata_Well": "B04", "Metadata_Field": "01"}}
    key2 = "B04_01_img2"
    write_image(str(tmp_path / image_path(key2)), numpy.zeros((1, 4, 4), dtype="uint16"))
    write_labels(str(tmp_path / labels_path(key2, "Nuclei")), numpy.zeros((4, 4), dtype="int32"))

    rows = element_rows(str(tmp_path), [1, 2], values, NAMING, ["Nuclei"],
                        errors={1: "OSError: disk full"})
    failed = [r for r in rows if r.image_number == 1]
    ok = [r for r in rows if r.image_number == 2]
    assert all(r.status == STATUS_FAILED and r.error == "OSError: disk full" for r in failed)
    assert all(r.status == STATUS_OK for r in ok)


def test_manifest_round_trips_as_dataframes(tmp_path):
    import anndata
    from scverse_export.h5ad import write_h5ad

    values = {1: {"Metadata_Well": "A02", "Metadata_Field": "03"}}
    key = "A02_03_img1"
    write_image(str(tmp_path / image_path(key)), numpy.zeros((2, 4, 4), dtype="uint16"))
    write_labels(str(tmp_path / labels_path(key, "Nuclei")), numpy.zeros((4, 4), dtype="int32"))
    rows = element_rows(str(tmp_path), [1], values, NAMING, ["Nuclei"])
    manifest = manifest_to_uns(rows, channel_axis_rows(["DNA", "Protein"]))

    table = Table(X=numpy.zeros((1, 1), dtype=numpy.float32), obs_names=["a"], var_names=["x"],
                  obsm={"spatial": numpy.zeros((1, 2))},
                  uns={"cellprofiler_mapping": manifest})
    path = tmp_path / "t.h5ad"
    write_h5ad(table, str(path))

    back = anndata.read_h5ad(str(path)).uns["cellprofiler_mapping"]
    assert set(back) == {"elements", "image_channels"}
    assert list(back["image_channels"]["channel"]) == ["DNA", "Protein"]
    assert list(back["image_channels"]["stack_index"]) == [0, 1]
    assert len(back["elements"]) == 2
    assert set(back["elements"]["element_type"]) == {ELEMENT_IMAGE, ELEMENT_LABELS}


def test_channel_axis_rows_number_from_zero():
    assert channel_axis_rows(["A", "B", "C"]) == [
        ChannelAxisRow("A", 0), ChannelAxisRow("B", 1), ChannelAxisRow("C", 2)]


# ---- region_key ------------------------------------------------------------------------------

def test_region_key_value_is_one_function_for_element_name_and_column():
    """The labels element's name and obs["region_key"] have to be the same string, or SpatialData
    rejects the table, so both come from here."""
    assert region_key_value("A02_03_img1", "Nuclei") == "A02_03_img1__Nuclei"


def test_region_key_column_matches_the_manifest(tmp_path):
    table = Table(
        X=numpy.zeros((3, 1), dtype=numpy.float32), obs_names=["a", "b", "c"], var_names=["x"],
        obs={"ImageNumber": numpy.array([1, 1, 2]),
             "Metadata_Well": numpy.array(["A02", "A02", "B04"], dtype=object),
             "Metadata_Field": numpy.array(["03", "03", "01"], dtype=object)})
    values = {1: {"Metadata_Well": "A02", "Metadata_Field": "03"},
              2: {"Metadata_Well": "B04", "Metadata_Field": "01"}}
    column = region_key_column(table, values, NAMING, "Cells")
    assert list(column) == ["A02_03_img1__Cells", "A02_03_img1__Cells", "B04_01_img2__Cells"]

    rows = element_rows(str(tmp_path), [1, 2], values, NAMING, ["Cells"])
    manifest_keys = {r.region_key_value for r in rows if r.element_type == ELEMENT_LABELS}
    assert set(column) == manifest_keys


# ---- subsetting ------------------------------------------------------------------------------

def _table():
    return Table(
        X=numpy.arange(12, dtype=numpy.float32).reshape(4, 3),
        obs_names=["a", "b", "c", "d"], var_names=["x", "y", "z"],
        obs={"ImageNumber": numpy.array([1, 1, 2, 2]),
             "qc_flag": numpy.array(["ok", "ok", "no_primary", "ok"], dtype=object)},
        var={"cp_name": numpy.array(["x", "y", "z"], dtype=object)},
        obsm={"spatial": numpy.arange(8, dtype=float).reshape(4, 2)},
        uns={"qc_summary": {"ok": 3, "no_primary": 1, "primaries_without_secondary": 2},
             "spatialdata_attrs": {"region": "Cells"}})


def test_subset_table_keeps_rows_and_shares_var():
    table = _table()
    subset = subset_table(table, numpy.array([True, False, True, False]))
    numpy.testing.assert_array_equal(subset.X, table.X[[0, 2]])
    assert subset.obs_names == ["a", "c"]
    assert subset.var_names == table.var_names
    numpy.testing.assert_array_equal(subset.var["cp_name"], table.var["cp_name"])
    numpy.testing.assert_array_equal(subset.obs["ImageNumber"], numpy.array([1, 2]))
    numpy.testing.assert_array_equal(subset.obsm["spatial"], table.obsm["spatial"][[0, 2]])


def test_subset_table_recounts_qc_summary():
    """A run-level count copied into every plate would overstate each of them."""
    subset = subset_table(_table(), numpy.array([True, True, False, False]))
    assert subset.uns["qc_summary"]["ok"] == 2
    assert "no_primary" not in subset.uns["qc_summary"]
    # not splittable per plate, so carried over rather than guessed at
    assert subset.uns["qc_summary"]["primaries_without_secondary"] == 2
    assert _table().uns["qc_summary"]["ok"] == 3, "the original must not be mutated"


def test_subset_table_leaves_other_uns_alone():
    subset = subset_table(_table(), numpy.array([True, False, False, False]))
    assert subset.uns["spatialdata_attrs"] == {"region": "Cells"}


# ---- selection -------------------------------------------------------------------------------

class FakeCtx:
    def __init__(self, channels, sources, objects):
        self.channels = channels
        self.channel_info = {c: type("I", (), {"source": s})() for c, s in zip(channels, sources)}
        self.objects = {o: None for o in objects}


def test_selected_channels_defaults_to_raw_images():
    """Derived images stay out by default: a pipeline can make many, and a folder of
    intermediates is not what anyone asked for."""
    ctx = FakeCtx(["DNA", "Protein", "LogDNA"], ["file", "file", "pipeline"], [])
    assert selected_channels(ctx, []) == ["DNA", "Protein"]


def test_selected_channels_honours_an_explicit_list_including_derived_images():
    ctx = FakeCtx(["DNA", "Protein", "LogDNA"], ["file", "file", "pipeline"], [])
    assert selected_channels(ctx, ["LogDNA", "DNA"]) == ["DNA", "LogDNA"]   # pipeline order


def test_selected_objects_defaults_to_every_object():
    """Unlike the table, which is role-scoped. Any segmentation is a valid labels element."""
    ctx = FakeCtx([], [], ["Nuclei", "Cells", "Spots"])
    assert selected_objects(ctx, []) == ["Nuclei", "Cells", "Spots"]
    assert selected_objects(ctx, ["Spots"]) == ["Spots"]


# ---- end to end ------------------------------------------------------------------------------

def test_two_plates_produce_two_self_consistent_folders(tmp_path):
    """The whole composition post_run performs, without CellProfiler: build the table, add
    region_key, split by plate, write the arrays, build each manifest, write each table. Then the
    checks that matter for the importer: every manifest row resolves to a real file, and every
    region_key in a table names a labels element that plate's manifest lists.
    """
    import anndata
    from scverse_export.export import build_export
    from scverse_export.h5ad import write_h5ad
    from scverse_export.introspect import ChannelInfo, Context, ObjectInfo
    from scverse_export.names import Feature
    from scverse_export.samples import stable_sample_naming
    from test_export import FakeMeasurements

    def feature(obj, cp_name, category, measurement, coltype="float"):
        return Feature(object=obj, cp_name=cp_name, category=category, measurement=measurement,
                       coltype=coltype, module_num=1, module_name="T")

    tags = ["Plate", "Well", "Field"]
    features = [feature("Cells", "AreaShape_Area", "AreaShape", "Area"),
                feature("Cells", "Location_Center_X", "Location", "Center_X"),
                feature("Cells", "Location_Center_Y", "Location", "Center_Y"),
                feature("Cells", "Parent_Nuclei", "Parent", "Nuclei", "integer"),
                feature("Nuclei", "AreaShape_Area", "AreaShape", "Area"),
                feature("Nuclei", "Location_Center_X", "Location", "Center_X"),
                feature("Nuclei", "Location_Center_Y", "Location", "Center_Y")]
    features += [feature("Image", f"Metadata_{t}", "Metadata", t, "varchar(128)") for t in tags]
    ctx = Context(
        channels=["DNA", "LogDNA"],
        channel_info={"DNA": ChannelInfo("DNA", 1, "NamesAndTypes", "file"),
                      "LogDNA": ChannelInfo("LogDNA", 2, "ImageMath", "pipeline")},
        objects={"Nuclei": ObjectInfo("Nuclei", 1, "IdentifyPrimaryObjects", role="primary"),
                 "Cells": ObjectInfo("Cells", 2, "IdentifySecondaryObjects", role="secondary"),
                 "Spots": ObjectInfo("Spots", 3, "IdentifyPrimaryObjects")},
        roles={"primary": "Nuclei", "secondary": "Cells"},
        features=features, metadata_tags=tags)

    # image sets 1 and 2 on plate P1, image set 3 on P2
    plates = {1: ("P1", "A01", "01"), 2: ("P1", "A01", "02"), 3: ("P2", "B02", "01")}
    a = numpy.array([10.0, 20.0])
    per_object = {}
    for obj, scale in (("Cells", 1.0), ("Nuclei", 0.5)):
        per_object[obj] = {"AreaShape_Area": {n: a * scale for n in plates},
                           "Location_Center_X": {n: a + 1 for n in plates},
                           "Location_Center_Y": {n: a + 2 for n in plates}}
    per_object["Cells"]["Parent_Nuclei"] = {n: numpy.array([1.0, 2.0]) for n in plates}
    image_data = {"Count_Cells": {n: 2 for n in plates}, "Count_Nuclei": {n: 2 for n in plates}}
    for i, tag in enumerate(tags):
        image_data[f"Metadata_{tag}"] = {n: plates[n][i] for n in plates}
    m = FakeMeasurements({**per_object, "Image": image_data}, sorted(plates))

    # --- what run() does, one cycle at a time ---
    naming = stable_sample_naming([f"Metadata_{t}" for t in tags])
    channels, objects = selected_channels(ctx, []), selected_objects(ctx, [])
    assert channels == ["DNA"], "derived images stay out by default"
    assert sorted(objects) == ["Cells", "Nuclei", "Spots"], "every object gets a label array"
    root = str(tmp_path / "cellprofiler_export")
    values = {n: {f"Metadata_{t}": plates[n][i] for i, t in enumerate(tags)} for n in plates}
    for n in sorted(plates):
        key = sample_key(values[n], n, naming)
        plate_root = os.path.join(root, plate_of(values[n]))
        write_image(os.path.join(plate_root, image_path(key)),
                    numpy.zeros((len(channels), 4, 4), dtype="uint16"))
        for obj in objects:
            write_labels(os.path.join(plate_root, labels_path(key, obj)),
                         numpy.zeros((4, 4), dtype="int32"))

    # --- what post_run() does ---
    export = build_export(ctx, m, exporter_settings={}, wants_mapping_uns=True, naming=naming)
    joined = export.joined
    joined.obs["region_key"] = region_key_column(joined, values, naming, "Cells")
    written = {}
    for plate, plate_images in plates_by_image(values).items():
        plate_root = os.path.join(root, plate)
        table = subset_table(joined, numpy.isin(numpy.asarray(joined.obs["ImageNumber"]),
                                                plate_images))
        rows = element_rows(plate_root, plate_images, values, naming, objects)
        table.uns.setdefault("cellprofiler_mapping", {}).update(
            manifest_to_uns(rows, channel_axis_rows(channels)))
        path = os.path.join(plate_root, table_path("cellprofiler"))
        os.makedirs(os.path.dirname(path), exist_ok=True)
        write_h5ad(table, path)
        written[plate] = path

    # --- the checks the importer depends on ---
    assert sorted(written) == ["P1", "P2"]
    for plate, path in written.items():
        adata = anndata.read_h5ad(path)
        manifest = adata.uns["cellprofiler_mapping"]["elements"]

        # every row of the manifest points at a file that exists and whose shape it reports
        assert set(manifest["status"]) == {STATUS_OK}, plate
        for _, row in manifest.iterrows():
            on_disk = os.path.join(root, plate, row["path"])
            assert os.path.exists(on_disk), on_disk
            assert row["shape"] == ",".join(str(d) for d in array_info(on_disk)[0])

        # every region_key names a labels element the manifest lists, which is what
        # spatialdata_attrs["region"] has to be derived from
        listed = set(manifest[manifest["element_type"] == ELEMENT_LABELS]["region_key_value"])
        assert set(adata.obs["region_key"]) <= listed, plate
        # only the base object is annotated, so the non-role Spots element has no rows
        assert not any(k.endswith("__Spots") for k in set(adata.obs["region_key"]))
        assert any(v.endswith("__Spots") for v in listed), "Spots is still exported as an element"

        # one image element per field of view on this plate, three labels elements each
        n_fovs = len(plates_by_image(values)[plate])
        assert sum(manifest["element_type"] == ELEMENT_IMAGE) == n_fovs
        assert sum(manifest["element_type"] == ELEMENT_LABELS) == n_fovs * 3
        assert adata.n_obs == n_fovs * 2, "two cells per field of view"

        # the channel axis matches what was written
        axis = adata.uns["cellprofiler_mapping"]["image_channels"]
        assert list(axis["channel"]) == channels
        assert list(axis["stack_index"]) == list(range(len(channels)))

    # the two plates' rows partition the run, and their keys never collide
    p1 = anndata.read_h5ad(written["P1"])
    p2 = anndata.read_h5ad(written["P2"])
    assert p1.n_obs + p2.n_obs == joined.X.shape[0]
    assert not (set(p1.obs["region_key"]) & set(p2.obs["region_key"]))
    assert not (set(p1.obs_names) & set(p2.obs_names))
