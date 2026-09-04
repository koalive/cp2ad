"""The shared table-building step. Needs neither cellprofiler_core nor the recorded fixtures: the
Measurements object CellProfiler passes in is replaced by a stand-in that answers the four methods
assemble.py actually calls.
"""
import numpy
import pytest

from scverse_export.advice import advice as _advice
from scverse_export.assemble import build_object_table, join_tables, provenance
from scverse_export.export import Export, build_export
from scverse_export.introspect import ChannelInfo, Context, ObjectInfo
from scverse_export.names import Feature
from scverse_export.preview import channel_report, mapping_to_uns, measurement_report, object_report


def _feature(obj, cp_name, category, measurement, image=None, coltype="float"):
    return Feature(object=obj, cp_name=cp_name, category=category, measurement=measurement,
                   image=image, coltype=coltype, module_num=1, module_name="Test")


class FakeMeasurements:
    """{object: {feature: {image_number: value}}}, answering the Measurements methods assemble.py
    uses. Enough to exercise the whole join without a CellProfiler install."""

    def __init__(self, data, image_numbers):
        self.data = data
        self.image_numbers = list(image_numbers)

    def get_image_numbers(self):
        return self.image_numbers

    def has_feature(self, obj, feature):
        return feature in self.data.get(obj, {})

    def get_feature_names(self, obj):
        return list(self.data.get(obj, {}).keys())

    def get_measurement(self, obj, feature, image_numbers):
        series = self.data[obj][feature]
        if isinstance(image_numbers, (list, tuple)):
            return [series.get(n) for n in image_numbers]
        return series.get(image_numbers)


@pytest.fixture
def ctx():
    features = [
        _feature("Cells", "AreaShape_Area", "AreaShape", "Area"),
        _feature("Cells", "Location_Center_X", "Location", "Center_X"),
        _feature("Cells", "Location_Center_Y", "Location", "Center_Y"),
        _feature("Cells", "Parent_Nuclei", "Parent", "Nuclei", coltype="integer"),
        _feature("Cells", "Intensity_MeanIntensity_DNA", "Intensity", "MeanIntensity", image="DNA"),
        _feature("Nuclei", "AreaShape_Area", "AreaShape", "Area"),
        _feature("Nuclei", "Location_Center_X", "Location", "Center_X"),
        _feature("Nuclei", "Location_Center_Y", "Location", "Center_Y"),
        _feature("Image", "Metadata_Well", "Metadata", "Well", coltype="varchar(128)"),
    ]
    return Context(
        channels=["DNA"],
        channel_info={"DNA": ChannelInfo("DNA", 1, "NamesAndTypes", "file")},
        objects={"Cells": ObjectInfo("Cells", 2, "IdentifySecondaryObjects", role="secondary"),
                 "Nuclei": ObjectInfo("Nuclei", 1, "IdentifyPrimaryObjects", role="primary")},
        roles={"primary": "Nuclei", "secondary": "Cells"},
        features=features, metadata_tags=["Well"])


@pytest.fixture
def measurements():
    def make():
        a = numpy.array([10.0, 20.0, 30.0])
        return FakeMeasurements({
            "Cells": {"AreaShape_Area": {1: a}, "Location_Center_X": {1: a + 1},
                      "Location_Center_Y": {1: a + 2},
                      "Parent_Nuclei": {1: numpy.array([1.0, 2.0, 3.0])},
                      "Intensity_MeanIntensity_DNA": {1: a / 100}},
            "Nuclei": {"AreaShape_Area": {1: a / 2}, "Location_Center_X": {1: a + 3},
                       "Location_Center_Y": {1: a + 4}},
            "Image": {"Count_Cells": {1: 3}, "Count_Nuclei": {1: 3}, "Metadata_Well": {1: "A02"}},
        }, [1])
    return make


SETTINGS = {"File name prefix": "test", "Overwrite existing files without warning?": "Yes"}


def test_returns_one_row_per_cell_with_the_per_object_tables(ctx, measurements):
    export = build_export(ctx, measurements(), exporter_settings=SETTINGS)
    assert isinstance(export, Export)
    assert export.joined.X.shape[0] == 3
    assert sorted(export.per_object) == ["Cells", "Nuclei"]
    assert all(v.startswith(("Cells__", "Nuclei__")) for v in export.joined.var_names)


def test_provenance_is_attached_to_every_table(ctx, measurements):
    """The per-object files are standalone, so each needs the pipeline provenance too."""
    export = build_export(ctx, measurements(), exporter_settings=SETTINGS)
    assert export.joined.uns["cellprofiler"] is export.provenance
    for table in export.per_object.values():
        assert table.uns["cellprofiler"] is export.provenance
    assert export.provenance["exporter"]["settings"] == SETTINGS


def test_mapping_uns_only_when_asked_and_only_on_the_joined_table(ctx, measurements):
    off = build_export(ctx, measurements(), exporter_settings=SETTINGS)
    assert "cellprofiler_mapping" not in off.joined.uns

    on = build_export(ctx, measurements(), exporter_settings=SETTINGS, wants_mapping_uns=True)
    assert set(on.joined.uns["cellprofiler_mapping"]) == {"channels", "objects", "measurements"}
    for table in on.per_object.values():
        assert "cellprofiler_mapping" not in table.uns


def test_advice_is_returned_not_logged(ctx, measurements):
    """Each module logs under its own logger name, which its tests assert on, so build_export
    hands the messages back instead of logging them itself."""
    export = build_export(ctx, measurements(), exporter_settings=SETTINGS)
    assert export.advice == _advice(ctx)
    assert any("Metadata" in message for message in export.advice), export.advice


@pytest.mark.parametrize("policy", ["flag", "drop"])
@pytest.mark.parametrize("wants_mapping_uns", [False, True])
def test_matches_the_inline_body_it_replaced(ctx, measurements, policy, wants_mapping_uns):
    """Guards the extraction from ExportToAnnData.post_run(). The reference below is a verbatim
    copy of the code that was inline there, so a divergence shows up here rather than as a changed
    export nobody noticed."""
    m = measurements()
    prov = provenance(ctx, m, SETTINGS)
    reference_per_object = {obj: build_object_table(ctx, m, obj) for obj in ctx.roles.values()}
    reference = join_tables(ctx, m, reference_per_object, policy=policy)
    reference.uns["cellprofiler"] = prov
    if wants_mapping_uns:
        reference.uns["cellprofiler_mapping"] = mapping_to_uns(
            channel_report(ctx), object_report(ctx), measurement_report(ctx))
    reference_advice = _advice(ctx)

    export = build_export(ctx, measurements(), policy=policy, exporter_settings=SETTINGS,
                          wants_mapping_uns=wants_mapping_uns)

    numpy.testing.assert_array_equal(export.joined.X, reference.X)
    assert export.joined.var_names == reference.var_names
    assert export.joined.obs_names == reference.obs_names
    assert sorted(export.joined.uns) == sorted(reference.uns)
    assert sorted(export.joined.obs) == sorted(reference.obs)
    for key in reference.obs:
        numpy.testing.assert_array_equal(numpy.asarray(export.joined.obs[key]),
                                         numpy.asarray(reference.obs[key]), err_msg=key)
    assert export.advice == reference_advice
    assert sorted(export.per_object) == sorted(reference_per_object)
    for obj, table in reference_per_object.items():
        numpy.testing.assert_array_equal(export.per_object[obj].X, table.X)


# ---- sample-key naming -----------------------------------------------------------------------

def test_naming_is_resolved_once_and_recorded(ctx, measurements):
    """Every object's table has to name its rows the same way, so the scheme is resolved once here
    rather than per table, and recorded in uns so a reader can tell what a row name means."""
    export = build_export(ctx, measurements(), exporter_settings=SETTINGS)
    assert export.naming.tags == ("Metadata_Well",)
    recorded = export.provenance["sample_naming"]
    assert recorded["tags"] == ["Metadata_Well"]
    assert recorded["parts"] == ["Well"]
    assert recorded["mode"] == "Automatic"
    assert recorded["note"]


def test_automatic_naming_uses_the_tags_the_pipeline_has(ctx, measurements):
    """Metadata_Well is the only usable tag, and with one image set it does identify it, so the
    name carries the well. The old all-or-nothing rule discarded it and named rows img1_1 because
    Plate and Site were missing."""
    export = build_export(ctx, measurements(), exporter_settings=SETTINGS)
    assert list(export.joined.obs_names) == ["A02_1", "A02_2", "A02_3"]
    assert export.naming.with_image_number is False


def test_image_number_is_added_when_the_tags_repeat_across_image_sets(ctx):
    """Two fields of view in one well, with no site tag to tell them apart. Checked against the
    real values rather than assumed, so the scheme adapts to what a run actually contains."""
    a = numpy.array([10.0, 20.0])
    m = FakeMeasurements({
        "Cells": {"AreaShape_Area": {1: a, 2: a}, "Location_Center_X": {1: a, 2: a},
                  "Location_Center_Y": {1: a, 2: a},
                  "Parent_Nuclei": {1: numpy.array([1.0, 2.0]), 2: numpy.array([1.0, 2.0])}},
        "Nuclei": {"AreaShape_Area": {1: a, 2: a}, "Location_Center_X": {1: a, 2: a},
                   "Location_Center_Y": {1: a, 2: a}},
        "Image": {"Count_Cells": {1: 2, 2: 2}, "Count_Nuclei": {1: 2, 2: 2},
                  "Metadata_Well": {1: "A02", 2: "A02"}},
    }, [1, 2])
    export = build_export(ctx, m, exporter_settings=SETTINGS)
    assert export.naming.with_image_number is True
    assert list(export.joined.obs_names) == ["A02_img1_1", "A02_img1_2", "A02_img2_1", "A02_img2_2"]
    assert len(set(export.joined.obs_names)) == 4


def test_manual_naming_overrides_detection(ctx, measurements):
    export = build_export(ctx, measurements(), exporter_settings=SETTINGS,
                          sample_tags=("Metadata_Well",))
    assert export.naming.mode == "Manual"
    assert export.provenance["sample_naming"]["mode"] == "Manual"


def test_manual_naming_with_a_missing_tag_falls_back_and_says_so(ctx, measurements):
    export = build_export(ctx, measurements(), exporter_settings=SETTINGS,
                          sample_tags=("Metadata_Plate", "Metadata_Well"))
    assert export.naming.tags == ()
    assert list(export.joined.obs_names) == ["img1_1", "img1_2", "img1_3"]
    assert "no Plate" in export.naming.note


def test_every_table_names_rows_the_same_way(ctx, measurements):
    export = build_export(ctx, measurements(), exporter_settings=SETTINGS)
    joined_names = list(export.joined.obs_names)
    for table in export.per_object.values():
        assert list(table.obs_names) == joined_names
