import dataclasses
import logging

import numpy
import pytest
import cellprofiler_core.measurement

from scverse_export.assemble import Table, build_object_table
from scverse_export.introspect import build_context
from scverse_export.names import Feature, to_cpm_names


def make_measurements(meas_arrays, image_metadata=None, image_numbers=(1,), without_counts=()):
    """Measurements with the probe's real per-object arrays replicated into each image set."""
    m = cellprofiler_core.measurement.Measurements()
    for n in image_numbers:
        m.next_image_set(n)
        for obj, feats in meas_arrays.items():
            for feat, arr in feats.items():
                m.add_measurement(obj, feat, arr)
            if obj not in without_counts:
                m.add_image_measurement(f"Count_{obj}", len(next(iter(feats.values()))))
        for k, v in (image_metadata or {}).items():
            m.add_image_measurement(k, v if not isinstance(v, (list, tuple)) else v[n - 1])
    return m


@pytest.fixture
def ctx(fake_pipeline):
    return build_context(fake_pipeline)


def test_object_table_shape_and_names(ctx, meas_arrays, cpm_columns):
    m = make_measurements(meas_arrays)
    t = build_object_table(ctx, m, "Cells")
    assert isinstance(t, Table)
    assert t.X.dtype == numpy.float32 and t.X.shape[0] == 289
    assert len(t.var_names) == t.X.shape[1] == len(set(t.var_names))
    cols = set(cpm_columns)
    for f in ctx.features:
        if f.object != "Cells":
            continue
        for name, backend in to_cpm_names(f, ctx.channels):
            if backend == "cp_measure" and name in t.var_names:
                assert name in cols, name
    assert "Area" in t.var_names and "Intensity_MeanIntensity__DNA" in t.var_names
    assert "Correlation_Pearson__DNA__PH3" in t.var_names and "Correlation_Pearson__PH3__DNA" in t.var_names


def test_values_are_the_cp_values(ctx, meas_arrays):
    m = make_measurements(meas_arrays)
    t = build_object_table(ctx, m, "Cells")
    j = t.var_names.index("Intensity_MeanIntensity__DNA")
    numpy.testing.assert_allclose(t.X[:, j], meas_arrays["Cells"]["Intensity_MeanIntensity_DNA"].astype("float32"))
    assert t.var["cp_name"][j] == "Intensity_MeanIntensity_DNA"
    assert t.var["channel"][j] == "DNA" and t.var["module_name"][j] == "MeasureObjectIntensity"


def test_obs_without_plate_metadata(ctx, meas_arrays):
    t = build_object_table(ctx, make_measurements(meas_arrays), "Cells")
    assert t.obs_names[:2] == ["img1_1", "img1_2"]
    assert list(t.obs["region"][:1]) == ["Cells"] and t.obs["label_id"][2] == 3 and t.obs["ImageNumber"][0] == 1
    assert t.uns["spatialdata_attrs"] == {"region": "Cells", "region_key": "region", "instance_key": "label_id"}
    assert t.obsm["spatial"].shape == (289, 2)
    numpy.testing.assert_allclose(t.obsm["spatial"][:, 0], meas_arrays["Cells"]["Location_Center_X"])


def test_obs_with_plate_metadata_and_two_sites(ctx, meas_arrays):
    md = {"Metadata_Plate": "P1", "Metadata_Well": "A01", "Metadata_Site": [1, 2]}
    m = make_measurements(meas_arrays, md, image_numbers=(1, 2))
    t = build_object_table(ctx, m, "Nuclei")
    assert t.X.shape[0] == 2 * 289
    assert t.obs_names[0] == "P1_A01_1_1" and t.obs_names[289] == "P1_A01_2_1"
    assert list(t.obs["Metadata_Site"][[0, 289]]) == [1, 2]
    assert t.obs["Metadata_Plate"][0] == "P1"


def test_location_and_orientation_go_to_obs_not_x(ctx, meas_arrays):
    """Position/orientation measurements must never enter X. They'd bias morphological similarity
    on where in the image, or how, an object happened to be imaged, not its biology."""
    t = build_object_table(ctx, make_measurements(meas_arrays), "Cells")
    for name in ("Center_X", "Center_Y", "Orientation"):
        assert name not in t.var_names, name
        assert name in t.obs, name
    numpy.testing.assert_allclose(t.obs["Center_X"], meas_arrays["Cells"]["Location_Center_X"].astype("float32"))
    numpy.testing.assert_allclose(t.obs["Orientation"],
                                  meas_arrays["Cells"]["AreaShape_Orientation"].astype("float32"))
    # obsm["spatial"] (squidpy/spatialdata-facing) is unaffected by the X/obs split
    numpy.testing.assert_allclose(t.obsm["spatial"][:, 0], meas_arrays["Cells"]["Location_Center_X"])


def test_identity_and_linkage_go_to_obs_not_x(ctx, meas_arrays):
    """An object's own arbitrary label, a Parent_* reference to another object's label, and a
    Children_*_Count summary are bookkeeping, not morphology, for the same reason as position and
    orientation, so they must never enter X either."""
    t = build_object_table(ctx, make_measurements(meas_arrays), "Cells")
    for name in ("Number_Object_Number", "Parent_Nuclei", "Children_Cytoplasm_Count"):
        assert name not in t.var_names, name
        assert name in t.obs, name
    numpy.testing.assert_allclose(t.obs["Number_Object_Number"],
                                  meas_arrays["Cells"]["Number_Object_Number"].astype("float32"))
    numpy.testing.assert_allclose(t.obs["Parent_Nuclei"],
                                  meas_arrays["Cells"]["Parent_Nuclei"].astype("float32"))
    numpy.testing.assert_allclose(t.obs["Children_Cytoplasm_Count"],
                                  meas_arrays["Cells"]["Children_Cytoplasm_Count"].astype("float32"))


def test_missing_feature_is_nan_and_counted(ctx, meas_arrays):
    arrays = {o: dict(f) for o, f in meas_arrays.items()}
    del arrays["Cells"]["AreaShape_Area"]
    t = build_object_table(ctx, make_measurements(arrays), "Cells")
    j = t.var_names.index("Area")
    assert numpy.isnan(t.X[:, j]).all()
    assert (t.obs["n_missing_features"] >= 1).all()


from scverse_export.assemble import JoinError, PolicyError, join_tables, provenance


def _tables(ctx, m):
    return {obj: build_object_table(ctx, m, obj) for obj in ctx.roles.values()}


def test_join_one_row_per_cell(ctx, meas_arrays):
    m = make_measurements(meas_arrays)
    t = join_tables(ctx, m, _tables(ctx, m))
    assert t.X.shape[0] == 289
    assert t.var_names[0].startswith(("Cells__", "Nuclei__", "Cytoplasm__"))
    assert "Cells__Intensity_MeanIntensity__DNA" in t.var_names and "Nuclei__Area" in t.var_names
    assert set(t.var["region"]) == {"Nuclei", "Cells", "Cytoplasm"}
    assert list(t.obs["qc_flag"][:3]) == ["ok", "ok", "ok"]
    assert t.uns["spatialdata_attrs"]["region"] == "Cells"
    assert t.uns["qc_summary"]["ok"] == 289
    jn = t.var_names.index("Nuclei__Area")
    numpy.testing.assert_allclose(t.X[:, jn], meas_arrays["Nuclei"]["AreaShape_Area"].astype("float32"))


def test_join_location_and_orientation_go_to_prefixed_obs(ctx, meas_arrays):
    """Joined table: each compartment's Location/orientation columns land in obs, prefixed the same
    way var_names are (`{object}__{name}`). They never appear in var/X."""
    m = make_measurements(meas_arrays)
    t = join_tables(ctx, m, _tables(ctx, m))
    for bad in ("Cells__Center_X", "Nuclei__Center_X", "Cells__Orientation", "Center_X", "Orientation"):
        assert bad not in t.var_names, bad
    for good in ("Cells__Center_X", "Cells__Center_Y", "Cells__Orientation",
                 "Nuclei__Center_X", "Nuclei__Center_Y", "Nuclei__Orientation"):
        assert good in t.obs, good
    numpy.testing.assert_allclose(t.obs["Nuclei__Center_X"],
                                  meas_arrays["Nuclei"]["Location_Center_X"].astype("float32"))
    numpy.testing.assert_allclose(t.obs["Nuclei__Orientation"],
                                  meas_arrays["Nuclei"]["AreaShape_Orientation"].astype("float32"))


def test_join_identity_and_linkage_go_to_prefixed_obs(ctx, meas_arrays):
    """Joined table: Number_Object_Number/Parent_*/Children_*_Count land in obs, prefixed like every
    other compartment column. They never appear in var/X: an arbitrary label or a reference to
    another row is not a measurement of this cell."""
    m = make_measurements(meas_arrays)
    t = join_tables(ctx, m, _tables(ctx, m))
    for bad in ("Cells__Number_Object_Number", "Cells__Parent_Nuclei", "Cells__Children_Cytoplasm_Count",
                "Number_Object_Number", "Parent_Nuclei"):
        assert bad not in t.var_names, bad
    for good in ("Cells__Number_Object_Number", "Cells__Parent_Nuclei", "Cells__Children_Cytoplasm_Count"):
        assert good in t.obs, good
    numpy.testing.assert_allclose(t.obs["Cells__Number_Object_Number"],
                                  meas_arrays["Cells"]["Number_Object_Number"].astype("float32"))
    numpy.testing.assert_allclose(t.obs["Cells__Parent_Nuclei"],
                                  meas_arrays["Cells"]["Parent_Nuclei"].astype("float32"))


def test_join_uses_parent_columns_not_row_order(ctx, meas_arrays):
    arrays = {o: dict(f) for o, f in meas_arrays.items()}
    # reverse the parent mapping: secondary i now claims primary (290 - i)
    n = len(arrays["Cells"]["Parent_Nuclei"])
    arrays["Cells"]["Parent_Nuclei"] = numpy.arange(n, 0, -1)
    m = make_measurements(arrays)
    t = join_tables(ctx, m, _tables(ctx, m))
    jn = t.var_names.index("Nuclei__Area")
    numpy.testing.assert_allclose(t.X[:, jn], arrays["Nuclei"]["AreaShape_Area"][::-1].astype("float32"))


def test_policy_flag_drop_error(ctx, meas_arrays):
    arrays = {o: dict(f) for o, f in meas_arrays.items()}
    arrays["Cells"]["Children_Cytoplasm_Count"] = arrays["Cells"]["Children_Cytoplasm_Count"].copy()
    arrays["Cells"]["Children_Cytoplasm_Count"][:5] = 0
    arrays["Cells"]["Parent_Nuclei"] = arrays["Cells"]["Parent_Nuclei"].copy()
    arrays["Cells"]["Parent_Nuclei"][5] = 0
    m = make_measurements(arrays)
    flagged = join_tables(ctx, m, _tables(ctx, m), policy="flag")
    assert list(flagged.obs["qc_flag"][:7]) == ["no_tertiary"] * 5 + ["no_primary", "ok"]
    assert flagged.uns["qc_summary"] == {"ok": 283, "no_tertiary": 5, "no_primary": 1,
                                         "primaries_without_secondary": 1}
    jc = flagged.var_names.index("Cytoplasm__Area")
    assert numpy.isnan(flagged.X[:5, jc]).all()
    dropped = join_tables(ctx, m, _tables(ctx, m), policy="drop")
    assert dropped.X.shape[0] == 283 and (dropped.obs["qc_flag"] == "ok").all()
    with pytest.raises(PolicyError):
        join_tables(ctx, m, _tables(ctx, m), policy="error")


def test_child_counts_for_non_role_objects(ctx, meas_arrays):
    m = make_measurements(meas_arrays)
    t = join_tables(ctx, m, _tables(ctx, m))
    assert "count_PH3" in t.obs   # Nuclei has Children_PH3_Count via RelateObjects; surfaced on the cell row
    assert t.obs["count_PH3"].sum() == meas_arrays["Nuclei"]["Children_PH3_Count"].sum()
    # the same data also appears under the generic per-compartment extrinsic-obs naming, since
    # Children_*_Count is excluded from X/var for every compartment, not only non-role children
    assert "Nuclei__Children_PH3_Count" in t.obs
    numpy.testing.assert_allclose(t.obs["Nuclei__Children_PH3_Count"],
                                  meas_arrays["Nuclei"]["Children_PH3_Count"].astype("float32"))
    assert "Children_PH3_Count" not in t.var_names and "Nuclei__Children_PH3_Count" not in t.var_names


def test_provenance(ctx, meas_arrays):
    m = make_measurements(meas_arrays, {"Metadata_Frame": 0, "Threshold_FinalThreshold_Nuclei": 0.1})
    p = provenance(ctx, m, {"prefix": "x"})
    assert p["roles"] == ctx.roles and p["channels"] == ["DNA", "PH3", "cellbody"]
    assert p["objects"]["Cells"]["module_name"] == "IdentifySecondaryObjects"
    assert p["image"]["ImageNumber"] == [1] and p["image"]["Threshold_FinalThreshold_Nuclei"] == [0.1]
    assert p["image"]["Metadata_Frame"] == [0]
    assert [mod["name"] for mod in p["modules"]][:2] == ["Images", "Metadata"]
    assert p["exporter"]["settings"] == {"prefix": "x"}
    # no Experiment measurements were recorded by make_measurements -> version/run_timestamp/
    # pipeline_text stay None (the try/except around get_feature_names("Experiment") swallows the
    # missing-group error and leaves `experiment` empty).
    assert p["version"] is None and p["run_timestamp"] is None and p["pipeline_text"] is None
    # provenance drops the collapsed "settings" dict from each module entry and keeps the
    # complete, ordered "setting_values" pairs instead.
    assert all("settings" not in mod and "setting_values" in mod for mod in p["modules"])


def test_join_primary_only_role(fake_pipeline, meas_arrays):
    # Manual mode with only the primary set: the joined table is that one table
    ctx = build_context(fake_pipeline, roles={"primary": "Nuclei"})
    m = make_measurements(meas_arrays)
    tables = {"Nuclei": build_object_table(ctx, m, "Nuclei")}
    t = join_tables(ctx, m, tables)
    assert t.X.shape[0] == 289
    assert all(v.startswith("Nuclei__") for v in t.var_names)
    assert (t.obs["qc_flag"] == "ok").all()
    assert t.uns["spatialdata_attrs"]["region"] == "Nuclei"
    assert "count_PH3" in t.obs   # Children_PH3_Count lives on the base object here


def test_join_two_image_sets(ctx, meas_arrays):
    m = make_measurements(meas_arrays, image_numbers=(1, 2))
    t = join_tables(ctx, m, _tables(ctx, m))
    assert t.X.shape[0] == 578
    assert list(t.obs["ImageNumber"][[0, 289]]) == [1, 2]
    jn = t.var_names.index("Nuclei__Area")
    numpy.testing.assert_allclose(t.X[289:, jn], meas_arrays["Nuclei"]["AreaShape_Area"].astype("float32"))


# ---- Count_ fallbacks and length mismatches (spec section 5: never a silent NaN / silent skip) ----

def test_missing_count_falls_back_to_longest_array(ctx, meas_arrays, caplog):
    with caplog.at_level(logging.WARNING, logger="scverse_export.assemble"):
        m = make_measurements(meas_arrays, without_counts=("Cells",))
        t = build_object_table(ctx, m, "Cells")
    assert t.X.shape[0] == 289
    assert any("Count_Cells" in r.getMessage() for r in caplog.records)


def test_feature_with_wrong_length_warns_and_stays_nan(ctx, meas_arrays, caplog):
    arrays = {o: dict(f) for o, f in meas_arrays.items()}
    arrays["Cells"]["Intensity_MeanIntensity_DNA"] = arrays["Cells"]["Intensity_MeanIntensity_DNA"][:10]
    with caplog.at_level(logging.WARNING, logger="scverse_export.assemble"):
        t = build_object_table(ctx, make_measurements(arrays), "Cells")
    j = t.var_names.index("Intensity_MeanIntensity__DNA")
    assert t.X.shape[0] == 289 and numpy.isnan(t.X[:, j]).all()
    assert any("Intensity_MeanIntensity_DNA" in r.getMessage() for r in caplog.records)


def test_feature_mismatch_only_affects_its_own_image_set(ctx, meas_arrays, caplog):
    """Regression test for the feature-major rewrite of build_object_table: a length mismatch in
    one image set's copy of a feature must not spill NaN into another image set's rows for that
    same feature -- each feature's series is written into its own offsets[k]:offsets[k+1] slice."""
    m = make_measurements(meas_arrays, image_numbers=(1, 2))
    full = meas_arrays["Cells"]["Intensity_MeanIntensity_DNA"]
    m.next_image_set(2)
    m.add_measurement("Cells", "Intensity_MeanIntensity_DNA", full[:10])
    with caplog.at_level(logging.WARNING, logger="scverse_export.assemble"):
        t = build_object_table(ctx, m, "Cells")
    j = t.var_names.index("Intensity_MeanIntensity__DNA")
    assert t.X.shape[0] == 2 * 289
    numpy.testing.assert_allclose(t.X[:289, j], full.astype("float32"))
    assert numpy.isnan(t.X[289:, j]).all()
    assert any("Intensity_MeanIntensity_DNA" in r.getMessage() for r in caplog.records)


def test_non_numeric_measurements_never_reach_x(ctx, meas_arrays):
    text = Feature(object="Cells", cp_name="Weird_Text", category="Weird", measurement="Text",
                   coltype="varchar(256)")
    number = dataclasses.replace(text, cp_name="Weird_Number", measurement="Number", coltype="float")
    ctx2 = dataclasses.replace(ctx, features=list(ctx.features) + [text, number])
    t = build_object_table(ctx2, make_measurements(meas_arrays), "Cells")
    assert "Weird_Text" not in t.var_names and "Weird_Number" in t.var_names


def test_obs_name_formats_coerced_metadata(ctx, meas_arrays):
    md = {"Metadata_Plate": "P1", "Metadata_Well": "A01", "Metadata_Site": 1.0}
    t = build_object_table(ctx, make_measurements(meas_arrays, md), "Cells")
    assert t.obs_names[0] == "P1_A01_1_1"


# ---- joins that have no Parent_ column (JUMP / FilterObjects pipelines) --------------------------

FILTER_OBJECTS_ENTRY = {
    "num": 99, "name": "FilterObjects", "enabled": True, "settings": {},
    "filter_pairs": [("NucleiIncludingEdges", "Nuclei"), ("CellsIncludingEdges", "Cells")],
}


def _without_parent_nuclei(meas_arrays):
    arrays = {o: dict(f) for o, f in meas_arrays.items()}
    del arrays["Cells"]["Parent_Nuclei"]
    return arrays


def test_join_jump_style_shared_labels(ctx, meas_arrays):
    """FilterObjects relabels the primary and secondary objects together, so Cells carries Parent_CellsIncludingEdges
    and no Parent_Nuclei; the shared label ids are the join key."""
    arrays = _without_parent_nuclei(meas_arrays)
    jump_ctx = dataclasses.replace(ctx, modules=list(ctx.modules) + [FILTER_OBJECTS_ENTRY])
    m = make_measurements(arrays)
    t = join_tables(jump_ctx, m, _tables(jump_ctx, m))
    assert t.X.shape[0] == 289
    assert t.uns["cellprofiler_join"] == {"primary": "shared_label_id", "tertiary": "Parent_Cells"}
    jn = t.var_names.index("Nuclei__Area")
    numpy.testing.assert_allclose(t.X[:, jn], meas_arrays["Nuclei"]["AreaShape_Area"].astype("float32"))


def test_join_missing_parent_raises(ctx, meas_arrays):
    arrays = _without_parent_nuclei(meas_arrays)
    m = make_measurements(arrays)
    with pytest.raises(JoinError, match="Parent_Nuclei"):
        join_tables(ctx, m, _tables(ctx, m))


def test_join_records_parent_columns(ctx, meas_arrays):
    m = make_measurements(meas_arrays)
    t = join_tables(ctx, m, _tables(ctx, m))
    assert t.uns["cellprofiler_join"] == {"primary": "Parent_Nuclei", "tertiary": "Parent_Cells"}


def test_multi_secondary_per_primary_is_flagged(ctx, meas_arrays):
    arrays = {o: dict(f) for o, f in meas_arrays.items()}
    arrays["Cells"]["Parent_Nuclei"] = arrays["Cells"]["Parent_Nuclei"].copy()
    arrays["Cells"]["Parent_Nuclei"][1] = arrays["Cells"]["Parent_Nuclei"][0]
    m = make_measurements(arrays)
    t = join_tables(ctx, m, _tables(ctx, m))
    assert list(t.obs["qc_flag"][:3]) == ["multi_secondary_per_primary", "multi_secondary_per_primary", "ok"]
    assert t.uns["qc_summary"]["multi_secondary_per_primary"] == 2
    assert t.uns["qc_summary"]["primaries_without_secondary"] == 1
    with pytest.raises(PolicyError, match="multi_secondary_per_primary"):
        join_tables(ctx, m, _tables(ctx, m), policy="error")


def test_all_missing_metadata_columns_are_dropped(ctx, meas_arrays):
    # CellProfiler stores NaN for per-image-set metadata that differs between channels
    # (e.g. Metadata_Dye); an all-NaN column carries no information and must not reach obs.
    md = {"Metadata_Plate": "P1", "Metadata_Well": "A01", "Metadata_Site": [1],
          "Metadata_Dye": float("nan")}
    m = make_measurements(meas_arrays, md)
    t = build_object_table(ctx, m, "Cells")
    assert "Metadata_Dye" not in t.obs
    assert t.obs["Metadata_Plate"][0] == "P1"        # partially/fully present columns stay
    assert t.obs_names[0] == "P1_A01_1_1"
