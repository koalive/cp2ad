import re
import pytest
from scverse_export.names import Feature, is_extrinsic, split_pair, to_cpm_names

CH = ["DNA", "PH3", "cellbody"]


def F(cp_name, category, measurement, image=None, image2=None, scale=None, other_object=None, obj="Cells"):
    return Feature(object=obj, cp_name=cp_name, category=category, measurement=measurement,
                   image=image, image2=image2, scale=scale, other_object=other_object)


@pytest.mark.parametrize("feat,expected", [
    (F("AreaShape_Area", "AreaShape", "Area"), [("Area", "cp_measure")]),
    (F("AreaShape_Zernike_3_1", "AreaShape", "Zernike_3_1"), [("Zernike_3_1", "cp_measure")]),
    (F("AreaShape_Center_X", "AreaShape", "Center_X"), [("Center_X", "cp_measure")]),
    (F("Location_Center_X", "Location", "Center_X"), [("Center_X", "cp_measure")]),
    (F("Location_Center_Z", "Location", "Center_Z"), [("Location_Center_Z", "cellprofiler")]),
    (F("Intensity_MeanIntensity_DNA", "Intensity", "MeanIntensity", image="DNA"),
     [("Intensity_MeanIntensity__DNA", "cp_measure")]),
    (F("Intensity_MeanIntensityEdge_PH3", "Intensity", "MeanIntensityEdge", image="PH3"),
     [("Intensity_MeanIntensityEdge__PH3", "cp_measure")]),
    (F("Location_CenterMassIntensity_X_DNA", "Location", "CenterMassIntensity_X", image="DNA"),
     [("Location_CenterMassIntensity_X__DNA", "cp_measure")]),
    (F("Texture_Contrast_DNA_3_00_256", "Texture", "Contrast", image="DNA", scale="3_00_256"),
     [("Contrast_3_00_256__DNA", "cp_measure")]),
    (F("Correlation_Correlation_DNA_PH3", "Correlation", "Correlation", image="DNA_PH3"),
     [("Correlation_Pearson__DNA__PH3", "cp_measure"), ("Correlation_Pearson__PH3__DNA", "cp_measure")]),
    (F("Correlation_Slope_DNA_PH3", "Correlation", "Slope", image="DNA_PH3"),
     [("Correlation_Slope__DNA__PH3", "cp_measure")]),
    (F("Correlation_Manders_DNA_PH3", "Correlation", "Manders", image="DNA_PH3"),
     [("Correlation_Manders_1__DNA__PH3", "cp_measure")]),
    (F("Correlation_Manders_PH3_DNA", "Correlation", "Manders", image="PH3_DNA"),
     [("Correlation_Manders_2__DNA__PH3", "cp_measure")]),
    (F("Correlation_RWC_cellbody_DNA", "Correlation", "RWC", image="cellbody_DNA"),
     [("Correlation_RWC_2__DNA__cellbody", "cp_measure")]),
    (F("Correlation_Costes_DNA_cellbody", "Correlation", "Costes", image="DNA_cellbody"),
     [("Correlation_Costes_1__DNA__cellbody", "cp_measure")]),
    (F("Correlation_K_DNA_PH3", "Correlation", "K", image="DNA_PH3"), [("Correlation_K_DNA_PH3", "cellprofiler")]),
    (F("Correlation_Overlap_DNA_PH3", "Correlation", "Overlap", image="DNA_PH3"),
     [("Correlation_Overlap_DNA_PH3", "cellprofiler")]),
    (F("RadialDistribution_FracAtD_DNA_1of4", "RadialDistribution", "FracAtD", image="DNA", scale="1of4"),
     [("RadialDistribution_FracAtD_1of4__DNA", "cp_measure")]),
    (F("RadialDistribution_ZernikeMagnitude_DNA_2_0", "RadialDistribution", "ZernikeMagnitude", image="DNA", scale="2_0"),
     [("RadialDistribution_ZernikeMagnitude_2_0__DNA", "cp_measure")]),
    (F("Granularity_1_DNA", "Granularity", "1", image="DNA"), [("Granularity_1__DNA", "cp_measure")]),
    (F("Neighbors_NumberOfNeighbors_5", "Neighbors", "NumberOfNeighbors", scale="5"),
     [("Neighbors_NumberOfNeighbors_5", "cellprofiler")]),
    (F("Parent_Nuclei", "Parent", "Nuclei"), [("Parent_Nuclei", "cellprofiler")]),
    (F("Children_Cytoplasm_Count", "Children", "Cytoplasm_Count"), [("Children_Cytoplasm_Count", "cellprofiler")]),
    (F("Number_Object_Number", "Number", "Object_Number"), [("Number_Object_Number", "cellprofiler")]),
])
def test_rewrite_rules(feat, expected):
    assert to_cpm_names(feat, CH) == expected


@pytest.mark.parametrize("feat,expected", [
    # position/orientation -> obs
    (F("Location_Center_X", "Location", "Center_X"), True),
    (F("Location_Center_Z", "Location", "Center_Z"), True),
    (F("Location_CenterMassIntensity_X_DNA", "Location", "CenterMassIntensity_X", image="DNA"), True),
    (F("AreaShape_Center_X", "AreaShape", "Center_X"), True),
    (F("AreaShape_Orientation", "AreaShape", "Orientation"), True),
    (F("AreaShape_BoundingBoxMinimum_X", "AreaShape", "BoundingBoxMinimum_X"), True),
    (F("AreaShape_BoundingBoxMaximum_Y", "AreaShape", "BoundingBoxMaximum_Y"), True),
    # identity/linkage -> obs
    (F("Number_Object_Number", "Number", "Object_Number"), True),
    (F("Parent_Nuclei", "Parent", "Nuclei"), True),
    (F("Children_Cytoplasm_Count", "Children", "Cytoplasm_Count"), True),
    (F("Neighbors_FirstClosestObjectNumber_5", "Neighbors", "FirstClosestObjectNumber"), True),
    (F("Neighbors_SecondClosestObjectNumber_5", "Neighbors", "SecondClosestObjectNumber"), True),
    # actual biology -> stays in X, including the Neighbors measurements treated as biology here
    (F("AreaShape_Area", "AreaShape", "Area"), False),
    (F("AreaShape_BoundingBoxArea", "AreaShape", "BoundingBoxArea"), False),
    (F("Intensity_MeanIntensity_DNA", "Intensity", "MeanIntensity", image="DNA"), False),
    (F("Neighbors_NumberOfNeighbors_5", "Neighbors", "NumberOfNeighbors"), False),
    (F("Neighbors_PercentTouching_5", "Neighbors", "PercentTouching"), False),
    (F("Neighbors_FirstClosestDistance_5", "Neighbors", "FirstClosestDistance"), False),
    (F("Neighbors_AngleBetweenNeighbors_5", "Neighbors", "AngleBetweenNeighbors"), False),
])
def test_is_extrinsic(feat, expected):
    assert is_extrinsic(feat) is expected


def test_split_pair_handles_underscores_in_channel_names():
    assert split_pair("DNA_PH3", CH) == ("DNA", "PH3")
    assert split_pair("Mito_Tubeness_DNA", ["DNA", "Mito_Tubeness"]) == ("Mito_Tubeness", "DNA")
    assert split_pair("DNA_Mito_Tubeness", ["DNA", "Mito_Tubeness"]) == ("DNA", "Mito_Tubeness")
    assert split_pair("DNA_XX", CH) is None


def test_every_value_verified_pair_is_reproduced_by_rules(mapping_cells, cpm_columns):
    """mapping_cells.json: cp_measure column -> list of CP features with identical values.
    Spurious all-zero matches (Z coordinates, empty granularity, moments) are excluded by the
    cp_measure-side filter below; every remaining CP feature must rewrite to that cp_measure column."""
    spurious_prefixes = ("CentralMoment", "SpatialMoment", "NormalizedMoment", "HuMoment", "InertiaTensor",
                         "Location_CenterMassIntensity_Z", "Location_MaxIntensity_Z", "FilledArea",
                         "Correlation_Costes", "Granularity")
    checked = 0
    for cpm_col, cp_feats in mapping_cells.items():
        if cpm_col.startswith(spurious_prefixes):
            continue
        for cp in cp_feats:
            if cp.startswith(("Location_Center_Z", "Location_CenterMassIntensity_Z", "Location_MaxIntensity_Z",
                              "Children_", "Granularity")):
                continue
            feat = feature_from_cp_name(cp)
            names = [n for n, backend in to_cpm_names(feat, CH) if backend == "cp_measure"]
            assert cpm_col in names, (cpm_col, cp, names)
            checked += 1
    assert checked >= 250


def feature_from_cp_name(cp):
    """Test-only parser for the probe's known CP names (the real code gets structure from the module API)."""
    cat, rest = cp.split("_", 1)
    if cat == "AreaShape":
        return F(cp, cat, rest)
    if cat == "Intensity":
        meas, ch = rest.rsplit("_", 1)
        return F(cp, cat, meas, image=ch)
    if cat == "Location":
        if rest in ("Center_X", "Center_Y", "Center_Z"):
            return F(cp, cat, rest)
        meas, ch = rest.rsplit("_", 1)
        return F(cp, cat, meas, image=ch)
    if cat == "Texture":
        m = re.match(r"(\w+?)_(DNA|PH3|cellbody)_(\d+_\d\d_\d+)$", rest)
        return F(cp, cat, m.group(1), image=m.group(2), scale=m.group(3))
    if cat == "Correlation":
        meas, pair = rest.split("_", 1)
        return F(cp, cat, meas, image=pair)
    if cat == "RadialDistribution":
        m = re.match(r"(\w+?)_(DNA|PH3|cellbody)_(.+)$", rest)
        return F(cp, cat, m.group(1), image=m.group(2), scale=m.group(3))
    return F(cp, cat, rest)
