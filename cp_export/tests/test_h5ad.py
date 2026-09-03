import numpy
import anndata

from cpexport.assemble import Table
from cpexport.h5ad import write_h5ad


def _table():
    return Table(
        X=numpy.array([[1.0, numpy.nan], [3.0, 4.0]], dtype=numpy.float32),
        obs_names=["P1_A01_1_1", "P1_A01_1_2"], var_names=["Area", "Intensity_MeanIntensity__DNA"],
        obs={"region": numpy.array(["Cells", "Cells"], dtype=object), "label_id": numpy.array([1, 2], dtype=numpy.int32),
             "Metadata_Site": numpy.array([1, 1]), "qc_flag": numpy.array(["ok", "no_primary"], dtype=object)},
        var={"cp_name": numpy.array(["AreaShape_Area", "Intensity_MeanIntensity_DNA"], dtype=object),
             "module_num": numpy.array([11, 10], dtype=numpy.int32),
             "channel": numpy.array(["", "DNA"], dtype=object),
             "note": numpy.array([b"a", b"b"], dtype="S1")},
        obsm={"spatial": numpy.array([[1.5, 2.5], [3.5, 4.5]])},
        uns={"spatialdata_attrs": {"region": "Cells", "region_key": "region", "instance_key": "label_id"},
             "cellprofiler": {"version": "4.2.8", "channels": ["DNA", "PH3"],
                              "modules": [{"num": 1, "name": "Images"},
                                          {"num": 2, "name": "X",
                                           "params": {"How to pick primary / secondary / tertiary objects": "Automatic"},
                                           "setting_values": [("Texture scale to measure", "3"),
                                                              ("Texture scale to measure", "5")]}],
                              "image": {"ImageNumber": [1, 2], "Count_Cells": [2.0], "Threshold_FinalThreshold_Nuclei": [0.1, None]}, "pipeline_text": None,
                              "qc": {"ok": 1, "no_primary": 1}}},
    )


def test_roundtrip_with_anndata(tmp_path):
    path = str(tmp_path / "t.h5ad")
    write_h5ad(_table(), path)
    a = anndata.read_h5ad(path)
    assert a.shape == (2, 2) and a.X.dtype == numpy.float32 and numpy.isnan(a.X[0, 1])
    assert list(a.obs_names) == ["P1_A01_1_1", "P1_A01_1_2"] and list(a.var_names) == ["Area", "Intensity_MeanIntensity__DNA"]
    assert list(a.obs["region"]) == ["Cells", "Cells"] and a.obs["label_id"].tolist() == [1, 2]
    assert list(a.obs["qc_flag"]) == ["ok", "no_primary"]
    assert a.var["cp_name"].tolist() == ["AreaShape_Area", "Intensity_MeanIntensity_DNA"]
    assert a.var["module_num"].tolist() == [11, 10] and a.var["channel"].tolist() == ["", "DNA"]
    assert a.var["note"].tolist() == ["a", "b"]
    assert a.obsm["spatial"].shape == (2, 2) and a.obsm["spatial"][1, 0] == 3.5
    assert a.uns["spatialdata_attrs"]["instance_key"] == "label_id"
    assert list(a.uns["cellprofiler"]["channels"]) == ["DNA", "PH3"]
    assert a.uns["cellprofiler"]["modules"]["0"]["name"] == "Images"      # lists of dicts are stored as "0","1",... groups
    assert (a.uns["cellprofiler"]["modules"]["1"]["params"]["How to pick primary | secondary | tertiary objects"]
            == "Automatic")  # '/' in a dict key is escaped to '|' (HDF5 path separator)
    assert a.uns["cellprofiler"]["modules"]["1"]["setting_values"].tolist() == [
        ["Texture scale to measure", "3"], ["Texture scale to measure", "5"]]  # pairs stay pairs
    assert list(a.uns["cellprofiler"]["image"]["ImageNumber"]) == [1, 2]
    assert numpy.isclose(a.uns["cellprofiler"]["image"]["Threshold_FinalThreshold_Nuclei"][0], 0.1)
    assert numpy.isnan(a.uns["cellprofiler"]["image"]["Threshold_FinalThreshold_Nuclei"][1])
    assert a.uns["cellprofiler"]["pipeline_text"] == ""
    assert a.uns["cellprofiler"]["qc"]["ok"] == 1


def test_empty_table(tmp_path):
    t = Table(X=numpy.zeros((0, 1), dtype=numpy.float32), obs_names=[], var_names=["Area"],
              obs={"region": numpy.array([], dtype=object)}, var={"cp_name": numpy.array(["AreaShape_Area"], dtype=object)})
    path = str(tmp_path / "e.h5ad")
    write_h5ad(t, path)
    assert anndata.read_h5ad(path).shape == (0, 1)
