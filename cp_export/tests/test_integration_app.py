"""Runs the real CellProfiler app headless with the plugin. Skipped when the app binary is absent."""
import os
import subprocess

import numpy
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PLUGIN_DIR = os.path.dirname(HERE)
SANDBOX = os.path.join(os.path.dirname(PLUGIN_DIR), "plugin_sandbox")
CP = "/Applications/CellProfiler.app/Contents/MacOS/cp"
pytestmark = pytest.mark.skipif(not os.path.exists(CP), reason="CellProfiler app not installed")


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    out = tmp_path_factory.mktemp("adata")
    cmd = [CP, "-c", "-r", "-p", os.path.join(HERE, "pipelines", "examplehuman_adata.cppipe"),
           "-i", os.path.join(SANDBOX, "probe", "img"), "-o", str(out),
           f"--plugins-directory={PLUGIN_DIR}", "--log-level=INFO"]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    log = res.stdout + res.stderr
    assert res.returncode == 0, log[-3000:]
    assert "Traceback" not in log, log[-3000:]
    return out


def test_joined_file(run, cpm_columns):
    import anndata
    a = anndata.read_h5ad(str(run / "cellprofiler.h5ad"))
    assert a.n_obs == 289
    assert set(a.var["region"]) == {"Nuclei", "Cells", "Cytoplasm"}
    assert (a.obs["qc_flag"] == "ok").all()
    assert list(a.obs_names[:2]) == ["img1_1", "img1_2"]
    cols = set(cpm_columns)
    matches = 0
    for name, region in zip(a.var_names, a.var["region"]):
        assert name.startswith(region + "__")
        matches += name[len(region) + 2:] in cols
    assert matches > 0
    assert a.uns["cellprofiler"]["roles"]["secondary"] == "Cells"
    assert a.uns["cellprofiler"]["role_detection"]["mode"] == "automatic"
    assert "IdentifySecondaryObjects" in [m["name"] for m in a.uns["cellprofiler"]["modules"].values()]
    assert "Version:5" in a.uns["cellprofiler"]["pipeline_text"]
    assert a.obsm["spatial"].shape == (289, 2)
    assert "count_PH3" in a.obs
    assert dict(a.uns["cellprofiler_join"]) == {"primary": "Parent_Nuclei", "tertiary": "Parent_Cells"}
    texture = next(m for m in a.uns["cellprofiler"]["modules"].values() if m["name"] == "MeasureTexture")
    scales = [v for t, v in texture["setting_values"].tolist() if t == "Texture scale to measure"]
    assert scales == ["3", "5", "10"]        # repeated setting texts survive in setting_values
    assert a.uns["cellprofiler"]["objects"]["Cells"]["source"] == "pipeline"


def test_per_object_files_match_cpm_columns(run, cpm_columns):
    import anndata
    cols = set(cpm_columns)
    for obj in ("Nuclei", "Cells", "Cytoplasm"):
        a = anndata.read_h5ad(str(run / f"cellprofiler_{obj}.h5ad"))
        assert a.n_obs == 289 and a.uns["spatialdata_attrs"]["region"] == obj
        cpm = [n for n in a.var_names if n in cols]
        assert cpm
        assert numpy.isfinite(a[:, "Area"].X).all()
