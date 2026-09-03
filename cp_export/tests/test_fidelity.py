import numpy
from cpexport.assemble import build_object_table, join_tables
from cpexport.introspect import build_context
from test_assemble import make_measurements


def test_renamed_columns_carry_cp_values(fake_pipeline, meas_arrays, mapping_cells):
    ctx = build_context(fake_pipeline)
    m = make_measurements(meas_arrays)
    tables = {obj: build_object_table(ctx, m, obj) for obj in ctx.roles.values()}
    t = join_tables(ctx, m, tables)
    checked = 0
    for cpm_col, cp_feats in mapping_cells.items():
        name = f"Cells__{cpm_col}"
        if name in t.var_names:
            j = t.var_names.index(name)
            cp = t.var["cp_name"][j]
            if cp not in cp_feats:
                continue
            numpy.testing.assert_allclose(t.X[:, j], meas_arrays["Cells"][cp].astype("float32"), rtol=1e-6)
            checked += 1
        elif name in t.obs:
            # Location/orientation columns (spec: never in X, biases similarity on position/rotation)
            # land in obs instead, under the same "{object}__{name}" prefix; check fidelity there.
            for cp in cp_feats:
                if cp in meas_arrays["Cells"]:
                    numpy.testing.assert_allclose(t.obs[name], meas_arrays["Cells"][cp].astype("float32"),
                                                  rtol=1e-6)
                    checked += 1
                    break
    assert checked >= 250
