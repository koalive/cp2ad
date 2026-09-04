"""Context + Measurements -> squidpy-shaped tables (plain numpy; no pandas at CellProfiler runtime)."""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

import numpy

from .introspect import Context, IMAGE, file_loaded_objects
from .names import is_extrinsic, to_cpm_names
from .samples import SampleNaming, obs_name, resolve_sample_naming

LOGGER = logging.getLogger(__name__)

VAR_COLUMNS = ("cp_name", "module_num", "module_name", "category", "measurement",
               "channel", "channel2", "other_object", "scale", "coltype", "parsed_by", "region")


@dataclass
class Table:
    X: numpy.ndarray
    obs_names: List[str]
    var_names: List[str]
    obs: Dict[str, numpy.ndarray] = field(default_factory=dict)
    var: Dict[str, numpy.ndarray] = field(default_factory=dict)
    obsm: Dict[str, numpy.ndarray] = field(default_factory=dict)
    uns: Dict[str, Any] = field(default_factory=dict)


def _str_array(values) -> numpy.ndarray:
    return numpy.array(["" if v is None else str(v) for v in values], dtype=object)


def _image_value(m, feat: str, n: int):
    try:
        return m.get_measurement(IMAGE, feat, n)
    except Exception:
        return None


def _object_count(m, obj: str, n: int) -> Optional[int]:
    """Count_<obj> for one image set, or None when it is absent or NaN (the caller falls back)."""
    v = _image_value(m, f"Count_{obj}", n)
    if v is None:
        return None
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return None if numpy.isnan(v) else int(v)


def _series(m, obj: str, feat: str, image_numbers) -> Optional[List[numpy.ndarray]]:
    """One per-object feature over every image set in a single get_measurement call (Measurements
    accepts a sequence of image numbers), as one float array per image set. None = never recorded."""
    if not m.has_feature(obj, feat):
        return None
    vals = m.get_measurement(obj, feat, list(image_numbers))
    if not isinstance(vals, (list, tuple)):
        vals = [vals]
    return [numpy.zeros(0) if v is None else numpy.asarray(v, dtype=numpy.float64).ravel() for v in vals]


NAN_COUNT_BLOCK = 65536


def _count_nan(X: numpy.ndarray) -> numpy.ndarray:
    """isnan(X).sum(axis=1) per row, without materialising a bool copy of the whole matrix (that is
    another 25% of X alive at peak); one row block at a time keeps the temporary small."""
    out = numpy.empty(X.shape[0], dtype=numpy.int32)
    for start in range(0, X.shape[0], NAN_COUNT_BLOCK):
        stop = min(start + NAN_COUNT_BLOCK, X.shape[0])
        out[start:stop] = numpy.isnan(X[start:stop]).sum(axis=1)
    return out


def _is_numeric(coltype: str) -> bool:
    """X is a float matrix, so only CellProfiler's numeric column types may enter it
    (varchar/blob measurements would blow up the cast)."""
    return str(coltype).strip().lower().startswith(("float", "integer"))


def _var_columns(ctx: Context, obj: str):
    """Expand the object's biology-describing features into (var_name, backend, feature) rows for
    X/var; Pearson yields two rows. is_extrinsic excludes position/orientation and identity/linkage
    measurements here; _extrinsic_columns sends those to obs instead, so they cannot bias
    morphological similarity in X."""
    rows = []
    seen = set()
    for f in ctx.features:
        if f.object != obj or not _is_numeric(f.coltype) or is_extrinsic(f):
            continue
        for name, backend in to_cpm_names(f, ctx.channels):
            if name in seen:
                continue
            seen.add(name)
            rows.append((name, backend, f))
    return rows


def _extrinsic_columns(ctx: Context, obj: str):
    """Like _var_columns, but for the measurements is_extrinsic excludes. Same (name, backend,
    feature) row shape, so build_object_table fills both in one pass over the object's features;
    these are destined for obs rather than X."""
    rows = []
    seen = set()
    for f in ctx.features:
        if f.object != obj or not _is_numeric(f.coltype) or not is_extrinsic(f):
            continue
        for name, backend in to_cpm_names(f, ctx.channels):
            if name in seen:
                continue
            seen.add(name)
            rows.append((name, backend, f))
    return rows


def _make_var(rows, obj: str) -> Dict[str, numpy.ndarray]:
    feats = [r[2] for r in rows]
    return {
        "cp_name": _str_array(f.cp_name for f in feats),
        "module_num": numpy.array([-1 if f.module_num is None else f.module_num for f in feats], dtype=numpy.int32),
        "module_name": _str_array(f.module_name for f in feats),
        "category": _str_array(f.category for f in feats),
        "measurement": _str_array(f.measurement for f in feats),
        "channel": _str_array(f.image for f in feats),
        "channel2": _str_array(f.image2 for f in feats),
        "other_object": _str_array(f.other_object for f in feats),
        "scale": _str_array(f.scale for f in feats),
        "coltype": _str_array(f.coltype for f in feats),
        "parsed_by": _str_array(f.parsed_by for f in feats),
        "region": _str_array([obj] * len(feats)),
    }


def metadata_features(ctx: Context, m) -> List[str]:
    """Every Metadata_* Image column this run has. ctx.metadata_tags only lists the *declared*
    ones; CellProfiler also adds Metadata_* Image measurements at runtime (Metadata_Frame,
    Metadata_Series, LoadData/Metadata module tags), so both sources are unioned and obs columns
    and sample keys see runtime-only tags too."""
    return sorted(set(f"Metadata_{t}" for t in ctx.metadata_tags)
                  | set(f for f in m.get_feature_names(IMAGE) if f.startswith("Metadata_")))


def metadata_by_image(m, md_feats: List[str], image_numbers) -> Dict[int, Dict[str, Any]]:
    """{image number: {tag: value}} for the given tags, which is what deciding whether those tags
    identify the image sets needs before any row is built."""
    return {int(n): {tag: _image_value(m, tag, n) for tag in md_feats} for n in image_numbers}


def resolve_naming(ctx: Context, m, requested: Optional[Sequence[str]] = None) -> SampleNaming:
    """The sample-key naming for this run, detected or as requested, verified against the real
    metadata values. Resolved once per export and passed to build_object_table, so every object's
    table names its rows the same way."""
    md_feats = metadata_features(ctx, m)
    image_numbers = [int(n) for n in m.get_image_numbers()]
    return resolve_sample_naming(md_feats, metadata_by_image(m, md_feats, image_numbers), requested)


def build_object_table(ctx: Context, m, obj: str, naming: Optional[SampleNaming] = None) -> Table:
    rows = _var_columns(ctx, obj)
    extrinsic_rows = _extrinsic_columns(ctx, obj)
    var_names = [r[0] for r in rows]
    extrinsic_names = [r[0] for r in extrinsic_rows]
    image_numbers = [int(n) for n in m.get_image_numbers()]
    md_feats = metadata_features(ctx, m)
    if naming is None:
        naming = resolve_naming(ctx, m)

    wanted = list(dict.fromkeys([f.cp_name for _, _, f in rows] + [f.cp_name for _, _, f in extrinsic_rows] +
                                 ["Location_Center_X", "Location_Center_Y"]))
    cols_by_feat: Dict[str, List[int]] = {}
    for j, (_, _, f) in enumerate(rows):
        cols_by_feat.setdefault(f.cp_name, []).append(j)
    extrinsic_cols_by_feat: Dict[str, List[int]] = {}
    for j, (_, _, f) in enumerate(extrinsic_rows):
        extrinsic_cols_by_feat.setdefault(f.cp_name, []).append(j)

    # ---- row layout: Count_<obj> per image set, falling back (only when it is missing) to the
    # longest wanted feature's array length for that one image set -- never reading every feature
    # for every image set just to find this out. ----
    counts: List[int] = []
    for n in image_numbers:
        count = _object_count(m, obj, n)
        if count is None:
            longest = 0
            for feat in wanted:
                if not m.has_feature(obj, feat):
                    continue
                v = m.get_measurement(obj, feat, n)
                longest = max(longest, 0 if v is None else numpy.atleast_1d(v).ravel().shape[0])
            if longest:
                LOGGER.warning("ExportToAnnData: no finite Count_%s in image set %d; falling back to the "
                               "longest %s measurement array (%d rows).", obj, n, obj, longest)
            count = longest
        counts.append(count)
    offsets = numpy.concatenate([[0], numpy.cumsum(counts)]).astype(numpy.int64)
    n_rows = int(offsets[-1])

    # ---- fill X feature-major: one get_measurement call per feature across every image set, one
    # float64 series alive at a time, written straight into the float32 X rows it belongs to, then
    # dropped -- peak memory stays ~= X (float32) + one feature's float64 series.
    #
    # extrinsic_rows (position/orientation, identifiers, linkage) get the same treatment, filled
    # into obs_extrinsic instead of X. None of them describe the object's own biology, so leaving
    # them in X would bias similarity on position, rotation, or an arbitrary CellProfiler label. ----
    X = numpy.full((n_rows, len(rows)), numpy.nan, dtype=numpy.float32)
    obs_extrinsic = numpy.full((n_rows, len(extrinsic_rows)), numpy.nan, dtype=numpy.float32)
    spatial_arr = numpy.full((n_rows, 2), numpy.nan)
    for feat in wanted:
        series = _series(m, obj, feat, image_numbers)
        if series is None:
            continue
        cols = cols_by_feat.get(feat, [])
        scols = extrinsic_cols_by_feat.get(feat, [])
        for k, n in enumerate(image_numbers):
            count = counts[k]
            if count == 0:
                continue
            a = series[k]
            if a.shape[0] != count:
                LOGGER.warning("ExportToAnnData: %s.%s has %d values in image set %d but the image set has "
                               "%d %s; that column stays NaN.", obj, feat, a.shape[0], n, count, obj)
                continue
            start, stop = int(offsets[k]), int(offsets[k + 1])
            for j in cols:
                X[start:stop, j] = a
            for j in scols:
                obs_extrinsic[start:stop, j] = a
            if feat == "Location_Center_X":
                spatial_arr[start:stop, 0] = a
            elif feat == "Location_Center_Y":
                spatial_arr[start:stop, 1] = a
        del series

    obs_names, label_id, imgnum, md_cols = [], [], [], {k: [] for k in md_feats}
    for k, n in enumerate(image_numbers):
        count = counts[k]
        if count == 0:
            continue
        md = {key: _image_value(m, key, n) for key in md_feats}
        for key in md_feats:
            md_cols[key].extend([md[key]] * count)
        labels = numpy.arange(1, count + 1)
        obs_names.extend(obs_name(md, n, int(l), naming) for l in labels)
        label_id.append(labels)
        imgnum.append(numpy.full(count, n, dtype=numpy.int32))

    if label_id:
        label_arr, img_arr = numpy.concatenate(label_id), numpy.concatenate(imgnum)
    else:
        label_arr = img_arr = numpy.zeros(0, dtype=numpy.int32)
    miss_arr = _count_nan(X)
    obs = {"region": _str_array([obj] * len(obs_names)), "label_id": label_arr.astype(numpy.int32),
           "ImageNumber": img_arr, "n_missing_features": miss_arr}
    for name, col in zip(extrinsic_names, obs_extrinsic.T):
        obs[name] = col
    for key in md_feats:
        vals = md_cols[key]
        if vals and all(v is None or (isinstance(v, (float, numpy.floating)) and numpy.isnan(v)) for v in vals):
            # CellProfiler stores NaN for image-set metadata that differs between channels
            # (e.g. Metadata_Dye, Metadata_FileLocation) -- an all-missing column carries nothing.
            LOGGER.debug("build_object_table: dropping all-missing Image metadata column %s", key)
            continue
        numeric = all(isinstance(v, (int, float, numpy.integer, numpy.floating)) for v in vals) and vals
        obs[key] = numpy.asarray(vals, dtype=float) if numeric else _str_array(vals)
        if numeric and numpy.all(numpy.mod(obs[key], 1) == 0):
            obs[key] = obs[key].astype(numpy.int64)
    return Table(X=X, obs_names=obs_names, var_names=var_names, obs=obs, var=_make_var(rows, obj),
                 obsm={"spatial": spatial_arr},
                 uns={"spatialdata_attrs": {"region": obj, "region_key": "region", "instance_key": "label_id"}})


POLICIES = ("flag", "drop", "error")
FLAG_OK, FLAG_NO_PRIMARY = "ok", "no_primary"
FLAG_NO_TERTIARY, FLAG_MULTI_TERTIARY = "no_tertiary", "multi_tertiary"
FLAG_MULTI_SECONDARY = "multi_secondary_per_primary"


class PolicyError(ValueError):
    pass


class JoinError(ValueError):
    pass


def _row_index(t: Table) -> Dict[tuple, int]:
    return {(int(n), int(l)): i for i, (n, l) in enumerate(zip(t.obs["ImageNumber"], t.obs["label_id"]))}


def _lookup(m, obj: str, feat: str, image_numbers, default=0) -> Optional[numpy.ndarray]:
    """A per-object integer column (Parent_*/Children_*) concatenated over image sets, in the row
    order build_object_table produced. None (never a short array) when the column does not exist."""
    series = _series(m, obj, feat, image_numbers)
    if series is None:
        return None
    out = numpy.concatenate(series) if series else numpy.zeros(0)
    return numpy.where(numpy.isnan(out), default, out)


def _relabelled_together(ctx: Context, a: str, b: str) -> bool:
    """True when one FilterObjects module emitted both objects: CellProfiler relabels the objects it
    filters together, so they share their label ids and no Parent_ column is written."""
    for mod in ctx.modules:
        if mod.get("name") == "FilterObjects":
            outputs = {dst for _, dst in mod.get("filter_pairs", [])}
            if a in outputs and b in outputs:
                return True
    return False


def _join_column(ctx: Context, m, child: str, parent: str, table: Table, image_numbers) -> tuple:
    """(values to look `parent` rows up by, join description) for each row of `child`'s table."""
    feat = f"Parent_{parent}"
    col = _lookup(m, child, feat, image_numbers)
    if col is None:
        if not _relabelled_together(ctx, child, parent):
            raise JoinError(
                f"{child} has no {feat} column and no FilterObjects module relabelled {child} and "
                f"{parent} together, so {child} rows cannot be matched to {parent} rows. Set "
                f"primary/secondary/tertiary roles explicitly.")
        # shared label ids: the objects came out of the same FilterObjects relabelling
        return table.obs["label_id"].astype(numpy.float64), "shared_label_id"
    if col.shape[0] != table.X.shape[0]:
        raise JoinError(f"{child}.{feat} has {col.shape[0]} values but {child} has {table.X.shape[0]} rows; "
                        f"the pipeline's Count_{child} and its measurement arrays disagree. Set "
                        f"primary/secondary/tertiary roles explicitly.")
    return col, feat


def join_tables(ctx: Context, m, tables: Dict[str, Table], policy: str = "flag") -> Table:
    if policy not in POLICIES:
        raise ValueError(f"policy must be one of {POLICIES}")
    roles = ctx.roles
    base_obj = roles.get("secondary", roles.get("primary"))
    base = tables[base_obj]
    n = base.X.shape[0]
    image_numbers = [int(x) for x in m.get_image_numbers()]
    order = [base_obj] + [roles[r] for r in ("primary", "tertiary") if r in roles and roles[r] != base_obj]

    # row of each compartment for each base row (-1 = missing)
    idx = {base_obj: numpy.arange(n)}
    flags = numpy.array([FLAG_OK] * n, dtype=object)
    join_keys: Dict[str, str] = {}
    primaries_without_secondary = 0
    if "primary" in roles and roles["primary"] != base_obj:
        prim = roles["primary"]; tp = tables[prim]; prim_rows = _row_index(tp)
        parent, join_keys["primary"] = _join_column(ctx, m, base_obj, prim, base, image_numbers)
        idx[prim] = numpy.array([prim_rows.get((int(img), int(p)), -1)
                                 for img, p in zip(base.obs["ImageNumber"], parent)])
        flags[idx[prim] < 0] = FLAG_NO_PRIMARY
        used = set(idx[prim][idx[prim] >= 0].tolist())
        primaries_without_secondary = tp.X.shape[0] - len(used)
        # primary side of the 1:1:1 check: one primary claimed by several secondaries (spec section 5)
        if n:
            _, inverse, counts = numpy.unique(idx[prim], return_inverse=True, return_counts=True)
            shared = (counts[inverse] > 1) & (idx[prim] >= 0) & (flags == FLAG_OK)
            flags[shared] = FLAG_MULTI_SECONDARY
    if "tertiary" in roles:
        tert = roles["tertiary"]; tc = tables[tert]
        parent, join_keys["tertiary"] = _join_column(ctx, m, tert, base_obj, tc, image_numbers)
        by_parent: Dict[tuple, list] = {}
        for i, (img, p) in enumerate(zip(tc.obs["ImageNumber"], parent)):
            by_parent.setdefault((int(img), int(p)), []).append(i)
        # cardinality comes from the secondary's own Children_<tertiary>_Count, not from re-deriving it by
        # counting Parent_<base_obj> hits: the two normally agree, but the count is the authoritative
        # CellProfiler-reported relationship and is what a pipeline-level 1:1:1 violation would perturb.
        counts = _lookup(m, base_obj, f"Children_{tert}_Count", image_numbers)
        if counts is not None and counts.shape[0] != n:
            raise JoinError(f"{base_obj}.Children_{tert}_Count has {counts.shape[0]} values but {base_obj} "
                            f"has {n} rows.")
        rows_c = []
        for i, (img, l) in enumerate(zip(base.obs["ImageNumber"], base.obs["label_id"])):
            hits = by_parent.get((int(img), int(l)), [])
            c = int(counts[i]) if counts is not None else len(hits)
            if c == 0:
                rows_c.append(-1)
                if flags[i] == FLAG_OK:
                    flags[i] = FLAG_NO_TERTIARY
            else:
                rows_c.append(hits[0] if hits else -1)
                if c > 1 and flags[i] == FLAG_OK:
                    flags[i] = FLAG_MULTI_TERTIARY
        idx[tert] = numpy.array(rows_c)

    summary = {k: int(v) for k, v in zip(*numpy.unique(flags, return_counts=True))}
    summary["primaries_without_secondary"] = int(primaries_without_secondary)
    keep = numpy.ones(n, dtype=bool)
    if policy == "error" and (flags != FLAG_OK).any():
        raise PolicyError("1:1:1 violations: " + ", ".join(f"{k}={v}" for k, v in summary.items() if k != FLAG_OK))
    if policy == "drop":
        keep = flags == FLAG_OK

    blocks, var_names, var = [], [], {k: [] for k in VAR_COLUMNS}
    for obj in order:
        t = tables[obj]; rows = idx[obj]
        block = numpy.full((n, t.X.shape[1]), numpy.nan, dtype=numpy.float32)
        ok = rows >= 0
        block[ok] = t.X[rows[ok]]
        blocks.append(block)
        var_names.extend(f"{obj}__{v}" for v in t.var_names)
        for k in VAR_COLUMNS:
            var[k].extend(t.var[k].tolist())
    X = numpy.hstack(blocks)[keep]
    # base's own extrinsic columns are re-added below with the same "{obj}__" prefix as every other
    # compartment, so they are excluded here rather than kept twice under two names.
    base_extrinsic_names = {r[0] for r in _extrinsic_columns(ctx, base_obj)}
    obs = {k: v[keep] for k, v in base.obs.items() if k not in base_extrinsic_names}
    obs["qc_flag"] = flags[keep]
    obs["n_missing_features"] = _count_nan(X)
    for obj in order:
        t = tables[obj]; rows = idx[obj]; ok = rows >= 0
        for name in (r[0] for r in _extrinsic_columns(ctx, obj)):
            mapped = numpy.full(n, numpy.nan, dtype=numpy.float32)
            mapped[ok] = t.obs[name][rows[ok]]
            obs[f"{obj}__{name}"] = mapped[keep]
    role_objs = set(roles.values())
    for f in ctx.features:
        if f.category == "Children" and f.measurement.endswith("_Count") and f.object in role_objs:
            child = f.measurement[: -len("_Count")]
            if child in role_objs:
                continue
            src = f.object
            counts = _lookup(m, src, f.cp_name, image_numbers)
            if counts is None:
                continue
            if src == base_obj:
                col = counts if counts.shape[0] == n else None
            else:  # child counts live on the primary/tertiary row: map through idx
                rows = idx.get(src)
                col = numpy.where(rows >= 0, counts[numpy.clip(rows, 0, None)], 0) \
                    if rows is not None and counts.shape[0] == tables[src].X.shape[0] else None
            if col is not None:
                obs[f"count_{child}"] = col[keep].astype(numpy.int32)
    var_arrays = {k: numpy.array(v, dtype=object) if k not in ("module_num",) else numpy.array(v, dtype=numpy.int32)
                  for k, v in var.items()}
    uns = {"spatialdata_attrs": {"region": base_obj, "region_key": "region", "instance_key": "label_id"},
           "qc_summary": summary, "cellprofiler_join": join_keys}
    return Table(X=X, obs_names=[o for o, k in zip(base.obs_names, keep) if k], var_names=var_names,
                 obs=obs, var=var_arrays, obsm={"spatial": base.obsm["spatial"][keep]}, uns=uns)


def provenance(ctx: Context, m, exporter_settings: Dict[str, Any],
               naming: Optional[SampleNaming] = None) -> Dict[str, Any]:
    image_numbers = [int(x) for x in m.get_image_numbers()]
    try:
        runtime_image_feats = set(m.get_feature_names(IMAGE))
    except Exception:
        runtime_image_feats = set()
    # ctx.features only lists *declared* Image measurements; CellProfiler also adds runtime-only
    # Image measurements (Metadata_Frame, Metadata_Series, LoadData/Metadata module tags) — union
    # both sources, mirroring build_object_table's md_feats handling.
    image_feats = sorted(set(f.cp_name for f in ctx.features if f.object == IMAGE) | runtime_image_feats)
    image = {"ImageNumber": image_numbers}
    for feat in image_feats:
        if not m.has_feature(IMAGE, feat):
            continue
        try:  # one call for every image set; falls back to one call each if a feature dislikes it
            vals = list(m.get_measurement(IMAGE, feat, image_numbers))
        except Exception:
            vals = [_image_value(m, feat, n) for n in image_numbers]
        image[feat] = [v.item() if isinstance(v, numpy.generic) else v for v in vals]
    experiment = {}
    try:
        experiment_feats = m.get_feature_names("Experiment")
    except Exception:
        experiment_feats = []
    for feat in experiment_feats:
        try:
            v = m.get_experiment_measurement(feat)
            experiment[feat] = v.item() if isinstance(v, numpy.generic) else v
        except Exception:
            pass
    relationships = []
    if hasattr(m, "get_relationship_groups"):
        for k in m.get_relationship_groups():
            r = m.get_relationships(k.module_number, k.relationship, k.object_name1, k.object_name2)
            relationships.append({"module_num": int(k.module_number), "relationship": k.relationship,
                                  "object1": k.object_name1, "object2": k.object_name2, "n": int(len(r))})
    from_files = set(file_loaded_objects(ctx))
    return {
        "version": experiment.get("CellProfiler_Version"),
        "run_timestamp": experiment.get("Run_Timestamp"),
        "pipeline_text": experiment.get("Pipeline_Pipeline"),
        # "settings" (a dict that keeps only the first occurrence of a repeated setting text, e.g.
        # one of several Texture scales) is a Context-internal aid for role detection; drop it here
        # and keep "setting_values" (the complete, ordered (text, value) pairs) as the exported record.
        "modules": [{k: v for k, v in mod.items() if k != "settings"} for mod in ctx.modules],
        "channels": list(ctx.channels),
        "objects": {k: {"module_num": v.module_num, "module_name": v.module_name, "role": v.role,
                        "source": "file" if k in from_files else "pipeline"}
                    for k, v in ctx.objects.items()},
        "roles": dict(ctx.roles),
        "role_detection": dict(ctx.role_note),
        # How obs names were built, so a reader can tell what a row name means without guessing,
        # and can see when the image number had to stand in for missing metadata.
        "sample_naming": ({"tags": list(naming.tags), "parts": list(naming.parts),
                           "with_image_number": bool(naming.with_image_number),
                           "mode": naming.mode, "note": naming.note}
                          if naming is not None else {}),
        "relationships": relationships,
        "image": image,
        "experiment": experiment,
        "exporter": {"name": "ExportToAnnData", "settings": dict(exporter_settings)},
    }
