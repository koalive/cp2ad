"""Read everything an exporter needs from a CellProfiler pipeline object (no measurements needed)."""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .names import Feature, split_pair

LOGGER = logging.getLogger(__name__)
IMAGE, EXPERIMENT = "Image", "Experiment"
OBJECTS_FILE = "ObjectsFileName_"   # Image column NamesAndTypes/LoadData writes per file-loaded object
LOAD_MODULES = ("NamesAndTypes", "LoadData", "LoadImages")
S_PRIMARY_OUT = "Name the primary objects to be identified"
S_SECONDARY_IN, S_SECONDARY_OUT = "Select the input objects", "Name the objects to be identified"
S_TERTIARY_BIG, S_TERTIARY_SMALL = "Select the larger identified objects", "Select the smaller identified objects"
S_TERTIARY_OUT = "Name the tertiary objects to be identified"
S_FILTER_IN, S_FILTER_OUT = "Select the objects to filter", "Name the output objects"
S_FILTER_EXTRA_IN, S_FILTER_EXTRA_OUT = "Select additional object to relabel", "Name the relabeled objects"


class RoleError(ValueError):
    pass


@dataclass
class ObjectInfo:
    name: str
    module_num: int
    module_name: str
    role: Optional[str] = None


@dataclass
class Context:
    channels: List[str]
    objects: Dict[str, ObjectInfo]
    roles: Dict[str, str]
    features: List[Feature]
    modules: List[dict] = field(default_factory=list)
    metadata_tags: List[str] = field(default_factory=list)
    role_note: Dict[str, object] = field(default_factory=dict)


def _settings_dict(module) -> Dict[str, str]:
    """{setting text: value}. Repeated texts (FilterObjects extra objects) keep the first occurrence;
    the extra-object pairs are read positionally in _filter_pairs. Settings with no text (hidden/internal
    settings some modules serialize as a bare ':' line, e.g. Images) are skipped: an empty-string key can't
    round-trip as an HDF5 dataset/group name when this ends up in provenance()."""
    out = {}
    for s in module.settings():
        text = getattr(s, "text", None)
        if text and text not in out:
            out[text] = getattr(s, "value_text", getattr(s, "value", None))
    return out


def _setting_values(module) -> List[tuple]:
    """(text, value) for every setting in order, repeats kept: MeasureTexture writes one
    "Texture scale to measure" per scale and _settings_dict collapses them."""
    out = []
    for s in module.settings():
        text = getattr(s, "text", None)
        if text:
            out.append((text, getattr(s, "value_text", getattr(s, "value", None))))
    return out


def _filter_pairs(settings_list) -> List[tuple]:
    """(input, output) pairs for FilterObjects: the main pair plus every additional relabel pair."""
    pairs, extra_in = [], None
    d = {}
    for s in settings_list:
        d.setdefault(s.text, s.value_text if hasattr(s, "value_text") else s.value)
    if S_FILTER_IN in d and S_FILTER_OUT in d:
        pairs.append((d[S_FILTER_IN], d[S_FILTER_OUT]))
    for s in settings_list:
        v = s.value_text if hasattr(s, "value_text") else s.value
        if s.text == S_FILTER_EXTRA_IN:
            extra_in = v
        elif s.text == S_FILTER_EXTRA_OUT and extra_in is not None:
            pairs.append((extra_in, v))
            extra_in = None
    return pairs


def _referencing_modules(name: str, objects: Dict[str, ObjectInfo], modules: List[dict]) -> int:
    """How many modules *other than the one that produced it* name this object in one of their
    settings, i.e. use it as an input: RelateObjects parent/child, FilterObjects source,
    MaskObjects, the MeasureObject* "Select objects to measure" lists (those hold ", "-joined
    names, so every value is split on ", " before comparing)."""
    producer = objects[name].module_num if name in objects else None
    hits = 0
    for m in modules:
        if m["num"] == producer:
            continue
        # build_context always sets both: setting_values keeps repeated setting texts, settings is
        # the collapsed first-occurrence dict a hand-built module entry may carry instead
        values = m.get("setting_values") or list(m.get("settings", {}).items())
        if any(name == part for _, v in values if v is not None for part in str(v).split(", ")):
            hits += 1
    return hits


def _relabel_pairs(modules: List[dict]) -> List[tuple]:
    """The FilterObjects (input, output) pairs that are a *rename*: only a module that relabels
    several objects together (the JUMP border-cleanup step, which relabels nucleus and cell in one
    go) renames what it filters. A single-pair FilterObjects makes a plain filtered subset: its
    output inherits its input's role and competes as an ordinary candidate, but it never stands in
    for its input -- substituting it would silently export the subset (e.g. PH3PosNuclei) instead of
    the object the pipeline is about (Nuclei)."""
    return [pair for m in modules
            if m["name"] == "FilterObjects" and len(m.get("filter_pairs", [])) > 1
            for pair in m["filter_pairs"]]


def _ordered(names, objects: Dict[str, ObjectInfo]) -> List[str]:
    """Candidate objects in pipeline order: producing module first, then name."""
    return sorted(set(names), key=lambda n: (objects[n].module_num if n in objects else 0, n))


def _most_related(candidates: List[str], objects: Dict[str, ObjectInfo], modules: List[dict]) -> str:
    """The candidate the rest of the pipeline uses most as an input; ties go to the object whose
    producing module comes first."""
    return min(candidates, key=lambda n: (-_referencing_modules(n, objects, modules),
                                          objects[n].module_num if n in objects else 0, n))


def _pick_chain_module(mods: List[dict], out_setting: str, fallback: str,
                       objects: Dict[str, ObjectInfo], modules: List[dict]):
    """(module, note) for the Identify{Secondary,Tertiary}Objects module whose chain becomes the
    roles. Several of them are no longer a refusal and no longer cost the pipeline its other
    compartments: the module whose output object the rest of the pipeline uses most as an input
    wins (ties -> the one produced first), and the whole chain is read off that one module."""
    named = [m for m in mods if m["settings"].get(out_setting)]
    if not named:
        raise RoleError("No %s module names its output objects (%s); set the object roles explicitly."
                        % (mods[0]["name"], out_setting))
    by_out = {m["settings"][out_setting]: m for m in named}
    if len(by_out) < 2:
        return named[0], {"fallback": None, "candidates": []}
    candidates = _ordered(by_out, objects)
    chosen = _most_related(candidates, objects, modules)
    LOGGER.warning("ExportToAnnData: several %s modules (%s); picked %s, the one the most other modules "
                   "use as an input, and took its whole object chain. Set the object roles explicitly "
                   "to override.", mods[0]["name"], ", ".join(candidates), chosen)
    return by_out[chosen], {"fallback": fallback, "candidates": candidates}


def _chain_roles(mod: dict, settings_by_role: Dict[str, str]) -> Dict[str, str]:
    """{role: object} read off one Identify module's settings; a module that leaves one of them
    empty is malformed and yields a RoleError rather than a KeyError."""
    roles = {role: mod["settings"].get(setting) for role, setting in settings_by_role.items()}
    blank = [setting for role, setting in settings_by_role.items() if not roles[role]]
    if blank:
        raise RoleError("%s (module %d) does not name: %s. Set the object roles explicitly."
                        % (mod["name"], mod["num"], ", ".join(blank)))
    return roles


def detect_roles(objects: Dict[str, ObjectInfo], modules: List[dict]):
    """(roles, note): primary / secondary / tertiary objects from module types and their
    input/output settings. In a standard pipeline primary = nucleus, secondary = cell,
    tertiary = cytoplasm.

    Priority: IdentifyTertiaryObjects names all three; else IdentifySecondaryObjects names
    primary + secondary; else the primary candidates. Where any of those is ambiguous the
    most-related candidate is picked automatically (never a refusal) and `note` records the
    choice for uns["cellprofiler"].
    """
    # follow FilterObjects renames: output -> ultimate origin role
    origin = {}
    for m in modules:
        s = m["settings"]
        if m["name"] == "IdentifyPrimaryObjects" and s.get(S_PRIMARY_OUT):
            origin[s[S_PRIMARY_OUT]] = "primary"
        elif m["name"] == "IdentifySecondaryObjects" and s.get(S_SECONDARY_OUT):
            origin[s[S_SECONDARY_OUT]] = "secondary"
        elif m["name"] == "IdentifyTertiaryObjects" and s.get(S_TERTIARY_OUT):
            origin[s[S_TERTIARY_OUT]] = "tertiary"
        elif m["name"] == "FilterObjects":
            # a filtered object takes its input's role: it is a candidate in its own right, whether
            # or not this module is a rename (_relabel_pairs decides that, further down)
            for src, dst in m.get("filter_pairs", []):
                if src in origin:
                    origin[dst] = origin[src]
    note = {"fallback": None, "candidates": []}
    tertiaries = [m for m in modules if m["name"] == "IdentifyTertiaryObjects"]
    secondaries = [m for m in modules if m["name"] == "IdentifySecondaryObjects"]
    if tertiaries:
        mod, note = _pick_chain_module(tertiaries, S_TERTIARY_OUT, "most_related_tertiary",
                                       objects, modules)
        roles = _chain_roles(mod, {"primary": S_TERTIARY_SMALL, "secondary": S_TERTIARY_BIG,
                                   "tertiary": S_TERTIARY_OUT})
    elif secondaries:
        mod, note = _pick_chain_module(secondaries, S_SECONDARY_OUT, "most_related_secondary",
                                       objects, modules)
        roles = _chain_roles(mod, {"primary": S_SECONDARY_IN, "secondary": S_SECONDARY_OUT})
    else:
        pool = [n for n, r in origin.items() if r == "primary" and n in objects] or list(objects)
        # a relabel-together FilterObjects renamed some of these: only the name that survives is a
        # real candidate, so a single renamed chain scores as one object and needs no fallback
        renamed = {src for src, dst in _relabel_pairs(modules) if dst in objects}
        candidates = _ordered([n for n in pool if n not in renamed] or pool, objects)
        if not candidates:
            raise RoleError("This pipeline produces no objects, so there is nothing to export. Add an "
                            "Identify*Objects module (or load objects) before ExportToAnnData.")
        chosen = candidates[0]
        if len(candidates) > 1:
            chosen = _most_related(candidates, objects, modules)
            note = {"fallback": "most_related_primary", "candidates": list(candidates)}
            LOGGER.warning("ExportToAnnData: several candidate primary objects (%s); picked %s, the one "
                           "the most other modules use as an input. Set the object roles explicitly to "
                           "override.", ", ".join(candidates), chosen)
        roles = {"primary": chosen}
    # FilterObjects may have renamed the objects the tertiary/secondary consumed -> keep the consumed names,
    # but if a *later* FilterObjects renamed an object we chose, prefer the renamed (filtered) one.
    for role, name in list(roles.items()):
        for src, dst in _relabel_pairs(modules):
            if src == name and dst in objects and role != "tertiary":
                roles[role] = dst
    missing = [n for n in roles.values() if n not in objects]
    if missing:
        raise RoleError("Role objects not produced by the pipeline: " + ", ".join(missing))
    return roles, note


def _channels(pipeline) -> List[str]:
    prov = pipeline.get_provider_dictionary("imagegroup")
    loaded = [n for n, lst in prov.items() if any(m.module_name in LOAD_MODULES for m, _ in lst)]
    derived = [n for n in prov if n not in loaded]
    return loaded + derived


def _objects(pipeline) -> Dict[str, ObjectInfo]:
    out = {}
    for name, lst in pipeline.get_provider_dictionary("objectgroup").items():
        module, _ = lst[-1]
        out[name] = ObjectInfo(name, module.module_num, module.module_name)
    return out


def _features(pipeline, objects, channels) -> List[Feature]:
    columns = {}
    for c in pipeline.get_measurement_columns():
        columns.setdefault((c[0], c[1]), c[2])
    feats: Dict[tuple, Feature] = {}
    targets = list(objects) + [IMAGE, EXPERIMENT]
    for module in pipeline.modules():
        # a third-party module with a sloppy get_categories/get_measurements must not break the
        # pipeline load: its columns simply take the name-splitting fallback below
        try:
            for obj in targets:
                for cat in module.get_categories(pipeline, obj):
                    for meas in module.get_measurements(pipeline, obj, cat):
                        images = list(module.get_measurement_images(pipeline, obj, cat, meas)) or [None]
                        others = list(module.get_measurement_objects(pipeline, obj, cat, meas)) or [None]
                        for img in images:
                            scales = list(module.get_measurement_scales(pipeline, obj, cat, meas, img)) or [None]
                            for other in others:
                                for sc in scales:
                                    name = "_".join(p for p in (cat, meas, other, img, sc) if p)
                                    if (obj, name) not in columns or (obj, name) in feats:
                                        continue
                                    img1, img2 = img, None
                                    if cat == "Correlation" and img:
                                        pair = split_pair(img, channels)
                                        if pair:
                                            img1, img2 = pair
                                    feats[(obj, name)] = Feature(
                                        object=obj, cp_name=name, category=cat, measurement=meas,
                                        module_num=module.module_num, module_name=module.module_name,
                                        image=img1, image2=img2, other_object=other, scale=sc,
                                        coltype=str(columns[(obj, name)]), parsed_by="api")
        except Exception as e:
            LOGGER.warning("ExportToAnnData: module %s: introspection failed (%s); its columns fall back "
                           "to name splitting.", getattr(module, "module_name", module), e)
    for (obj, name), coltype in columns.items():
        if (obj, name) not in feats:
            cat, _, meas = name.partition("_")
            feats[(obj, name)] = Feature(object=obj, cp_name=name, category=cat, measurement=meas,
                                         coltype=str(coltype), parsed_by="fallback")
    return list(feats.values())


def build_context(pipeline, roles: Optional[Dict[str, str]] = None) -> Context:
    channels = _channels(pipeline)
    objects = _objects(pipeline)
    modules = []
    for m in pipeline.modules():
        entry = {"num": m.module_num, "name": m.module_name, "enabled": getattr(m, "enabled", True),
                 "settings": _settings_dict(m), "setting_values": _setting_values(m)}
        if m.module_name == "FilterObjects":
            entry["filter_pairs"] = _filter_pairs(m.settings())
        modules.append(entry)
    if roles:
        unknown = [v for v in roles.values() if v not in objects]
        if unknown:
            raise RoleError("Unknown object(s) for explicit roles: " + ", ".join(unknown))
        # exactly the roles the user set: predictable, and it cannot pair an explicit object with
        # an auto-detected one the pipeline never relates to it
        resolved = dict(roles)
        note = {"mode": "manual", "fallback": None, "candidates": []}
    else:
        resolved, detected = detect_roles(objects, modules)
        note = {"mode": "automatic", "fallback": detected["fallback"],
                "candidates": detected["candidates"]}
    for role, name in resolved.items():
        objects[name].role = role
    features = _features(pipeline, objects, channels)
    tags = sorted({f.cp_name[len("Metadata_"):] for f in features
                   if f.object == IMAGE and f.cp_name.startswith("Metadata_")})
    return Context(channels=channels, objects=objects, roles=resolved, features=features,
                   modules=modules, metadata_tags=tags, role_note=note)


def file_loaded_objects(ctx: Context) -> List[str]:
    """Objects read from label files (NamesAndTypes/LoadData "Objects") rather than made by an
    Identify* module: they carry no parent/child links unless RelateObjects makes them."""
    return sorted({f.cp_name[len(OBJECTS_FILE):] for f in ctx.features
                   if f.object == IMAGE and f.cp_name.startswith(OBJECTS_FILE)})
