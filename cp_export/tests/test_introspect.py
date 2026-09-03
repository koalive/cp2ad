import logging

import pytest
from cpexport.introspect import Context, ObjectInfo, RoleError, _filter_pairs, build_context, detect_roles
from cpexport.names import to_cpm_names
from conftest import FakeModule, FakePipeline, FakeSetting


def test_channels_in_namesandtypes_order(fake_pipeline):
    ctx = build_context(fake_pipeline)
    assert ctx.channels == ["DNA", "PH3", "cellbody"]


def test_objects_and_roles_examplehuman(fake_pipeline):
    ctx = build_context(fake_pipeline)
    assert set(ctx.objects) == {"Nuclei", "PH3", "Cells", "Cytoplasm"}
    assert ctx.objects["Cells"].module_name == "IdentifySecondaryObjects"
    assert ctx.roles == {"primary": "Nuclei", "secondary": "Cells", "tertiary": "Cytoplasm"}
    assert ctx.role_note == {"mode": "automatic", "fallback": None, "candidates": []}


def test_explicit_roles_override(fake_pipeline):
    ctx = build_context(fake_pipeline, roles={"primary": "PH3", "secondary": "Cells", "tertiary": "Cytoplasm"})
    assert ctx.roles["primary"] == "PH3"
    assert ctx.role_note == {"mode": "manual", "fallback": None, "candidates": []}


def test_manual_roles_are_exact(fake_pipeline):
    """Manual mode exports exactly the roles that were set -- no auto-filled companions, which could
    pair the chosen object with one the pipeline never relates to it."""
    ctx = build_context(fake_pipeline, roles={"primary": "PH3"})
    assert ctx.roles == {"primary": "PH3"}
    assert ctx.role_note["mode"] == "manual"


def test_explicit_role_unknown_object_raises(fake_pipeline):
    with pytest.raises(RoleError):
        build_context(fake_pipeline, roles={"primary": "Nope"})


def test_every_declared_column_becomes_a_feature(fake_pipeline, probe):
    ctx = build_context(fake_pipeline)
    declared = {(c[0], c[1]) for c in probe["measurement_columns"]}
    got = {(f.object, f.cp_name) for f in ctx.features}
    assert declared == got


def test_api_parse_coverage(fake_pipeline):
    ctx = build_context(fake_pipeline)
    per_obj = [f for f in ctx.features if f.object not in ("Image", "Experiment")]
    assert per_obj and all(f.parsed_by == "api" for f in per_obj)
    fallback = {f.cp_name.split("_")[0] for f in ctx.features if f.parsed_by == "fallback"}
    assert fallback <= {"ExecutionTime", "ModuleError", "Group", "Granularity", "Pipeline", "CellProfiler",
                        "Run", "Modification"}


def test_structured_fields(fake_pipeline):
    ctx = build_context(fake_pipeline)
    by = {(f.object, f.cp_name): f for f in ctx.features}
    t = by[("Cells", "Texture_Contrast_DNA_3_00_256")]
    assert (t.category, t.measurement, t.image, t.scale, t.module_name) == \
        ("Texture", "Contrast", "DNA", "3_00_256", "MeasureTexture")
    c = by[("Cells", "Correlation_Manders_PH3_DNA")]
    assert (c.image, c.image2) == ("PH3", "DNA")
    n = by[("Cells", "Neighbors_NumberOfNeighbors_5")]
    assert (n.scale, n.other_object) == ("5", None)
    assert by[("Cells", "Parent_Nuclei")].coltype == "integer"


def test_cpm_names_are_squidpy_columns(fake_pipeline, cpm_columns):
    ctx = build_context(fake_pipeline)
    cols = set(cpm_columns)
    for f in ctx.features:
        if f.object != "Cells":
            continue
        for name, backend in to_cpm_names(f, ctx.channels):
            if backend == "cp_measure":
                assert name in cols, (f.cp_name, name)


def test_metadata_tags(fake_pipeline):
    # ExampleHuman fixture has metadata extraction off; Metadata_Frame/Series are runtime-only
    assert build_context(fake_pipeline).metadata_tags == []


def _mods(*names):
    mods = []
    for i, spec in enumerate(names):
        name, settings = spec[0], spec[1]
        entry = {"num": i + 1, "name": name, "enabled": True, "settings": settings,
                 "setting_values": list(settings.items())}
        if len(spec) > 2:
            entry["filter_pairs"] = spec[2]
        mods.append(entry)
    return mods


def test_detect_roles_jump_style_filterobjects_chain():
    objects = {
        "NucleiIncludingEdges": ObjectInfo("NucleiIncludingEdges", 1, "IdentifyPrimaryObjects"),
        "CellsIncludingEdges": ObjectInfo("CellsIncludingEdges", 2, "IdentifySecondaryObjects"),
        "Nuclei": ObjectInfo("Nuclei", 3, "FilterObjects"),
        "Cells": ObjectInfo("Cells", 3, "FilterObjects"),
        "Cytoplasm": ObjectInfo("Cytoplasm", 4, "IdentifyTertiaryObjects"),
    }
    modules = _mods(
        ("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "NucleiIncludingEdges"}),
        ("IdentifySecondaryObjects", {"Select the input objects": "NucleiIncludingEdges",
                                      "Name the objects to be identified": "CellsIncludingEdges"}),
        ("FilterObjects", {"Select the objects to filter": "NucleiIncludingEdges", "Name the output objects": "Nuclei",
                           "Select additional object to relabel": "CellsIncludingEdges",
                           "Name the relabeled objects": "Cells"},
         [("NucleiIncludingEdges", "Nuclei"), ("CellsIncludingEdges", "Cells")]),
        ("IdentifyTertiaryObjects", {"Select the larger identified objects": "Cells",
                                     "Select the smaller identified objects": "Nuclei",
                                     "Name the tertiary objects to be identified": "Cytoplasm"}),
    )
    assert detect_roles(objects, modules)[0] == {"primary": "Nuclei", "secondary": "Cells",
                                                 "tertiary": "Cytoplasm"}


def test_detect_roles_secondary_then_filterobjects_rename():
    """No tertiary module: primary/secondary come from IdentifySecondaryObjects, then a later
    FilterObjects renames both -- exercises the origin-tracking + rename-preference pass in
    detect_roles (previously reachable only alongside a tertiary, whose branch alone decided
    the answer)."""
    objects = {
        "NucleiIncludingEdges": ObjectInfo("NucleiIncludingEdges", 1, "IdentifyPrimaryObjects"),
        "CellsIncludingEdges": ObjectInfo("CellsIncludingEdges", 2, "IdentifySecondaryObjects"),
        "Nuclei": ObjectInfo("Nuclei", 3, "FilterObjects"),
        "Cells": ObjectInfo("Cells", 3, "FilterObjects"),
    }
    modules = _mods(
        ("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "NucleiIncludingEdges"}),
        ("IdentifySecondaryObjects", {"Select the input objects": "NucleiIncludingEdges",
                                      "Name the objects to be identified": "CellsIncludingEdges"}),
        ("FilterObjects", {"Select the objects to filter": "NucleiIncludingEdges", "Name the output objects": "Nuclei",
                           "Select additional object to relabel": "CellsIncludingEdges",
                           "Name the relabeled objects": "Cells"},
         [("NucleiIncludingEdges", "Nuclei"), ("CellsIncludingEdges", "Cells")]),
    )
    assert detect_roles(objects, modules)[0] == {"primary": "Nuclei", "secondary": "Cells"}


def _chain_objects(*specs):
    return {name: ObjectInfo(name, num, mod) for name, num, mod in specs}


def test_detect_roles_several_tertiaries_picks_the_most_related_chain(caplog):
    """Two complete Identify chains: the tertiary object the rest of the pipeline measures wins and
    its own smaller/larger settings become primary/secondary -- the other chain is not exported, but
    no compartment of the chosen chain is lost."""
    objects = _chain_objects(("A1", 1, "IdentifyPrimaryObjects"), ("B1", 2, "IdentifySecondaryObjects"),
                             ("C1", 3, "IdentifyTertiaryObjects"), ("A2", 4, "IdentifyPrimaryObjects"),
                             ("B2", 5, "IdentifySecondaryObjects"), ("C2", 6, "IdentifyTertiaryObjects"))
    modules = _mods(
        ("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "A1"}),
        ("IdentifySecondaryObjects", {"Select the input objects": "A1",
                                      "Name the objects to be identified": "B1"}),
        ("IdentifyTertiaryObjects", {"Select the larger identified objects": "B1",
                                     "Select the smaller identified objects": "A1",
                                     "Name the tertiary objects to be identified": "C1"}),
        ("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "A2"}),
        ("IdentifySecondaryObjects", {"Select the input objects": "A2",
                                      "Name the objects to be identified": "B2"}),
        ("IdentifyTertiaryObjects", {"Select the larger identified objects": "B2",
                                     "Select the smaller identified objects": "A2",
                                     "Name the tertiary objects to be identified": "C2"}),
        ("MeasureObjectIntensity", {"Select objects to measure": "C2"}),
    )
    with caplog.at_level(logging.WARNING, logger="cpexport.introspect"):
        roles, note = detect_roles(objects, modules)
    assert roles == {"primary": "A2", "secondary": "B2", "tertiary": "C2"}
    assert note == {"fallback": "most_related_tertiary", "candidates": ["C1", "C2"]}
    assert any("IdentifyTertiaryObjects" in r.getMessage() for r in caplog.records)


def test_detect_roles_several_secondaries_picks_the_most_related_chain():
    """No tertiary module at all: the most-related secondary output wins and its input is the primary."""
    objects = _chain_objects(("A1", 1, "IdentifyPrimaryObjects"), ("B1", 2, "IdentifySecondaryObjects"),
                             ("A2", 3, "IdentifyPrimaryObjects"), ("B2", 4, "IdentifySecondaryObjects"))
    modules = _mods(
        ("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "A1"}),
        ("IdentifySecondaryObjects", {"Select the input objects": "A1",
                                      "Name the objects to be identified": "B1"}),
        ("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "A2"}),
        ("IdentifySecondaryObjects", {"Select the input objects": "A2",
                                      "Name the objects to be identified": "B2"}),
        ("MeasureObjectIntensity", {"Select objects to measure": "B2, A2"}),
    )
    roles, note = detect_roles(objects, modules)
    assert roles == {"primary": "A2", "secondary": "B2"}
    assert note == {"fallback": "most_related_secondary", "candidates": ["B1", "B2"]}


def test_detect_roles_primary_only():
    objects = {"Nuclei": ObjectInfo("Nuclei", 1, "IdentifyPrimaryObjects")}
    modules = _mods(("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "Nuclei"}))
    roles, note = detect_roles(objects, modules)
    assert roles == {"primary": "Nuclei"} and note == {"fallback": None, "candidates": []}


def _two_primaries():
    return {"A": ObjectInfo("A", 1, "IdentifyPrimaryObjects"), "B": ObjectInfo("B", 2, "IdentifyPrimaryObjects")}


def test_detect_roles_ambiguous_primary_picks_the_most_related(caplog):
    """Two primaries, no secondary/tertiary chain: the one the rest of the pipeline uses as an input
    wins (B is referenced by RelateObjects and MeasureObjectIntensity, A by nothing)."""
    modules = _mods(
        ("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "A"}),
        ("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "B"}),
        ("RelateObjects", {"Parent objects": "B", "Child objects": "Speckles"}),
        ("MeasureObjectIntensity", {"Select objects to measure": "B, Speckles"}),
    )
    with caplog.at_level(logging.WARNING, logger="cpexport.introspect"):
        roles, note = detect_roles(_two_primaries(), modules)
    assert roles == {"primary": "B"}
    assert note == {"fallback": "most_related_primary", "candidates": ["A", "B"]}
    assert any("A, B" in r.getMessage() for r in caplog.records)


def test_detect_roles_ambiguous_primary_tie_goes_to_the_first_module():
    """Neither primary is used as an input anywhere: the earlier producing module wins."""
    modules = _mods(("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "A"}),
                    ("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "B"}))
    roles, note = detect_roles(_two_primaries(), modules)
    assert roles == {"primary": "A"} and note["candidates"] == ["A", "B"]


def test_detect_roles_single_pair_filterobjects_is_not_a_rename():
    """ExamplePercentPositive shape: FilterObjects makes PH3PosNuclei, a *subset* of Nuclei. The
    subset is scored like any other candidate, but it never stands in for Nuclei -- substituting it
    would silently export the 20 positive nuclei instead of all of them."""
    objects = _chain_objects(("Nuclei", 1, "IdentifyPrimaryObjects"), ("PH3", 2, "IdentifyPrimaryObjects"),
                             ("PH3PosNuclei", 4, "FilterObjects"))
    modules = _mods(
        ("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "Nuclei"}),
        ("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "PH3"}),
        ("RelateObjects", {"Parent objects": "Nuclei", "Child objects": "PH3"}),
        ("FilterObjects", {"Select the objects to filter": "Nuclei", "Name the output objects": "PH3PosNuclei"},
         [("Nuclei", "PH3PosNuclei")]),
        ("MeasureObjectIntensity", {"Select objects to measure": "Nuclei, PH3PosNuclei"}),
    )
    roles, note = detect_roles(objects, modules)
    assert roles == {"primary": "Nuclei"}        # RelateObjects + FilterObjects + MeasureObjectIntensity
    assert note == {"fallback": "most_related_primary",
                    "candidates": ["Nuclei", "PH3", "PH3PosNuclei"]}


def test_detect_roles_filtered_object_wins_when_the_pipeline_measures_it():
    """The other half of the same rule: a filtered subset is a candidate in its own right, so a
    pipeline that only measures FilteredNuclei exports FilteredNuclei, not Nuclei."""
    objects = _chain_objects(("Nuclei", 1, "IdentifyPrimaryObjects"), ("FilteredNuclei", 2, "FilterObjects"))
    modules = _mods(
        ("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "Nuclei"}),
        ("FilterObjects", {"Select the objects to filter": "Nuclei", "Name the output objects": "FilteredNuclei"},
         [("Nuclei", "FilteredNuclei")]),
        ("MeasureObjectIntensity", {"Select objects to measure": "FilteredNuclei"}),
        ("MeasureObjectSizeShape", {"Select objects to measure": "FilteredNuclei"}),
    )
    roles, note = detect_roles(objects, modules)
    assert roles == {"primary": "FilteredNuclei"}
    assert note == {"fallback": "most_related_primary", "candidates": ["Nuclei", "FilteredNuclei"]}


def test_detect_roles_relabelled_chain_is_one_candidate(caplog):
    """A relabel-together FilterObjects renames NucleiIncludingEdges -> Nuclei. Both names must not
    count as two candidates: that is one chain, so there is no ambiguity and no fallback."""
    objects = _chain_objects(("NucleiIncludingEdges", 1, "IdentifyPrimaryObjects"),
                             ("MaskedSpots", 2, "MaskObjects"), ("Nuclei", 3, "FilterObjects"),
                             ("Spots", 3, "FilterObjects"))
    modules = _mods(
        ("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "NucleiIncludingEdges"}),
        ("MaskObjects", {"Select objects to be masked": "MaskedSpots"}),
        ("FilterObjects", {"Select the objects to filter": "NucleiIncludingEdges",
                           "Name the output objects": "Nuclei",
                           "Select additional object to relabel": "MaskedSpots",
                           "Name the relabeled objects": "Spots"},
         [("NucleiIncludingEdges", "Nuclei"), ("MaskedSpots", "Spots")]),
    )
    with caplog.at_level(logging.WARNING, logger="cpexport.introspect"):
        roles, note = detect_roles(objects, modules)
    assert roles == {"primary": "Nuclei"}
    assert note == {"fallback": None, "candidates": []}
    assert not caplog.records


def _tertiary_pipeline(tertiary_settings):
    objects = _chain_objects(("A", 1, "IdentifyPrimaryObjects"), ("B", 2, "IdentifySecondaryObjects"))
    return objects, _mods(
        ("IdentifyPrimaryObjects", {"Name the primary objects to be identified": "A"}),
        ("IdentifySecondaryObjects", {"Select the input objects": "A",
                                      "Name the objects to be identified": "B"}),
        ("IdentifyTertiaryObjects", tertiary_settings),
    )


def test_detect_roles_malformed_chain_module_raises_roleerror_not_keyerror():
    """A module that leaves one of the settings the roles are read from empty must not crash the
    pipeline load with a KeyError."""
    no_output = _tertiary_pipeline({"Select the larger identified objects": "B",
                                    "Select the smaller identified objects": "A"})
    with pytest.raises(RoleError, match="output objects"):
        detect_roles(*no_output)
    no_smaller = _tertiary_pipeline({"Select the larger identified objects": "B",
                                     "Name the tertiary objects to be identified": "C"})
    with pytest.raises(RoleError, match="does not name"):
        detect_roles(*no_smaller)


def test_detect_roles_without_objects_raises():
    with pytest.raises(RoleError, match="no objects"):
        detect_roles({}, _mods(("Images", {})))


def test_filter_pairs_reads_main_and_additional_pairs():
    settings = [
        FakeSetting("Select the objects to filter", "NucleiIncludingEdges"),
        FakeSetting("Name the output objects", "Nuclei"),
        FakeSetting("Select additional object to relabel", "CellsIncludingEdges"),
        FakeSetting("Name the relabeled objects", "Cells"),
        FakeSetting("Select additional object to relabel", "Cyto0"),
        FakeSetting("Name the relabeled objects", "Cyto"),
    ]
    assert _filter_pairs(settings) == [
        ("NucleiIncludingEdges", "Nuclei"), ("CellsIncludingEdges", "Cells"), ("Cyto0", "Cyto"),
    ]


class _FilterObjectsModule:
    """Minimal fake module exposing the ordered, repeated-text settings() list _filter_pairs needs;
    conftest.FakeModule collapses repeated setting texts into a dict, which would hide the
    additional-object pairs."""

    module_num, module_name, enabled = 1, "FilterObjects", True

    def __init__(self, settings_list):
        self._settings = settings_list

    def settings(self):
        return self._settings

    def get_categories(self, pipeline, object_name):
        return []


def test_build_context_populates_filter_pairs():
    settings = [
        FakeSetting("Select the objects to filter", "NucleiIncludingEdges"),
        FakeSetting("Name the output objects", "Nuclei"),
        FakeSetting("Select additional object to relabel", "CellsIncludingEdges"),
        FakeSetting("Name the relabeled objects", "Cells"),
        FakeSetting("Select additional object to relabel", "Cyto0"),
        FakeSetting("Name the relabeled objects", "Cyto"),
    ]
    module = _FilterObjectsModule(settings)
    pipeline = FakePipeline(
        [module], columns=[], providers_images={},
        providers_objects={"Nuclei": [(None, 1, "Name the output objects")]},
    )
    # Explicit, non-empty roles bypass detect_roles (which this minimal pipeline can't satisfy);
    # this test only cares that build_context recorded the module's filter_pairs.
    ctx = build_context(pipeline, roles={"primary": "Nuclei"})
    assert ctx.modules[0]["filter_pairs"] == [
        ("NucleiIncludingEdges", "Nuclei"), ("CellsIncludingEdges", "Cells"), ("Cyto0", "Cyto"),
    ]


class _BrokenModule(FakeModule):
    """A third-party module whose introspection API raises, as a sloppy plugin's would."""

    def get_categories(self, pipeline, object_name):
        raise RuntimeError("third-party plugin bug")


def _pipeline_with_broken_module(probe, broken_num):
    intro = probe["introspection"]
    modules = []
    for mod in probe["modules"]:
        cls = _BrokenModule if mod["num"] == broken_num else FakeModule
        modules.append(cls(mod["num"], mod["name"], mod["settings"],
                           intro.get(f"{mod['num']}:{mod['name']}", []), mod["enabled"]))
    return FakePipeline(modules, probe["measurement_columns"], probe["providers_images"],
                        probe["providers_objects"])


def test_broken_module_introspection_falls_back(probe, caplog):
    with caplog.at_level(logging.WARNING, logger="cpexport.introspect"):
        ctx = build_context(_pipeline_with_broken_module(probe, 10))  # MeasureObjectIntensity
    intensity = [f for f in ctx.features if f.object == "Cells" and f.category == "Intensity"]
    assert intensity and all(f.parsed_by == "fallback" for f in intensity)
    assert any("introspection failed" in r.getMessage() for r in caplog.records)


def test_setting_values_keep_repeated_texts(fake_pipeline):
    ctx = build_context(fake_pipeline)
    texture = next(m for m in ctx.modules if m["name"] == "MeasureTexture")
    assert [v for t, v in texture["setting_values"] if t == "Texture scale to measure"] == ["3", "5", "10"]
