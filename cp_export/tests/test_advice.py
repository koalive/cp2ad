"""Spec section 5 advice messages."""
import dataclasses

from cpexport.advice import advice
from cpexport.assemble import provenance
from cpexport.introspect import build_context
from cpexport.names import Feature
from test_assemble import make_measurements

OBJECTS_FROM_FILE = Feature(object="Image", cp_name="ObjectsFileName_Cells", category="ObjectsFileName",
                            measurement="Cells", coltype="varchar(256)")
SPREADSHEET_MODULE = {"num": 99, "name": "ExportToSpreadsheet", "enabled": True, "settings": {},
                      "setting_values": []}


def test_metadata_advice_names_the_missing_tags(fake_pipeline):
    messages = advice(build_context(fake_pipeline))
    assert len(messages) == 1 and "Metadata_Plate/Well/Site" in messages[0]


def test_role_fallback_advice_names_the_choice_and_the_override(fake_pipeline):
    ctx = build_context(fake_pipeline)
    ctx = dataclasses.replace(ctx, roles={"primary": "PH3"},
                              role_note={"mode": "automatic", "fallback": "most_related_primary",
                                         "candidates": ["Nuclei", "PH3"]})
    first = advice(ctx)[0]
    assert "Nuclei, PH3" in first and "PH3 was picked" in first
    assert "How to pick primary / secondary / tertiary objects" in first and "Manual" in first


def test_chain_fallback_advice_names_the_module_type_and_the_chain(fake_pipeline):
    ctx = build_context(fake_pipeline)
    ctx = dataclasses.replace(ctx, roles={"primary": "Nuclei", "secondary": "Cells", "tertiary": "Cytoplasm"},
                              role_note={"mode": "automatic", "fallback": "most_related_tertiary",
                                         "candidates": ["Cytoplasm", "Rings"]})
    first = advice(ctx)[0]
    assert "IdentifyTertiaryObjects" in first and "Cytoplasm, Rings" in first
    assert "Cytoplasm was picked" in first and "Manual" in first


def test_file_loaded_objects_advice(fake_pipeline):
    ctx = build_context(fake_pipeline)
    ctx = dataclasses.replace(ctx, features=list(ctx.features) + [OBJECTS_FROM_FILE])
    text = "\n".join(advice(ctx))
    assert "Objects Cells are loaded from files" in text and "RelateObjects" in text


def test_file_loaded_objects_are_marked_in_provenance(fake_pipeline, meas_arrays):
    ctx = build_context(fake_pipeline)
    ctx = dataclasses.replace(ctx, features=list(ctx.features) + [OBJECTS_FROM_FILE])
    objects = provenance(ctx, make_measurements(meas_arrays), {})["objects"]
    assert objects["Cells"]["source"] == "file" and objects["Nuclei"]["source"] == "pipeline"


def test_exporttospreadsheet_advice(fake_pipeline):
    ctx = build_context(fake_pipeline)
    ctx = dataclasses.replace(ctx, modules=list(ctx.modules) + [SPREADSHEET_MODULE])
    text = "\n".join(advice(ctx))
    for expected in ("Add image metadata columns: Yes", "Add image file and folder names: Yes",
                     "Export all measurement types: Yes", "aggregation", "Representation of NaN: NaN",
                     "one file per object"):
        assert expected in text, expected
