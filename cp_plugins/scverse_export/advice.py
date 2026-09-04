"""Spec section 5 advice: what to change in the pipeline so the exported object carries maximum
information. Plain strings, so the module can raise them as one ValidationError at pipeline-load
time and log the same text during the run."""
from __future__ import annotations

from typing import List

from .introspect import Context, file_loaded_objects

ROLE_FALLBACK = (
    "Several objects could be the primary object (%s) and the pipeline has no "
    "IdentifySecondaryObjects/IdentifyTertiaryObjects chain to resolve them, so %s was picked "
    "automatically: it is the one the most other modules use as an input. Only its measurements are "
    "exported. To pick another, set \"How to pick primary / secondary / tertiary objects\" to Manual "
    "and name the objects there.")
CHAIN_FALLBACK = (
    "The pipeline has several %s modules and their output objects (%s) describe different object "
    "chains; %s was picked automatically, because it is the one the most other modules use as an "
    "input, and its chain became the primary/secondary/tertiary roles. The other chains are not "
    "exported. To pick another, set \"How to pick primary / secondary / tertiary objects\" to Manual "
    "and name the objects there.")
# fallback -> (module type that was ambiguous, role its output object took)
CHAIN = {"most_related_tertiary": ("IdentifyTertiaryObjects", "tertiary"),
         "most_related_secondary": ("IdentifySecondaryObjects", "secondary")}
METADATA = (
    "No Metadata_%s found. Add a Metadata module extracting Plate, Well and Site from file names "
    "(e.g. ^(?P<Plate>.*)_(?P<Well>[A-P][0-9]{2})_s(?P<Site>[0-9])) or LoadData columns "
    "Metadata_Plate/Well/Site; until then obs names are img<ImageNumber>_<label> and are not "
    "unique across runs.")
FILE_OBJECTS = (
    "Objects %s are loaded from files; parent/child links exist only if RelateObjects is used. "
    "Add a RelateObjects module if these objects should join onto the cell rows.")
SPREADSHEET = (
    "ExportToSpreadsheet is in this pipeline. For its CSVs to round-trip to the same object this "
    "module writes, set: Add image metadata columns: Yes; Add image file and folder names: Yes; "
    "Export all measurement types: Yes; no per-image aggregation; Representation of NaN: NaN; "
    "one file per object.")


def advice(ctx: Context) -> List[str]:
    """Every applicable advice message, most important first."""
    out = []
    fallback = ctx.role_note.get("fallback")
    candidates = ", ".join(ctx.role_note.get("candidates", []))
    if fallback == "most_related_primary":
        out.append(ROLE_FALLBACK % (candidates, ctx.roles.get("primary")))
    elif fallback in CHAIN:
        module_name, role = CHAIN[fallback]
        out.append(CHAIN_FALLBACK % (module_name, candidates, ctx.roles.get(role)))
    missing = [t for t in ("Plate", "Well", "Site") if t not in ctx.metadata_tags]
    if missing:
        out.append(METADATA % "/".join(missing))
    loaded = file_loaded_objects(ctx)
    if loaded:
        out.append(FILE_OBJECTS % ", ".join(loaded))
    if any(mod.get("name") == "ExportToSpreadsheet" for mod in ctx.modules):
        out.append(SPREADSHEET)
    return out
