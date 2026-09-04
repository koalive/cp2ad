"""The table-building step both export modules share.

ExportToAnnData and ExportForSpatialData need the same per-cell table: gather provenance, build one
table per role object, join them into one row per cell, attach the uns entries. Only what happens
to the result differs, since one writes a single .h5ad and the other writes one plate folder of
images, label arrays, and a table.

Building stops here and writing stays in the modules on purpose. Each module logs under its own
logger name, which its tests assert on, and the two have different ideas about where files go, so
a shared writer would have to invent a path protocol that suits neither.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

from .advice import advice as _advice
from .assemble import Table, build_object_table, join_tables, provenance, resolve_naming
from .introspect import Context
from .preview import channel_report, mapping_to_uns, measurement_report, object_report
from .samples import SampleNaming


@dataclass
class Export:
    """Everything assembled and ready to write, plus what the caller should tell the user.

    `advice` is returned rather than logged so each module logs it under its own logger. `naming`
    is the sample-key scheme this run resolved to, which ExportForSpatialData also names its
    images, label arrays, and coordinate systems from.
    """
    joined: Table
    per_object: Dict[str, Table] = field(default_factory=dict)
    provenance: Dict[str, Any] = field(default_factory=dict)
    advice: List[str] = field(default_factory=list)
    naming: SampleNaming = field(default_factory=SampleNaming)


def build_export(ctx: Context, m, policy: str = "flag",
                 exporter_settings: Dict[str, Any] = None,
                 wants_mapping_uns: bool = False,
                 sample_tags: Optional[Sequence[str]] = None) -> Export:
    """One row per cell across every image set in this run, with the per-object tables it was
    joined from. `uns["cellprofiler"]` carries the pipeline provenance on every table;
    `uns["cellprofiler_mapping"]` carries the channel, object and measurement tables when asked
    for, on the joined table only.

    `sample_tags` names the Metadata tags that form the sample key; None detects them. Resolved
    once here and passed down, so every object's table names its rows the same way.
    """
    naming = resolve_naming(ctx, m, sample_tags)
    prov = provenance(ctx, m, exporter_settings or {}, naming=naming)
    per_object = {obj: build_object_table(ctx, m, obj, naming=naming) for obj in ctx.roles.values()}
    joined = join_tables(ctx, m, per_object, policy=policy)
    joined.uns["cellprofiler"] = prov
    if wants_mapping_uns:
        joined.uns["cellprofiler_mapping"] = mapping_to_uns(
            channel_report(ctx), object_report(ctx), measurement_report(ctx))
    for table in per_object.values():
        table.uns["cellprofiler"] = prov
    return Export(joined=joined, per_object=per_object, provenance=prov,
                  advice=_advice(ctx), naming=naming)
