"""Identifiers built from an image set's Metadata tags: the sample key for one field of view, and
the obs name for one object in it.

Kept separate from names.py, which is about measurement names, because these name image sets and
rows rather than features. ExportForSpatialData names its images, label arrays, and coordinate
systems from the same sample key that ExportToAnnData's obs names are built on, so obs_name()
composes on sample_key() rather than formatting the tags a second time. Two conventions that have
to agree are better as one function called twice.

Pure numpy, no CellProfiler and no other module in this package, so anything here can import it.
"""
from __future__ import annotations

from typing import Any, Dict, Sequence

import numpy

PLATE, WELL, SITE = "Metadata_Plate", "Metadata_Well", "Metadata_Site"
SAMPLE_TAGS = (PLATE, WELL, SITE)


def format_tag(v) -> str:
    """One Metadata value as it should read in a key. Metadata values are coerced to float in obs,
    so a Site read back as 1.0 still has to spell "1"."""
    if v is None:
        return ""
    if isinstance(v, (float, numpy.floating)) and numpy.isfinite(v) and float(v).is_integer():
        return str(int(v))
    return str(v)


def has_sample_tags(md_feats: Sequence[str]) -> bool:
    """Whether Plate, Well and Site are all available, which is what decides between a readable
    sample key and the image-number fallback."""
    return all(f in md_feats for f in SAMPLE_TAGS)


def sample_key(md: Dict[str, Any], image_number: int, tags_ok: bool) -> str:
    """One field of view's identity: `<Plate>_<Well>_<Site>`, or `img<n>` when the pipeline does
    not extract all three tags. Unique across plates when the tags are there, which is what lets
    two exports concatenate without renaming anything."""
    if tags_ok:
        return "_".join(format_tag(md[tag]) for tag in SAMPLE_TAGS)
    return f"img{image_number}"


def obs_name(md: Dict[str, Any], image_number: int, label: int, tags_ok: bool) -> str:
    """One object's row name: its field of view's sample key plus the integer label CellProfiler
    gave it inside that image set."""
    return f"{sample_key(md, image_number, tags_ok)}_{label}"
