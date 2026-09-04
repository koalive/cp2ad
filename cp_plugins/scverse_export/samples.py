"""Identifiers built from an image set's Metadata tags: the sample key naming one field of view,
and the obs name naming one object in it.

Kept separate from names.py, which is about measurement names, because these name image sets and
rows rather than features. ExportForSpatialData names its images, label arrays, and coordinate
systems from the same sample key ExportToAnnData's obs names are built on, so obs_name() composes
on sample_key() rather than formatting tags a second time. Two conventions that have to agree are
better as one function called twice.

Which tags form the key is detected from the pipeline by default, the same way object roles are,
because pipelines spell these differently: Plate/Well/Site is the classic triple, but Opera and
Harmony write Row/Column/Field, and plenty of pipelines have only some of them. Detecting beats
demanding one spelling, and beats the old all-or-nothing behavior that fell back to the image
number whenever any of Plate, Well and Site was missing, discarding well information the pipeline
did have.

Whatever the tags, the resulting key has to identify each image set uniquely, or rows from two
fields of view collide. That is checked against the real values at run time rather than assumed,
and the image number is appended when the check fails.

Pure numpy, no CellProfiler and no other module in this package, so anything here can import it.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, Optional, Sequence, Tuple

import numpy

PREFIX = "Metadata_"
PLATE, WELL, SITE = "Metadata_Plate", "Metadata_Well", "Metadata_Site"
SAMPLE_TAGS = (PLATE, WELL, SITE)       # the classic triple, and the Manual default

AUTOMATIC, MANUAL = "Automatic", "Manual"

# Detection preferences, most specific first. Each group contributes at most one part of the key,
# so a pipeline carrying both Well and Row/Column uses Well rather than all three.
PLATE_TAGS = ("Plate", "Barcode", "PlateID", "PlateName")
WELL_TAGS = ("Well", "WellName")
ROW_TAGS = ("Row", "WellRow")
COLUMN_TAGS = ("Column", "Col", "WellColumn")
SITE_TAGS = ("Site", "Field", "FieldIndex", "Position")

# Deliberately never auto-detected as part of a key. Frame and Series index a z plane or a
# timepoint within one field of view, not the field itself, and Channel and FileLocation vary
# within an image set. Including any of them would either not disambiguate fields or would split
# one field across several keys. When they do make image sets non-unique, the image number gets
# appended instead, which is visible in the key rather than a silent guess.
NEVER_AUTO = ("Frame", "Series", "Channel", "FileLocation", "ChannelName", "Z", "T")


@dataclass(frozen=True)
class SampleNaming:
    """How to build sample keys for one run.

    tags is the Metadata_* names in key order, empty when nothing usable was found.
    with_image_number appends `img<n>`, either because there are no tags or because the tags do
    not tell the image sets apart. mode and note record why, for uns and the logs.
    """
    tags: Tuple[str, ...] = ()
    with_image_number: bool = True
    mode: str = AUTOMATIC
    note: str = ""

    @property
    def parts(self) -> Tuple[str, ...]:
        """The key's components in order, as a person would describe them."""
        return tuple(t[len(PREFIX):] if t.startswith(PREFIX) else t for t in self.tags) + \
            (("ImageNumber",) if self.with_image_number else ())


def format_tag(v) -> str:
    """One Metadata value as it should read in a key. Metadata values are coerced to float in obs,
    so a Site read back as 1.0 still has to spell "1"."""
    if v is None:
        return ""
    if isinstance(v, (float, numpy.floating)) and numpy.isfinite(v) and float(v).is_integer():
        return str(int(v))
    return str(v)


def qualify(tag: str) -> str:
    """A tag name with the Metadata_ prefix, accepting it either way round, so a user typing
    "Well,Field" and one typing "Metadata_Well,Metadata_Field" get the same thing."""
    tag = tag.strip()
    return tag if tag.startswith(PREFIX) else PREFIX + tag


def parse_tags(text: str) -> Tuple[str, ...]:
    """A comma-separated setting value as qualified tag names, ignoring blanks."""
    return tuple(qualify(part) for part in text.split(",") if part.strip())


def has_sample_tags(md_feats: Sequence[str]) -> bool:
    """Whether all of Plate, Well and Site are available. Kept because it is what the Manual
    default assumes, and because advice.py reports on the same three."""
    return all(f in md_feats for f in SAMPLE_TAGS)


def detect_sample_tags(md_feats: Iterable[str]) -> Tuple[Tuple[str, ...], str]:
    """(tags, note) for the tags that best identify a field of view in this pipeline.

    At most one plate part, one well part, and one site part, preferring the whole-well tag over
    Row plus Column. Returns empty tags when the pipeline carries none of them.
    """
    available = set(md_feats)

    def first(candidates):
        for name in candidates:
            if PREFIX + name in available:
                return PREFIX + name
        return None

    plate = first(PLATE_TAGS)
    well = first(WELL_TAGS)
    row, column = first(ROW_TAGS), first(COLUMN_TAGS)
    site = first(SITE_TAGS)

    tags = []
    if plate:
        tags.append(plate)
    if well:
        tags.append(well)
    elif row and column:
        tags.extend([row, column])
    if site:
        tags.append(site)

    if not tags:
        return (), ("no Plate, Well, Row/Column or Site tag found, so the image number is the only "
                    "thing telling image sets apart")
    return tuple(tags), "detected from the pipeline's Metadata tags: " + ", ".join(
        t[len(PREFIX):] for t in tags)


def keys_identify_image_sets(values: Dict[int, Dict[str, Any]], tags: Sequence[str]) -> bool:
    """Whether these tags give every image set a different key. False when two fields of view
    share one key, which would collide their rows, and when any key is blank."""
    if not tags:
        return False
    seen = set()
    for image_number, md in values.items():
        key = "_".join(format_tag(md.get(tag)) for tag in tags)
        if not key.strip("_") or key in seen:
            return False
        seen.add(key)
    return True


def resolve_sample_naming(md_feats: Sequence[str], values: Dict[int, Dict[str, Any]],
                          requested: Optional[Sequence[str]] = None) -> SampleNaming:
    """Pick the sample-key tags for this run and check them against the real values.

    `requested` is the Manual tag list; None means detect. Either way the tags have to identify
    the image sets, and the image number is appended when they do not, so a key is never
    ambiguous even if the tags were a poor choice.
    """
    if requested is not None:
        tags = tuple(requested)
        if not tags:
            return SampleNaming(tags=(), with_image_number=True, mode=MANUAL,
                                note="no tags set, so the image number names the field of view")
        missing = [t for t in tags if t not in md_feats]
        if missing:
            note = ("set manually, but this pipeline has no " +
                    ", ".join(t[len(PREFIX):] for t in missing) +
                    "; the image number is used instead")
            return SampleNaming(tags=(), with_image_number=True, mode=MANUAL, note=note)
        note = "set manually: " + ", ".join(t[len(PREFIX):] for t in tags)
        mode = MANUAL
    else:
        tags, note = detect_sample_tags(md_feats)
        mode = AUTOMATIC

    if not tags:
        return SampleNaming(tags=(), with_image_number=True, mode=mode, note=note)

    if keys_identify_image_sets(values, tags):
        return SampleNaming(tags=tags, with_image_number=False, mode=mode, note=note)
    return SampleNaming(
        tags=tags, with_image_number=True, mode=mode,
        note=note + ", plus the image number because those tags do not tell every image set apart")


def sample_key(md: Dict[str, Any], image_number: int, naming: SampleNaming) -> str:
    """One field of view's identity, e.g. `P1_A01_1` from Plate/Well/Site, `A02_03` from
    Well/Field, or `img3` when the pipeline has nothing usable. Unique within a run by
    construction, and unique across plates too when a plate tag is part of it, which is what lets
    two exports concatenate without renaming anything."""
    parts = [format_tag(md.get(tag)) for tag in naming.tags]
    if naming.with_image_number:
        parts.append(f"img{image_number}")
    return "_".join(parts)


def obs_name(md: Dict[str, Any], image_number: int, label: int, naming: SampleNaming) -> str:
    """One object's row name: its field of view's sample key, then the integer label CellProfiler
    gave it inside that image set."""
    return f"{sample_key(md, image_number, naming)}_{label}"
