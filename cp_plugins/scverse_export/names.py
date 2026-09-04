"""CellProfiler feature identity -> cp_measure / squidpy column names.

Rules are spec section 2c; they operate on the *structured* identity a CellProfiler module
reports through get_categories/get_measurements/get_measurement_images/get_measurement_scales,
never on regex over the CP name.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

CP_MEASURE = "cp_measure"
CELLPROFILER = "cellprofiler"

# Correlation measurements cp_measure emits with a _1/_2 direction suffix (canonical channel order).
_DIRECTIONAL = {"Manders", "RWC", "Costes"}
_NOT_IN_CPM = {"K", "KS", "Overlap"}

# AreaShape measurements that are pixel coordinates or an angle in the image frame, not a shape
# descriptor. MeasureObjectSizeShape reports these alongside Area/Eccentricity/etc. under the same
# category, so category alone can't separate them the way it does for Location.
_SPATIAL_AREASHAPE = {"Center_X", "Center_Y", "Orientation",
                      "BoundingBoxMinimum_X", "BoundingBoxMinimum_Y",
                      "BoundingBoxMaximum_X", "BoundingBoxMaximum_Y"}

# Neighbors measurements that name another object instead of measuring or relating to this one: a
# label with no biological content, the same problem as Parent_<object>. The rest of the category
# (NumberOfNeighbors, PercentTouching, {First,Second}ClosestDistance, AngleBetweenNeighbors) stays
# in X. Crowding, density, and the alignment between neighboring cells count as biology here, not
# imaging artifacts.
_NEIGHBOR_IDENTIFIERS = {"FirstClosestObjectNumber", "SecondClosestObjectNumber"}


@dataclass(frozen=True)
class Feature:
    object: str
    cp_name: str
    category: str
    measurement: str
    module_num: Optional[int] = None
    module_name: Optional[str] = None
    image: Optional[str] = None
    image2: Optional[str] = None
    other_object: Optional[str] = None
    scale: Optional[str] = None
    coltype: str = "float"
    parsed_by: str = "api"


def split_pair(pair: str, channels: Sequence[str]) -> Optional[Tuple[str, str]]:
    """Split 'A_B' into two known channel names; channel names may themselves contain '_'."""
    for a in sorted(channels, key=len, reverse=True):
        if pair.startswith(a + "_"):
            b = pair[len(a) + 1:]
            if b in channels:
                return a, b
    return None


def _cp_only(f: Feature) -> List[Tuple[str, str]]:
    return [(f.cp_name, CELLPROFILER)]


def is_extrinsic(f: Feature) -> bool:
    """True for measurements that don't describe the object's own morphology or intensity. Left in
    X, they'd bias downstream similarity on something other than biology.

    Two categories qualify:

    - Position and orientation: the whole Location category (an object's own center and, per
      channel, its intensity-weighted/max-intensity center, all absolute pixel coordinates) plus
      the handful of AreaShape measurements that are themselves a coordinate or an angle in the
      image frame (Center_X/Y, the bounding-box corners, Orientation).
    - Identity and linkage: Number_Object_Number (the object's own arbitrary label, already carried
      as obs["label_id"]), Parent_<object> (an arbitrary label referencing another row, not a
      measurement of this one), Children_<object>_Count (join_tables already surfaces this in obs
      as count_<child>, so it stays out of X rather than duplicating that value under a second,
      inconsistent name), and Neighbors_{First,Second}ClosestObjectNumber (another object's label,
      the same problem as Parent).

    Two objects with identical biology should land at the same point in feature space no matter
    where they were imaged, how the sample was rotated, or what label CellProfiler happened to
    assign them or their neighbors. ExportToAnnData keeps these measurements out of X/var and
    reports them in obs instead, prefixed per object the same way var_names are
    (`{object}__{name}`) once objects are joined.
    """
    if f.category == "Location":
        return True
    if f.category == "AreaShape":
        return f.measurement in _SPATIAL_AREASHAPE
    if f.category in ("Number", "Parent"):
        return True
    if f.category == "Children":
        return f.measurement.endswith("_Count")
    if f.category == "Neighbors":
        return f.measurement in _NEIGHBOR_IDENTIFIERS
    return False


def to_cpm_names(f: Feature, channels: Sequence[str]) -> List[Tuple[str, str]]:
    cat, meas, img = f.category, f.measurement, f.image
    if cat == "AreaShape":
        return [(meas, CP_MEASURE)]
    if cat == "Location":
        if meas in ("Center_X", "Center_Y"):
            return [(meas, CP_MEASURE)]
        if img and meas.startswith(("CenterMassIntensity_", "MaxIntensity_")):
            return [(f"Location_{meas}__{img}", CP_MEASURE)]
        return _cp_only(f)
    if cat == "Intensity" and img:
        return [(f"Intensity_{meas}__{img}", CP_MEASURE)]
    if cat == "Texture" and img and f.scale:
        return [(f"{meas}_{f.scale}__{img}", CP_MEASURE)]
    if cat == "RadialDistribution" and img and f.scale:
        return [(f"RadialDistribution_{meas}_{f.scale}__{img}", CP_MEASURE)]
    if cat == "Granularity" and img:
        return [(f"Granularity_{meas}__{img}", CP_MEASURE)]
    if cat == "Correlation" and img:
        pair = (img, f.image2) if f.image2 else split_pair(img, channels)
        if pair is None or meas in _NOT_IN_CPM:
            return _cp_only(f)
        a, b = pair
        if meas == "Correlation":
            return [(f"Correlation_Pearson__{a}__{b}", CP_MEASURE), (f"Correlation_Pearson__{b}__{a}", CP_MEASURE)]
        if meas == "Slope":
            return [(f"Correlation_Slope__{a}__{b}", CP_MEASURE)]
        if meas in _DIRECTIONAL:
            canonical = channels.index(a) < channels.index(b)
            first, second = (a, b) if canonical else (b, a)
            return [(f"Correlation_{meas}_{1 if canonical else 2}__{first}__{second}", CP_MEASURE)]
        return _cp_only(f)
    return _cp_only(f)
