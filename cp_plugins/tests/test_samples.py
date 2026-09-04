"""Sample keys and obs names. Needs neither cellprofiler_core nor the recorded fixtures."""
import itertools

import numpy
import pytest

from scverse_export.samples import PLATE, SAMPLE_TAGS, format_tag, has_sample_tags, obs_name, sample_key


def md(plate="P1", well="A01", site=1):
    return {"Metadata_Plate": plate, "Metadata_Well": well, "Metadata_Site": site}


@pytest.mark.parametrize("value,expected", [
    (None, ""),
    ("P1", "P1"),
    ("A01", "A01"),
    (1, "1"),
    (1.0, "1"),                          # obs coerces metadata to float; a Site of 1.0 spells "1"
    (2.5, "2.5"),
    (numpy.float64(3.0), "3"),
    (numpy.float32(4.0), "4"),
    (numpy.int64(7), "7"),
    (float("nan"), "nan"),               # not finite, so it is not an integer either
    (float("inf"), "inf"),
    (0.0, "0"),
    ("", ""),
])
def test_format_tag(value, expected):
    assert format_tag(value) == expected


def test_sample_key_from_the_three_tags():
    assert sample_key(md(), 3, True) == "P1_A01_1"
    assert sample_key(md(site=1.0), 3, True) == "P1_A01_1"
    assert sample_key(md(plate="Plate2", well="H12", site=4), 9, True) == "Plate2_H12_4"


def test_sample_key_falls_back_to_the_image_number():
    """No Plate/Well/Site means no readable key, so the image number stands in. It is unique
    within a run but not across runs, which is what makes the tags worth having."""
    assert sample_key(md(), 3, False) == "img3"
    assert sample_key({}, 7, False) == "img7"


def test_obs_name_composes_on_sample_key():
    """The point of this module: element names and row names come from one function, so they
    cannot drift apart."""
    for tags_ok in (True, False):
        key = sample_key(md(), 3, tags_ok)
        assert obs_name(md(), 3, 5, tags_ok) == f"{key}_5"


def test_obs_name_exact_strings():
    assert obs_name(md(), 3, 5, True) == "P1_A01_1_5"
    assert obs_name(md(), 3, 5, False) == "img3_5"


def test_obs_name_matches_the_implementation_it_replaced():
    """Guards the extraction from assemble._obs_name. The reference below is a verbatim copy of
    the code that was removed, so a change in either branch shows up here rather than silently
    renaming every row in every existing export."""
    def reference(md_, n, label, tags_ok):
        def fmt(v):
            if v is None:
                return ""
            if isinstance(v, (float, numpy.floating)) and numpy.isfinite(v) and float(v).is_integer():
                return str(int(v))
            return str(v)
        if tags_ok:
            return "_".join([fmt(md_["Metadata_Plate"]), fmt(md_["Metadata_Well"]),
                             fmt(md_["Metadata_Site"]), str(label)])
        return f"img{n}_{label}"

    values = [None, "P1", "A01", 1, 1.0, 2.5, numpy.float64(3.0), numpy.float32(4.0),
              float("nan"), numpy.int64(7), "", 0]
    for plate, well, site in itertools.product(values, values[:5], values):
        tags = {"Metadata_Plate": plate, "Metadata_Well": well, "Metadata_Site": site}
        for n, label, tags_ok in [(1, 1, True), (3, 42, True), (1, 1, False), (7, 9, False)]:
            assert obs_name(tags, n, label, tags_ok) == reference(tags, n, label, tags_ok)


@pytest.mark.parametrize("feats,expected", [
    (list(SAMPLE_TAGS), True),
    (list(SAMPLE_TAGS) + ["Metadata_Frame"], True),
    (["Metadata_Plate", "Metadata_Well"], False),
    (["Metadata_Well", "Metadata_Site"], False),
    (["Metadata_Site"], False),
    ([], False),
])
def test_has_sample_tags(feats, expected):
    assert has_sample_tags(feats) is expected


def test_sample_tags_order_is_plate_well_site():
    """Order is load-bearing: it is the order in every existing obs name, and reordering would
    rename every row."""
    assert SAMPLE_TAGS == ("Metadata_Plate", "Metadata_Well", "Metadata_Site")
    assert PLATE == "Metadata_Plate"
