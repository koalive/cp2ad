"""Sample keys and obs names. Needs neither cellprofiler_core nor the recorded fixtures."""
import itertools

import numpy
import pytest

from scverse_export.samples import (AUTOMATIC, MANUAL, PLATE, SAMPLE_TAGS, SampleNaming,
                                    detect_sample_tags, format_tag, has_sample_tags,
                                    keys_identify_image_sets, obs_name, parse_tags, qualify,
                                    resolve_sample_naming, sample_key)

CLASSIC = SampleNaming(tags=SAMPLE_TAGS, with_image_number=False)
IMAGE_ONLY = SampleNaming(tags=(), with_image_number=True)


def md(plate="P1", well="A01", site=1):
    return {"Metadata_Plate": plate, "Metadata_Well": well, "Metadata_Site": site}


# ---- formatting ------------------------------------------------------------------------------

@pytest.mark.parametrize("value,expected", [
    (None, ""), ("P1", "P1"), ("A01", "A01"), (1, "1"),
    (1.0, "1"),                          # obs coerces metadata to float; a Site of 1.0 spells "1"
    (2.5, "2.5"), (numpy.float64(3.0), "3"), (numpy.float32(4.0), "4"), (numpy.int64(7), "7"),
    (float("nan"), "nan"),               # not finite, so not an integer either
    (float("inf"), "inf"), (0.0, "0"), ("", ""),
])
def test_format_tag(value, expected):
    assert format_tag(value) == expected


@pytest.mark.parametrize("given,expected", [
    ("Well", "Metadata_Well"),
    ("Metadata_Well", "Metadata_Well"),
    ("  Well  ", "Metadata_Well"),
])
def test_qualify_accepts_the_tag_either_way_round(given, expected):
    assert qualify(given) == expected


def test_parse_tags():
    assert parse_tags("Well,Field") == ("Metadata_Well", "Metadata_Field")
    assert parse_tags("Metadata_Plate, Well , Site") == (PLATE, "Metadata_Well", "Metadata_Site")
    assert parse_tags("") == ()
    assert parse_tags(" , ") == ()


# ---- detection -------------------------------------------------------------------------------

@pytest.mark.parametrize("available,expected", [
    # the classic triple
    (["Metadata_Plate", "Metadata_Well", "Metadata_Site"],
     ("Metadata_Plate", "Metadata_Well", "Metadata_Site")),
    # Opera and Harmony: Field is the site, and Well beats Row plus Column
    (["Metadata_Well", "Metadata_Row", "Metadata_Column", "Metadata_Field"],
     ("Metadata_Well", "Metadata_Field")),
    # no whole-well tag, so Row and Column together stand in for it
    (["Metadata_Row", "Metadata_Column", "Metadata_Field"],
     ("Metadata_Row", "Metadata_Column", "Metadata_Field")),
    # a row without a column is not a well, so neither is used
    (["Metadata_Row", "Metadata_Site"], ("Metadata_Site",)),
    # plate synonyms
    (["Metadata_Barcode", "Metadata_Well"], ("Metadata_Barcode", "Metadata_Well")),
    # well with no site at all
    (["Metadata_Plate", "Metadata_Well"], ("Metadata_Plate", "Metadata_Well")),
    # nothing usable
    ([], ()),
    (["Metadata_Frame", "Metadata_Series", "Metadata_Channel"], ()),
])
def test_detect_sample_tags(available, expected):
    tags, note = detect_sample_tags(available)
    assert tags == expected
    assert note


def test_detection_never_uses_frame_series_or_channel():
    """Frame and Series index a z plane or timepoint inside one field of view, and Channel varies
    within an image set. Using any of them would split one field across several keys."""
    tags, _ = detect_sample_tags(["Metadata_Well", "Metadata_Frame", "Metadata_Series",
                                  "Metadata_Channel", "Metadata_FileLocation"])
    assert tags == ("Metadata_Well",)


def test_detection_matches_the_example_pipeline():
    """The repo's test pipeline: Well plus Field are present, Plate and Site are not. The old
    all-or-nothing rule threw the well away and named rows img1_1; this keeps it."""
    available = ["Metadata_Column", "Metadata_Field", "Metadata_Frame", "Metadata_Row",
                 "Metadata_Series", "Metadata_Well"]
    tags, _ = detect_sample_tags(available)
    assert tags == ("Metadata_Well", "Metadata_Field")
    values = {1: {"Metadata_Well": "A02", "Metadata_Field": "03"},
              2: {"Metadata_Well": "B04", "Metadata_Field": "01"},
              3: {"Metadata_Well": "C07", "Metadata_Field": "05"},
              4: {"Metadata_Well": "G05", "Metadata_Field": "03"}}
    naming = resolve_sample_naming(available, values)
    assert naming.with_image_number is False
    assert sample_key(values[1], 1, naming) == "A02_03"
    assert obs_name(values[1], 1, 1, naming) == "A02_03_1"


# ---- uniqueness ------------------------------------------------------------------------------

def test_keys_identify_image_sets():
    unique = {1: {"Metadata_Well": "A01"}, 2: {"Metadata_Well": "A02"}}
    clashing = {1: {"Metadata_Well": "A01"}, 2: {"Metadata_Well": "A01"}}
    assert keys_identify_image_sets(unique, ("Metadata_Well",)) is True
    assert keys_identify_image_sets(clashing, ("Metadata_Well",)) is False
    assert keys_identify_image_sets(unique, ()) is False


def test_blank_values_do_not_count_as_identifying():
    """A tag that exists but is empty for every image set names nothing."""
    values = {1: {"Metadata_Well": None}, 2: {"Metadata_Well": ""}}
    assert keys_identify_image_sets(values, ("Metadata_Well",)) is False


def test_image_number_is_appended_when_tags_do_not_separate_image_sets():
    """Two fields of view in one well with no site tag. Without the image number their rows would
    collide, so it goes in, and the key shows that it happened."""
    available = ["Metadata_Well"]
    values = {1: {"Metadata_Well": "A01"}, 2: {"Metadata_Well": "A01"}}
    naming = resolve_sample_naming(available, values)
    assert naming.tags == ("Metadata_Well",)
    assert naming.with_image_number is True
    assert sample_key(values[1], 1, naming) == "A01_img1"
    assert sample_key(values[2], 2, naming) == "A01_img2"
    assert "do not tell every image set apart" in naming.note


# ---- resolution ------------------------------------------------------------------------------

def test_resolve_automatic():
    values = {1: md(), 2: md(site=2)}
    naming = resolve_sample_naming(list(SAMPLE_TAGS), values)
    assert naming.mode == AUTOMATIC
    assert naming.tags == SAMPLE_TAGS
    assert naming.with_image_number is False


def test_resolve_manual():
    values = {1: {"Metadata_Well": "A01"}, 2: {"Metadata_Well": "A02"}}
    naming = resolve_sample_naming(["Metadata_Well", "Metadata_Row"], values,
                                   requested=("Metadata_Well",))
    assert naming.mode == MANUAL
    assert naming.tags == ("Metadata_Well",)
    assert "set manually" in naming.note


def test_resolve_manual_with_a_tag_the_pipeline_lacks():
    """Naming a tag that is not there is a mistake worth reporting rather than crashing on, so it
    falls back to the image number and says which tag was missing."""
    values = {1: {"Metadata_Well": "A01"}}
    naming = resolve_sample_naming(["Metadata_Well"], values, requested=(PLATE, "Metadata_Well"))
    assert naming.tags == ()
    assert naming.with_image_number is True
    assert "no Plate" in naming.note


def test_resolve_with_no_usable_tags():
    naming = resolve_sample_naming([], {1: {}})
    assert naming.tags == () and naming.with_image_number is True
    assert sample_key({}, 3, naming) == "img3"


# ---- keys and row names ----------------------------------------------------------------------

def test_sample_key_and_obs_name():
    assert sample_key(md(), 3, CLASSIC) == "P1_A01_1"
    assert sample_key(md(site=1.0), 3, CLASSIC) == "P1_A01_1"
    assert obs_name(md(), 3, 5, CLASSIC) == "P1_A01_1_5"
    assert sample_key({}, 7, IMAGE_ONLY) == "img7"
    assert obs_name({}, 7, 9, IMAGE_ONLY) == "img7_9"


def test_obs_name_composes_on_sample_key():
    """The point of this module: element names and row names come from one function, so they
    cannot drift apart."""
    for naming in (CLASSIC, IMAGE_ONLY, SampleNaming(tags=("Metadata_Well",), with_image_number=True)):
        assert obs_name(md(), 3, 5, naming) == f"{sample_key(md(), 3, naming)}_5"


def test_classic_triple_still_produces_the_same_strings_as_before():
    """Detection changed what a pipeline missing Plate or Site is named, on purpose. A pipeline
    carrying all three has to be unaffected, or every row in every existing export is renamed.
    The reference is a verbatim copy of the code this replaced."""
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
        for n, label in [(1, 1), (3, 42), (7, 9)]:
            assert obs_name(tags, n, label, CLASSIC) == reference(tags, n, label, True)
            assert obs_name(tags, n, label, IMAGE_ONLY) == reference(tags, n, label, False)


# ---- reporting -------------------------------------------------------------------------------

def test_parts_describes_the_key_in_readable_terms():
    assert CLASSIC.parts == ("Plate", "Well", "Site")
    assert IMAGE_ONLY.parts == ("ImageNumber",)
    assert SampleNaming(tags=("Metadata_Well",), with_image_number=True).parts == ("Well", "ImageNumber")


@pytest.mark.parametrize("feats,expected", [
    (list(SAMPLE_TAGS), True),
    (list(SAMPLE_TAGS) + ["Metadata_Frame"], True),
    (["Metadata_Plate", "Metadata_Well"], False),
    ([], False),
])
def test_has_sample_tags(feats, expected):
    assert has_sample_tags(feats) is expected


def test_sample_tags_order_is_plate_well_site():
    """Order is load-bearing: it is the order in every existing obs name built from the triple."""
    assert SAMPLE_TAGS == ("Metadata_Plate", "Metadata_Well", "Metadata_Site")
    assert PLATE == "Metadata_Plate"


# ---- what "Auto-configure from this pipeline" pins -------------------------------------------

@pytest.mark.parametrize("available,expected_setting", [
    (["Well", "Field", "Row", "Column", "Frame", "Series"], "Well,Field"),
    (["Plate", "Well", "Site"], "Plate,Well,Site"),
    (["Row", "Column", "Field"], "Row,Column,Field"),
    (["Barcode", "Well"], "Barcode,Well"),
    (["Frame", "Series"], ""),
    ([], ""),
])
def test_autoconfig_pins_a_tag_list_that_reproduces_automatic(available, expected_setting):
    """The Auto-configure button writes detected tags into the Manual setting and switches to
    Manual, so the export stops depending on the pipeline at run time. Pinning is only useful if
    it reproduces what Automatic would have chosen, which is what this checks: the same two steps
    the button runs, then what the run resolves from the value it wrote.
    """
    from scverse_export.samples import PREFIX
    md_feats = [PREFIX + tag for tag in available]

    detected, _ = detect_sample_tags(md_feats)
    setting_value = ",".join(t[len(PREFIX):] for t in detected)
    assert setting_value == expected_setting

    values = {1: {PREFIX + tag: f"{tag}1" for tag in available}}
    pinned = resolve_sample_naming(md_feats, values, requested=parse_tags(setting_value))
    automatic = resolve_sample_naming(md_feats, values)
    assert pinned.mode == MANUAL
    assert pinned.tags == automatic.tags
    assert pinned.with_image_number == automatic.with_image_number


def test_manual_with_no_tags_names_rows_by_image_number():
    """What the button writes when a pipeline has nothing usable. An empty list is a real answer,
    not a broken setting, so it gets its own message rather than a dangling one."""
    naming = resolve_sample_naming(["Metadata_Frame"], {1: {}}, requested=())
    assert naming.tags == () and naming.with_image_number is True
    assert naming.mode == MANUAL
    assert naming.note == "no tags set, so the image number names the field of view"
    assert sample_key({}, 4, naming) == "img4"
