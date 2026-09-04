"""Where an export module sits in the pipeline. Needs neither cellprofiler_core nor the recorded
fixtures, because placement_advice only calls get_measurement_columns and modules().

Being last used to be a hard ValidationError in both modules, which was wrong: ExportToAnnData
followed by ExportForSpatialData is a normal pipeline, and so is either followed by SaveImages. The
check is advice now, and it only fires for a module that actually makes measurements.
"""
import pytest

from scverse_export.advice import (PLACEMENT_TABLE_ONLY, PLACEMENT_WITH_PIXELS, placement_advice)


class FakeModule:
    """Enough of a module for placement_advice: a number, a name, and its declared columns."""

    def __init__(self, module_num, module_name, columns=()):
        self.module_num = module_num
        self.module_name = module_name
        self._columns = list(columns)

    def get_measurement_columns(self, pipeline):
        return self._columns


class SloppyModule(FakeModule):
    """A third-party module whose get_measurement_columns raises. Validation must survive it."""

    def get_measurement_columns(self, pipeline):
        raise RuntimeError("this module's override is broken")


class FakePipeline:
    def __init__(self, modules):
        self._modules = list(modules)

    def modules(self):
        return self._modules


FEATURE = [("Cells", "Texture_Contrast", "float")]


def test_no_advice_when_last():
    exporter = FakeModule(1, "ExportToAnnData")
    assert placement_advice(exporter, FakePipeline([exporter])) == []


@pytest.mark.parametrize("name", ["ExportForSpatialData", "ExportToSpreadsheet", "SaveImages",
                                  "DisplayDataOnImage", "CreateBatchFiles"])
def test_no_advice_for_a_later_module_that_makes_no_measurements(name):
    """The case that used to be a blocking error. A module declaring no columns after an exporter
    is ordinary, and chaining the two export modules is the obvious example."""
    exporter = FakeModule(1, "ExportToAnnData")
    assert placement_advice(exporter, FakePipeline([exporter, FakeModule(2, name)])) == []


def test_advice_names_the_later_measuring_modules():
    exporter = FakeModule(1, "ExportToAnnData")
    pipeline = FakePipeline([exporter,
                             FakeModule(2, "MeasureTexture", FEATURE),
                             FakeModule(3, "MeasureGranularity", FEATURE)])
    messages = placement_advice(exporter, pipeline)
    assert len(messages) == 1
    assert "#2 MeasureTexture" in messages[0]
    assert "#3 MeasureGranularity" in messages[0]


def test_a_measuring_module_before_the_exporter_is_fine():
    """Only what comes after matters, which is the whole point."""
    exporter = FakeModule(2, "ExportToAnnData")
    pipeline = FakePipeline([FakeModule(1, "MeasureTexture", FEATURE), exporter])
    assert placement_advice(exporter, pipeline) == []


def test_a_broken_later_module_is_skipped_not_raised():
    """A third-party module with a sloppy override must not break validation for everyone else."""
    exporter = FakeModule(1, "ExportToAnnData")
    pipeline = FakePipeline([exporter, SloppyModule(2, "Sloppy")])
    assert placement_advice(exporter, pipeline) == []

    pipeline = FakePipeline([exporter, SloppyModule(2, "Sloppy"),
                             FakeModule(3, "MeasureTexture", FEATURE)])
    assert "#3 MeasureTexture" in placement_advice(exporter, pipeline)[0]


def test_a_module_not_in_the_pipeline_gets_no_advice():
    """Validation can run against a module that is not part of the pipeline it was handed."""
    assert placement_advice(FakeModule(1, "ExportToAnnData"),
                            FakePipeline([FakeModule(2, "MeasureTexture", FEATURE)])) == []


@pytest.mark.parametrize("consequence,expected", [
    (PLACEMENT_TABLE_ONLY, "features may not be exported"),
    (PLACEMENT_WITH_PIXELS, "written before they ran"),
])
def test_the_consequence_differs_between_the_two_modules(consequence, expected):
    """ExportToAnnData only risks missing features. ExportForSpatialData also writes pixels at its
    own position, so a later module changing an image leaves the arrays on disk stale, which is a
    worse and quieter failure and says so."""
    exporter = FakeModule(1, "Exporter")
    pipeline = FakePipeline([exporter, FakeModule(2, "MeasureTexture", FEATURE)])
    assert expected in placement_advice(exporter, pipeline, consequence)[0]


def test_advice_says_how_to_fix_it():
    exporter = FakeModule(1, "Exporter")
    pipeline = FakePipeline([exporter, FakeModule(2, "MeasureTexture", FEATURE)])
    assert "Move this module last" in placement_advice(exporter, pipeline)[0]


def test_advice_is_two_sentences():
    """It appears in a validation popup next to other advice, so length matters. The old version
    explained at the user that a warning is not an error, which is meta-talk they can see for
    themselves from the fact that nothing stopped."""
    exporter = FakeModule(1, "Exporter")
    pipeline = FakePipeline([exporter, FakeModule(2, "MeasureTexture", FEATURE)])
    message = placement_advice(exporter, pipeline)[0]
    assert len(message) < 160, message
    assert "advice rather than an error" not in message
