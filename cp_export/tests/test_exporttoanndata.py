"""Unit tests for the ExportToAnnData module lifecycle. No CellProfiler app run here —
see test_integration_app.py for that."""
import logging

import pytest
from cellprofiler_core.preferences import ABSOLUTE_FOLDER_NAME
from cellprofiler_core.setting import ValidationError

from exporttoanndata import ROLE_MANUAL, ExportToAnnData

SETTINGS_ORDER = [
    "Output file location", "File name prefix", "Also write one file per object?",
    "When an object is not exactly one primary + one secondary + one tertiary",
    "How to pick primary / secondary / tertiary objects", "Primary objects (e.g. nuclei)",
    "Secondary objects (e.g. cells)", "Tertiary objects (e.g. cytoplasm)",
    "Overwrite existing files without warning?",
]


def _module():
    return ExportToAnnData()


class _Workspace:
    def __init__(self, pipeline):
        self.pipeline = pipeline


class _NoObjectPipeline:
    """A pipeline that identifies nothing: build_context can find no object to export."""

    def modules(self, exclude_disabled=True):
        return []

    def get_measurement_columns(self, terminating_module=None):
        return []

    def get_provider_dictionary(self, groupname, module=None):
        return {}


class _TinyPipeline:
    """Just enough of the pipeline API for validate_module's 'must be last' check."""

    def __init__(self, mods):
        self._mods = mods

    def modules(self):
        return self._mods


def test_roles_automatic_is_none():
    assert _module()._roles() is None


def test_roles_manual_all_none_is_none():
    m = _module()
    m.role_mode.value = ROLE_MANUAL
    assert m._roles() is None


def test_roles_manual_primary_only():
    m = _module()
    m.role_mode.value = ROLE_MANUAL
    m.primary_object.value = "PH3"
    assert m._roles() == {"primary": "PH3"}


def test_validate_module_raises_when_not_last():
    m = _module()
    pipeline = _TinyPipeline([m, object()])
    with pytest.raises(ValidationError):
        m.validate_module(pipeline)


def test_validate_module_raises_on_role_error(fake_pipeline, monkeypatch):
    m = _module()
    m.role_mode.value = ROLE_MANUAL
    m.primary_object.value = "Nope"
    orig_modules = fake_pipeline.modules
    monkeypatch.setattr(fake_pipeline, "modules",
                        lambda exclude_disabled=True: list(orig_modules(exclude_disabled)) + [m])
    with pytest.raises(ValidationError):
        m.validate_module(fake_pipeline)


def test_validate_module_warnings_flags_missing_metadata(fake_pipeline):
    m = _module()
    with pytest.raises(ValidationError, match="Metadata_Plate/Well/Site"):
        m.validate_module_warnings(fake_pipeline)


def test_validate_module_warnings_swallows_role_errors(fake_pipeline):
    """A bad role setting is validate_module's ValidationError; the warning pass must not also
    let the RoleError escape as a traceback out of the GUI validator."""
    m = _module()
    m.role_mode.value = ROLE_MANUAL
    m.primary_object.value = "Nope"
    assert m.validate_module_warnings(fake_pipeline) is None


def test_prepare_run_raises_when_overwrite_off_and_file_exists(fake_pipeline, tmp_path):
    m = _module()
    m.directory.dir_choice = ABSOLUTE_FOLDER_NAME
    m.directory.custom_path = str(tmp_path)
    m.wants_overwrite.value = False
    (tmp_path / "cellprofiler.h5ad").write_bytes(b"")
    with pytest.raises(ValueError):
        m.prepare_run(_Workspace(fake_pipeline))


def test_prepare_run_returns_true_when_overwrite_on(fake_pipeline, tmp_path):
    m = _module()
    m.directory.dir_choice = ABSOLUTE_FOLDER_NAME
    m.directory.custom_path = str(tmp_path)
    m.wants_overwrite.value = True
    (tmp_path / "cellprofiler.h5ad").write_bytes(b"")
    assert m.prepare_run(_Workspace(fake_pipeline)) is True


def test_prepare_run_refuses_cleanly_when_roles_cannot_be_detected(caplog, tmp_path):
    """A pipeline with no objects: prepare_run logs the reason and returns False, so CellProfiler
    aborts the run instead of showing a RoleError traceback."""
    m = _module()
    m.directory.dir_choice = ABSOLUTE_FOLDER_NAME
    m.directory.custom_path = str(tmp_path)
    with caplog.at_level(logging.ERROR, logger="exporttoanndata"):
        assert m.prepare_run(_Workspace(_NoObjectPipeline())) is False
    assert any(r.levelno == logging.ERROR and "no objects" in r.getMessage() for r in caplog.records)


def test_settings_order_is_frozen():
    assert [s.text for s in _module().settings()] == SETTINGS_ORDER


def test_apply_autoconfig_sets_manual_roles_and_explains(fake_pipeline):
    import exporttoanndata as e
    module = e.ExportToAnnData()
    text, warned = module.apply_autoconfig(fake_pipeline)
    assert module.role_mode.value == e.ROLE_MANUAL
    assert (module.primary_object.value, module.secondary_object.value, module.tertiary_object.value) == \
        ("Nuclei", "Cells", "Cytoplasm")
    assert 'primary = "Nuclei"' in text and "IdentifyPrimaryObjects" in text and "Manual" in text
    assert warned and "Metadata_Plate" in text  # ExampleHuman has no plate metadata -> advice repeated


def test_apply_autoconfig_ignores_previous_manual_values(fake_pipeline):
    import exporttoanndata as e
    module = e.ExportToAnnData()
    module.role_mode.value = e.ROLE_MANUAL
    module.primary_object.value = "PH3"
    module.apply_autoconfig(fake_pipeline)
    assert module.primary_object.value == "Nuclei"


def test_apply_autoconfig_no_objects_raises_roleerror():
    import exporttoanndata as e
    from conftest import FakePipeline
    from cpexport.introspect import RoleError
    import pytest as _pytest
    module = e.ExportToAnnData()
    with _pytest.raises(RoleError):
        module.apply_autoconfig(FakePipeline([], [], {}, {}))


def test_autoconfig_button_visible_but_not_persisted():
    import exporttoanndata as e
    module = e.ExportToAnnData()
    assert module.autoconfig not in module.settings()
    assert module.autoconfig in module.visible_settings()
