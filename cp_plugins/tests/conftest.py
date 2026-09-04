"""Shared fixtures. Headless must be set before any other cellprofiler_core import."""
import json
import os
import re
import sys

import cellprofiler_core.preferences

cellprofiler_core.preferences.set_headless()

import numpy
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PLUGIN_DIR = os.path.dirname(HERE)                      # cp_plugins/
SANDBOX = os.path.join(os.path.dirname(PLUGIN_DIR), "plugin_sandbox")
sys.path.insert(0, PLUGIN_DIR)                          # so `import scverse_export` and `import exporttoanndata` work


@pytest.fixture(scope="session")
def probe():
    with open(os.path.join(SANDBOX, "probe", "out", "full.json")) as fd:
        return json.load(fd)


@pytest.fixture(scope="session")
def cpm_columns():
    with open(os.path.join(SANDBOX, "cpm_alignment", "cpm_columns.json")) as fd:
        return json.load(fd)


@pytest.fixture(scope="session")
def mapping_cells():
    with open(os.path.join(SANDBOX, "cpm_alignment", "mapping_cells.json")) as fd:
        return json.load(fd)


@pytest.fixture(scope="session")
def meas_arrays():
    out = {}
    for obj in ("Nuclei", "Cells", "Cytoplasm", "PH3"):
        z = numpy.load(os.path.join(SANDBOX, "probe", "out", "arrays", f"meas_{obj}.npz"))
        out[obj] = {k: z[k] for k in z.files}
    return out


class FakeSetting:
    def __init__(self, text, value):
        self.text, self.value, self.value_text = text, value, value


def cppipe_settings(path):
    """Ordered (text, value) setting pairs per module_num, read back from the saved pipeline.

    probe/out/full.json stores each module's settings as a dict, which silently drops repeated
    setting texts (MeasureTexture writes one "Texture scale to measure" per scale). The real
    Module.settings() keeps them, so parse the .cppipe -- one "Text:value" line per setting.
    """
    out, num = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        header = re.match(r"^\w+:\[module_num:(\d+)\|", line)
        if header:
            num = int(header.group(1))
            out[num] = []
        elif num and line.startswith("    ") and ":" in line:
            text, _, value = line[4:].partition(":")
            out[num].append((text, value))
    return out


class FakeModule:
    """Replays the introspection API results recorded in full.json for one module."""

    def __init__(self, num, name, settings, entries, enabled=True, setting_pairs=None):
        self.module_num, self.module_name, self.enabled = num, name, enabled
        pairs = setting_pairs if setting_pairs is not None else list(settings.items())
        self._settings = [FakeSetting(k, v) for k, v in pairs]
        self._entries = entries  # list of {obj, cat, meas, images, objects, scales{img: [..]}}

    def settings(self):
        return self._settings

    def is_object_identification_module(self):
        return self.module_name in ("IdentifyPrimaryObjects", "IdentifySecondaryObjects", "IdentifyTertiaryObjects")

    def get_categories(self, pipeline, object_name):
        return sorted({e["cat"] for e in self._entries if e["obj"] == object_name})

    def get_measurements(self, pipeline, object_name, category):
        seen = []
        for e in self._entries:
            if e["obj"] == object_name and e["cat"] == category and e["meas"] not in seen:
                seen.append(e["meas"])
        return seen

    def _find(self, object_name, category, measurement):
        for e in self._entries:
            if (e["obj"], e["cat"], e["meas"]) == (object_name, category, measurement):
                return e
        return None

    def get_measurement_images(self, pipeline, object_name, category, measurement):
        e = self._find(object_name, category, measurement)
        return list(e["images"]) if e else []

    def get_measurement_objects(self, pipeline, object_name, category, measurement):
        e = self._find(object_name, category, measurement)
        return list(e["objects"]) if e else []

    def get_measurement_scales(self, pipeline, object_name, category, measurement, image_name):
        e = self._find(object_name, category, measurement)
        if not e:
            return []
        return list(e["scales"].get(str(image_name), []))


class FakePipeline:
    def __init__(self, modules, columns, providers_images, providers_objects):
        self._modules, self._columns = modules, columns
        self._prov = {"imagegroup": providers_images, "objectgroup": providers_objects}

    def modules(self, exclude_disabled=True):
        return [m for m in self._modules if m.enabled or not exclude_disabled]

    def get_measurement_columns(self, terminating_module=None):
        return [tuple(c) for c in self._columns]

    def get_provider_dictionary(self, groupname, module=None):
        # value shape mirrors cellprofiler_core: {name: [(module, setting), ...]}
        by_num = {m.module_num: m for m in self._modules}
        return {name: [(by_num[num], FakeSetting(text, name)) for _, num, text in lst]
                for name, lst in self._prov[groupname].items()}


@pytest.fixture(scope="session")
def fake_pipeline(probe):
    intro = probe["introspection"]  # keys "num:name"
    ordered = cppipe_settings(os.path.join(SANDBOX, "probe", "full.cppipe"))
    modules = []
    for mod in probe["modules"]:
        key = f"{mod['num']}:{mod['name']}"
        modules.append(FakeModule(mod["num"], mod["name"], mod["settings"], intro.get(key, []), mod["enabled"],
                                  setting_pairs=ordered.get(mod["num"])))
    return FakePipeline(modules, probe["measurement_columns"], probe["providers_images"], probe["providers_objects"])
