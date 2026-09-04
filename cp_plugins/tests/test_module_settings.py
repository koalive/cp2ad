"""Invariants about the two top-level module files themselves.

These are checked by reading the source rather than by importing it, because importing needs
cellprofiler_core. That is a real limitation, but it also catches the class of bug that put this
file here: an assignment whose type is wrong only blows up at GUI time, where an exception in a
hook like on_activated is a freeze rather than a message.
"""
import ast
import pathlib

import pytest

MODULE_DIR = pathlib.Path(__file__).resolve().parent.parent
MODULE_FILES = sorted(p for p in MODULE_DIR.glob("*.py"))
LIST_SETTINGS = ("ImageListSubscriber", "LabelListSubscriber")


def _tree(path):
    return ast.parse(path.read_text())


def _setting_types(tree):
    """{attribute name: setting class} for every self.<name> = SomeSetting(...) in the file."""
    types = {}
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign) and isinstance(node.value, ast.Call)
                and isinstance(node.value.func, ast.Name)):
            for target in node.targets:
                if isinstance(target, ast.Attribute) and isinstance(target.value, ast.Name) \
                        and target.value.id == "self":
                    types[target.attr] = node.value.func.id
    return types


@pytest.mark.parametrize("path", MODULE_FILES, ids=lambda p: p.name)
def test_file_defines_exactly_one_module_class(path):
    """CellProfiler's loader takes the first Module subclass in a top-level file and warns if it
    finds none, so one class per file is the contract, not a style choice."""
    tree = _tree(path)
    classes = [n for n in tree.body if isinstance(n, ast.ClassDef)]
    assert len(classes) == 1, f"{path.name} defines {[c.name for c in classes]}"
    bases = [getattr(b, "id", getattr(b, "attr", None)) for b in classes[0].bases]
    assert bases == ["Module"], bases


@pytest.mark.parametrize("path", MODULE_FILES, ids=lambda p: p.name)
def test_list_settings_are_assigned_strings_not_lists(path):
    """The bug this file exists for.

    ImageListSubscriber and LabelListSubscriber read a list back out of .value, but their setter
    takes the saved string form and calls .split(", ") on it. Assigning a list therefore raises
    AttributeError: 'list' object has no attribute 'split'. In on_activated that is one exception
    per repaint, which presents as the GUI hanging rather than as an error anyone can read.

    So every assignment to a list setting has to produce a string. A str.join call is the
    idiomatic way and the only form allowed here; anything else needs a look.
    """
    tree = _tree(path)
    types = _setting_types(tree)
    offenders = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if not (isinstance(target, ast.Attribute) and target.attr == "value"
                    and isinstance(target.value, ast.Attribute)):
                continue
            name = target.value.attr
            if types.get(name) not in LIST_SETTINGS:
                continue
            rhs = node.value
            joined = (isinstance(rhs, ast.Call) and isinstance(rhs.func, ast.Attribute)
                      and rhs.func.attr == "join")
            literal_str = isinstance(rhs, ast.Constant) and isinstance(rhs.value, str)
            if not (joined or literal_str):
                offenders.append(f"line {node.lineno}: self.{name}.value = "
                                 f"{ast.dump(rhs)[:60]}")
    assert not offenders, (
        f"{path.name}: list settings must be assigned a string, e.g. \", \".join(names). "
        f"Offenders: {offenders}")


@pytest.mark.parametrize("path", MODULE_FILES, ids=lambda p: p.name)
def test_gui_hooks_cannot_raise_out(path):
    """on_activated runs on every repaint of the module, so anything it calls that is a
    convenience rather than a necessity has to be wrapped. This is what turns a mistake into a
    message instead of a freeze."""
    tree = _tree(path)
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "on_activated":
            calls = [n for n in ast.walk(node) if isinstance(n, ast.Call)]
            if len(calls) <= 1:          # only the assignment of the pipeline; nothing to guard
                return
            handlers = [n for n in ast.walk(node) if isinstance(n, ast.Try)]
            assert handlers, (f"{path.name}: on_activated calls {len(calls)} things but has no "
                              "try/except; an exception here is one per repaint")
            return


def test_the_list_setting_contract_this_guards():
    """Why the test above exists, as executable documentation.

    A faithful copy of ListSubscriber's value property from cellprofiler_core v4.2.8.1
    (_list_subscriber.py). Kept here because the real class cannot be imported without
    cellprofiler_core, and the asymmetry is surprising enough to be worth pinning: the getter
    hands back a list, the setter wants the string.
    """
    class ListSubscriber:
        def __init__(self, value=""):
            self.__value = []
            self.value = value

        def _get(self):
            return self.__value

        def _set(self, value):
            self.__value = [] if len(value) == 0 else value.split(", ")

        value = property(_get, _set)

        @property
        def value_text(self):
            return ", ".join(map(str, self.__value))

    setting = ListSubscriber()
    assert setting.value == []

    with pytest.raises(AttributeError, match="split"):
        setting.value = ["DNA", "Tubulin"]              # what the freeze was

    setting.value = ", ".join(["DNA", "Tubulin"])        # what it has to be
    assert setting.value == ["DNA", "Tubulin"]
    assert setting.value_text == "DNA, Tubulin"

    setting.value = ""                                   # and the empty case round-trips
    assert setting.value == []
