"""Custom model file paths and model expressions share the loading API."""

import os
import shutil
from pathlib import Path

import pytest

from sasmodels import core
from sasmodels.models import sphere


def write_sphere(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(sphere.__file__, path)
    return path


@pytest.mark.parametrize("name", [
    "plain", "lab+x", "me@org", "my.custom.models",
    pytest.param("lab*x", marks=pytest.mark.skipif(
        os.name == "nt", reason="Windows does not allow '*' in filenames")),
])
@pytest.mark.parametrize("in_filename", [False, True])
@pytest.mark.parametrize("relative", [False, True])
def test_load_model_file(tmp_path, monkeypatch, name, in_filename, relative):
    path = tmp_path / (name + ".py") if in_filename else tmp_path / name / "mysphere.py"
    write_sphere(path)
    monkeypatch.chdir(tmp_path)
    argument = str(path.relative_to(tmp_path) if relative else path)

    info = core.load_model_info(argument)

    assert Path(info.filename) == path
    assert info.name == "sphere"
    assert info.parameters.defaults == core.load_model_info("sphere").parameters.defaults


def test_load_missing_model_file(tmp_path):
    path = tmp_path / "missing.py"
    with pytest.raises(FileNotFoundError) as error:
        core.load_model_info(str(path))
    assert error.value.filename == str(path)


@pytest.mark.parametrize("expression, operation", [
    ("sphere+{path}", "+"), ("{path}+sphere", "+"),
    ("sphere*{path}", "*"), ("{path}*sphere", "*"),
    ("{path}@hardsphere", "@"),
])
def test_expression_with_model_file(tmp_path, expression, operation):
    path = write_sphere(tmp_path / "mysphere.py")
    info = core.load_model_info(expression.format(path=path))
    assert info.composition[0] == ("product" if operation == "@" else "mixture")
    assert len(info.composition[1]) == 2
    assert any(part.filename == str(path) for part in info.composition[1])


def test_sasview_sum_plugin_file(tmp_path):
    path = tmp_path / "lab+x" / "sum_model.py"
    path.parent.mkdir()
    path.write_text(
        "from sasmodels.core import load_model_info\n"
        "from sasmodels.sasview_model import make_model_from_info\n"
        "model_info = load_model_info('sphere+cylinder')\n"
        "Model = make_model_from_info(model_info)\n",
        encoding="utf-8",
    )

    info = core.load_model_info(str(path))

    assert info.composition[0] == "mixture"
    assert [part.id for part in info.composition[1]] == ["sphere", "cylinder"]


def test_user_relative_model_file(tmp_path, monkeypatch):
    path = write_sphere(tmp_path / "lab+x" / "mysphere.py")
    monkeypatch.setenv("USERPROFILE" if os.name == "nt" else "HOME", str(tmp_path))
    info = core.load_model_info(os.path.join("~", "lab+x", "mysphere.py"))
    assert Path(info.filename) == path


def test_custom_model_prefix(tmp_path, monkeypatch):
    path = write_sphere(tmp_path / "mysphere.py")
    monkeypatch.setattr(core, "CUSTOM_MODEL_PATH", str(tmp_path))
    info = core.load_model_info("custom.mysphere")
    assert Path(info.filename) == path
