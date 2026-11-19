"""Custom model file paths and model expressions share the loading API."""

import shutil
from pathlib import Path

import pytest

from sasmodels import core
from sasmodels.models import sphere


def write_sphere(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(sphere.__file__, path)
    return path


@pytest.mark.parametrize("directory", ["lab+x", "me@org", "my.custom.models"])
def test_load_model_file(tmp_path, directory):
    path = write_sphere(tmp_path / directory / "mysphere.py")

    info = core.load_model_info(str(path))

    assert Path(info.filename) == path
    assert info.name == "sphere"
    assert info.parameters.defaults == core.load_model_info("sphere").parameters.defaults


def test_expression_with_model_file(tmp_path):
    path = write_sphere(tmp_path / "mysphere.py")
    info = core.load_model_info("sphere+" + str(path))
    assert info.composition[0] == "mixture"
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
