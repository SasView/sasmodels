#!/usr/bin/env python3
"""Create plugins.rst with an alphabetically sorted toctree of RST files."""

from pathlib import Path

directory = Path(__file__).resolve().parent / "model"
output = directory / "plugins.rst"

rst_files = sorted(
    (
        path
        for path in directory.glob("*.rst")
        if path.name != output.name
    ),
    key=lambda path: path.stem.casefold(),
)

contents = """.. _model-plugins:

*******
Plugins
*******

.. toctree::

"""

contents += "".join(f"    {path.name}\n" for path in rst_files)
output.write_text(contents, encoding="utf-8")
