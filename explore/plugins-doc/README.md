Intersphinx Exploration
=======================

An intersphinx playground for exploring the dynamic build of plugin docs help for saaview

Build using:

```shell
cd sasmodes/explore/plugins-doc
mkdir model
python -m sasmodels.genhelp ../../sasmodels/models/*ellipsoid*.py
python gen_plugins_toc.py
python -m sphinx -b html . html
```
