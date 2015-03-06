import sys
# make sure sasmodels is on the path
sys.path.append('..')

from os import mkdir
from os.path import basename, exists, join as joinpath
from sasmodels.core import load_model_definition


TEMPLATE="""\
..
    Generated from doc/gentoc.py -- DO NOT EDIT --

.. %(label)s:

%(bar)s
%(title)s
%(bar)s

.. toctree::

"""

MODEL_TOC_PATH = "ref/models"

def _make_category(category_name, label, title, parent=None):
    file = open(joinpath(MODEL_TOC_PATH, category_name+".rst"), "w")
    file.write(TEMPLATE%{'label':label, 'title':title, 'bar':'*'*len(title)})
    if parent:
        _add_subcategory(parent, category_name)
    return file

def _add_subcategory(file, category_name):
    file.write("    %s.rst\n"%category_name)

def _add_model(file, model_name):
    file.write("    ../../model/%s.rst\n"%model_name)

def generate_toc(model_files):
    if not model_files:
        print >>sys.stderr, "gentoc needs a list of model files"

    # find all categories
    category = {}
    for item in model_files:
        # assume model is in sasmodels/models/name.py, and ignore the full path
        model_name = basename(item)[:-3]
        if model_name.startswith('_'): continue
        model_definition = load_model_definition(model_name)
        if not hasattr(model_definition, 'category'):
            print >>sys.stderr, "Missing category for",item
        else:
            category.setdefault(model_definition.category,[]).append(model_name)

    # Check category names
    for k,v in category.items():
        if len(v) == 1:
            print >>sys.stderr, "Category %s contains only %s"%(k,v[0])

    # Generate category files
    # Assume that top-level groups are:
    #    shape, shape-independent, structure-factor and custom
    # Supports a two-level category structure.

    if not exists(MODEL_TOC_PATH): mkdir(MODEL_TOC_PATH)
    model_toc = _make_category(
        'index',  'Models', 'Model Functions')
    shape_toc = _make_category(
        'shape',  'Shapes', 'Shape Functions', model_toc)
    free_toc = _make_category(
        'shape-independent',  'Shape-independent',
        'Shape-Independent Functions', model_toc)
    struct_toc = _make_category(
        'structure-factor',  'Structure-factor', 'Structure Factors', model_toc)
    custom_toc = _make_category(
        'custom-models',  'Custom-models', 'Custom Models', model_toc)
    model_toc.close()

    # remember to top level categories
    cat_files = {
        'shape':shape_toc,
        'shape-independent':free_toc,
        'structure-factor': struct_toc,
        'custom': custom_toc,
        }

    # Process the model lists
    for k,v in sorted(category.items()):
        if ':' in k:
            cat,subcat = k.split(':')
            cat_file = cat_files[cat]
            label = "-".join((cat,subcat))
            filename = label
            title = subcat.capitalize()+" Functions"
            sub_toc = _make_category(filename, label, title, cat_file)
            for model in sorted(v):
                _add_model(sub_toc, model)
            sub_toc.close()
        else:
            if k not in cat_files:
                print >>sys.stderr, "Unknown category %s containing"%cat, v
            else:
                cat_file = cat_files[k]
                for model in sorted(v):
                    _add_model(cat_file, model)

    # Close the top-level category files
    for f in cat_files.values(): f.close()


if __name__ == "__main__":
    generate_toc(sys.argv[1:])