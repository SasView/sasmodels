"""
List all parameters used along with the models which use them.

Usage:

    python -m sasmodels.list_pars [-v]

If '-v' is given, then list the models containing the parameter in
addition to just the parameter name.
"""
import sys

from .core import load_model_definition
from .generate import make_info
from .compare import MODELS, columnize

def find_pars():
    partable = {}
    for name in sorted(MODELS):
        definition = load_model_definition(name)
        info = make_info(definition)
        for p in info['parameters']:
            pname = p[0]
            partable.setdefault(pname, [])
            partable[pname].append(name)
    return partable

def list_pars(names_only=True):
    partable = find_pars()
    if names_only:
        print(columnize(list(sorted(partable.keys()))))
    else:
        for k, v in sorted(partable.items()):
            print("%s: %s"%(k, ", ".join(v)))

def main():
    if len(sys.argv) == 2 and sys.argv[1] == '-v':
        verbose = True
    elif len(sys.argv) == 1:
        verbose = False
    else:
        print(__doc__)
        sys.exit(1)
    list_pars(names_only=not verbose)

if __name__ == "__main__":
    main()