"""
List all parameters used along with the models which use them.

Usage:

    python -m sasmodels.list_pars [-v]

If '-v' is given, then list the models containing the parameter in
addition to just the parameter name.
"""
from __future__ import print_function

import sys

from .core import load_model_info
from .compare import MODELS, columnize

def find_pars():
    """
    Find all parameters in all models.

    Returns the reference table *{parameter: [model, model, ...]}*
    """
    partable = {}
    for name in sorted(MODELS):
        model_info = load_model_info(name)
        for p in model_info['parameters']:
            partable.setdefault(p.name, [])
            partable[p.name].append(name)
    return partable

def list_pars(names_only=True):
    """
    Print all parameters in all models.

    If *names_only* then only print the parameter name, not the models it
    occurs in.
    """
    partable = find_pars()
    if names_only:
        print(columnize(list(sorted(partable.keys()))))
    else:
        for k, v in sorted(partable.items()):
            print("%s: %s"%(k, ", ".join(v)))

def main():
    """
    Program to list the parameters used across all models.
    """
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
