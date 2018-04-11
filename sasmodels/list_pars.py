"""
List all parameters used along with the models which use them.

Usage:

    python -m sasmodels.list_pars [-v]

If '-v' is given, then list the models containing the parameter in
addition to just the parameter name.
"""
from __future__ import print_function

import argparse

from .core import load_model_info, list_models
from .compare import columnize

def find_pars(type=None):
    """
    Find all parameters in all models.

    Returns the reference table *{parameter: [model, model, ...]}*
    """
    partable = {}
    for name in list_models():
        model_info = load_model_info(name)
        for p in model_info.parameters.kernel_parameters:
            if type is None or p.type == type:
                partable.setdefault(p.name, [])
                partable[p.name].append(name)
    return partable

def list_pars(names_only=True, type=None):
    """
    Print all parameters in all models.

    If *names_only* then only print the parameter name, not the models it
    occurs in.
    """
    partable = find_pars(type)
    if names_only:
        print(columnize(list(sorted(partable.keys()))))
    else:
        for k, v in sorted(partable.items()):
            print("%s: %s"%(k, ", ".join(v)))

def main():
    """
    Program to list the parameters used across all models.
    """
    parser = argparse.ArgumentParser(
        description="Find all parameters in all models",
        )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help="list models which use this argument")
    parser.add_argument(
        'type', default="any", nargs='?',
        metavar="volume|orientation|sld|none|any",
        choices=['volume', 'orientation', 'sld', None, 'any'],
        type=lambda v: None if v == 'any' else '' if v == 'none' else v,
        help="only list arguments of the given type")
    args = parser.parse_args()

    list_pars(names_only=not args.verbose, type=args.type)

if __name__ == "__main__":
    main()
