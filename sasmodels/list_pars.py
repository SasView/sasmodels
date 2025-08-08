"""
List all parameters used along with the models which use them.

Usage:

    python -m sasmodels.list_pars [-v]

If '-v' is given, then list the models containing the parameter in
addition to just the parameter name.
"""

import argparse

from .compare import columnize
from .core import list_models, load_model_info


def find_pars(kind=None):
    """
    Find all parameters in all models.

    Returns the reference table *{parameter: [model, model, ...]}*
    """
    partable = {}
    for name in list_models():
        model_info = load_model_info(name)
        for p in model_info.parameters.kernel_parameters:
            if kind is None or p.type == kind:
                partable.setdefault(p.name, [])
                partable[p.name].append(name)
    return partable

def list_pars(names_only=True, kind=None):
    """
    Print all parameters in all models.

    If *names_only* then only print the parameter name, not the models it
    occurs in.
    """
    partable = find_pars(kind)
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
        'kind', default="any", nargs='?',
        metavar="volume|orientation|sld|none|any",
        choices=['volume', 'orientation', 'sld', None, 'any'],
        type=lambda v: None if v == 'any' else None if v == 'none' else v,
        help="only list arguments of the given kind")
    args = parser.parse_args()

    list_pars(names_only=not args.verbose, kind=args.kind)

if __name__ == "__main__":
    main()
