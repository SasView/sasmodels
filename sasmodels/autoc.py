"""
Automatically translate python models to C
"""
from __future__ import print_function

import ast
import inspect
from functools import reduce

import numpy as np

from . import codegen
from . import special

# pylint: disable=unused-import
try:
    from types import ModuleType
    #from .modelinfo import ModelInfo  # circular import
except ImportError:
    pass
# pylint: enable=unused-import

DEPENDENCY = {
    'core_shell_kernel': ['lib/core_shell.c'],
    'fractal_sq': ['lib/fractal_sq.c'],
    'gfn4': ['lib/gfn.c'],
    'polevl': ['lib/polevl.c'],
    'p1evl': ['lib/polevl.c'],
    'sas_2J1x_x': ['lib/polevl.c', 'lib/sas_J1.c'],
    'sas_3j1x_x': ['lib/sas_3j1x_x.c'],
    'sas_erf': ['lib/polevl.c', 'lib/sas_erf.c'],
    'sas_erfc': ['lib/polevl.c', 'lib/sas_erf.c'],
    'sas_gamma': ['lib/sas_gamma.c'],
    'sas_J0': ['lib/polevl.c', 'lib/sas_J0.c'],
    'sas_J1': ['lib/polevl.c', 'lib/sas_J1.c'],
    'sas_JN': ['lib/polevl.c', 'lib/sas_J0.c', 'lib/sas_J1.c', 'lib/sas_JN.c'],
    'sas_Si': ['lib/Si.c'],
}

DEFINES = frozenset("M_PI M_PI_2 M_PI_4 M_SQRT1_2 M_E NAN INFINITY M_PI_180 M_4PI_3".split())

def convert(info, module):
    # type: ("ModelInfo", ModuleType) -> bool
    """
    convert Iq, Iqxy and form_volume to c
    """
    # Check if there is already C code
    if info.source or info.c_code is not None:
        return

    public_methods = "Iq", "Iqxy", "form_volume"

    tagged = [] # type: List[str]
    translate = [] # type: List[Callable]
    for function_name in public_methods:
        function = getattr(info, function_name)
        if callable(function):
            if getattr(function, 'vectorized', None):
                return  # Don't try to translate vectorized code
            tagged.append(function_name)
            translate.append((function_name, function))
    if not translate:
        # nothing to translate---maybe Iq, etc. are already C snippets?
        return

    libs = []  # type: List[str]
    code = {}  # type: Dict[str, str]
    depends = {}  # type: Dict[str, List[str]]
    while translate:
        function_name, function = translate.pop(0)
        filename = function.__code__.co_filename
        offset = function.__code__.co_firstlineno
        refs = function.__code__.co_names
        depends[function_name] = set(refs)
        source = inspect.getsource(function)
        for name in refs:
            if name in tagged or name in DEFINES:
                continue
            tagged.append(name)
            obj = getattr(module, name, None)
            if obj is None:
                pass # ignore unbound variables for now
                #raise ValueError("global %s is not defined" % name)
            elif callable(obj):
                if getattr(special, name, None):
                    # special symbol: look up depenencies
                    libs.extend(DEPENDENCY.get(name, []))
                else:
                    # not special: add function to translate stack
                    translate.append((name, obj))
            elif isinstance(obj, float):
                code[name] = "const double %s = %.15g\n"%(name, obj)
            elif isinstance(obj, int):
                code[name] = "const int %s = %d\n"%(name, obj)
            elif isinstance(obj, (list, tuple, np.ndarray)):
                vals = ", ".join("%.15g"%v for v in obj)
                code[name] = "const double %s[] = {%s};\n" %(name, vals)
            elif isinstance(obj, special.Gauss):
                libs.append('lib/gauss%d.c'%obj.n)
                source = (source.replace(name+'.n', 'GAUSS_N')
                          .replace(name+'.z', 'GAUSS_Z')
                          .replace(name+'.w', 'GAUSS_W'))
            else:
                raise TypeError("Could not convert global %s of type %s"
                                % (name, str(type(obj))))

        tree = ast.parse(source)
        snippet = codegen.to_source(tree) #, filename, offset)
        code[function_name] = snippet

    # remove duplicates from the dependecy list
    unique_libs = []
    for filename in libs:
        if filename not in unique_libs:
            unique_libs.append(filename)

    info.source = unique_libs
    info.c_code = "\n".join(code[k] for k in ordered_dag(depends) if k in code)

    info.Iq = info.Iqxy = info.form_volume = None

    print("source", info.source)
    print(info.c_code)

    raise RuntimeError("not yet converted...")


# Modified from the following:
#
#    http://code.activestate.com/recipes/578272-topological-sort/
#    Copyright (C) 2012 Sam Denton
#    License: MIT
def ordered_dag(dag):
    # type: (Dict[T, Set[T]]) -> Iterator[T]
    dag = dag.copy()

    # make leaves depend on the empty set
    leaves = reduce(set.union, dag.values()) - set(dag.keys())
    dag.update({node: set() for node in leaves})
    while True:
        leaves = set(node for node, links in dag.items() if not links)
        if not leaves:
            break
        for node in leaves:
            yield node
        dag = {node: (links-leaves)
               for node, links in dag.items() if node not in leaves}
    if dag:
        raise ValueError("Cyclic dependes exists amongst these items:\n%s"
                            % ", ".join(str(node) for node in dag.keys()))
