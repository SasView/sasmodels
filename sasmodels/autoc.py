"""
Automatically translate python models to C
"""
from __future__ import print_function

import ast
import inspect

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

def convert(info, module):
    # type: ("ModelInfo", ModuleType) -> bool
    """
    convert Iq, Iqxy and form_volume to c
    """
    # Check if there is already C code
    if info.source or info.c_code is not None:
        print("a", info.source, info.c_code)
        return

    public_methods = "Iq", "Iqxy", "form_volume"

    tagged = []
    translate = []
    for function_name in public_methods:
        function = getattr(info, function_name)
        if callable(function):
            translate.append(function)
            tagged.append(function_name)
    if not translate:
        print("b")
        return  # nothing to translate---maybe Iq, etc. are already C snippets?

    deps = []
    code = []
    while translate:
        function = translate.pop(0)
        source = inspect.getsource(function)
        tree = ast.parse(source)
        filename = function.__code__.co_filename
        offset = function.__code__.co_firstlineno
        refs = function.__code__.co_names
        snippet = codegen.to_source(tree) #, filename, offset)
        code.insert(0, snippet)
        for name in refs:
            obj = getattr(module, name, None)
            if obj is None:
                pass # ignore unbound variables for now
                #raise ValueError("global %s is not defined" % name)
            elif getattr(special, name, None):
                # special symbol: look up depenencies
                deps.extend(DEPENDENCY.get(name, []))
            elif callable(obj):
                if name not in tagged:
                    translate.append(obj)
                    tagged.append(name)
            elif isinstance(obj, float):
                code.insert(0, "const double %s = %.15g\n"%(name, obj))
            elif isinstance(obj, int):
                code.insert(0, "const int %s = %d\n"%(name, obj))
            elif isinstance(obj, (list, tuple, np.ndarray)):
                vals = ", ".join("%.15g"%v for v in obj)
                code.insert(0, "const double %s[] = {%s};\n" %(name, vals))
            elif isinstance(obj, special.Gauss):
                deps.append('lib/gauss%d.c'%obj.n)
            else:
                raise TypeError("Could not convert global %s of type %s"
                                % (name, str(type(obj))))

    # remove duplicates from the dependecy list
    unique_deps = []
    for dep in deps:
        if dep not in unique_deps:
            unique_deps.append(dep)

    info.source = unique_deps
    info.c_code = "\n".join(code)

    info.Iq = info.Iqxy = info.form_volume = None

    print("source", info.source)
    print(info.c_code)

    raise RuntimeError("not yet converted...")