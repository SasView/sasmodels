"""
Core model handling routines.
"""
from __future__ import print_function

__all__ = [
    "list_models", "load_model", "load_model_info",
    "build_model", "precompile_dlls",
    ]

import os
import re
from os.path import basename, dirname, join as joinpath
from glob import glob

import numpy as np # type: ignore

from . import generate
from . import modelinfo
from . import product
from . import mixture
from . import kernelpy
from . import kerneldll
from . import custom

if os.environ.get("SAS_OPENCL", "").lower() == "none":
    HAVE_OPENCL = False
else:
    try:
        from . import kernelcl
        HAVE_OPENCL = True
    except Exception:
        HAVE_OPENCL = False

CUSTOM_MODEL_PATH = os.environ.get('SAS_MODELPATH', "")
if CUSTOM_MODEL_PATH == "":
    path = joinpath(os.path.expanduser("~"), ".sasmodels", "custom_models")
    if not os.path.isdir(path):
        os.makedirs(path)
    CUSTOM_MODEL_PATH = path

try:
    from typing import List, Union, Optional, Any
    from .kernel import KernelModel
    from .modelinfo import ModelInfo
except ImportError:
    pass

# TODO: refactor composite model support
# The current load_model_info/build_model does not reuse existing model
# definitions when loading a composite model, instead reloading and
# rebuilding the kernel for each component model in the expression.  This
# is fine in a scripting environment where the model is built when the script
# starts and is thrown away when the script ends, but may not be the best
# solution in a long-lived application.  This affects the following functions:
#
#    load_model
#    load_model_info
#    build_model

KINDS = ("all", "py", "c", "double", "single", "opencl", "1d", "2d",
         "nonmagnetic", "magnetic")
def list_models(kind=None):
    # type: (str) -> List[str]
    """
    Return the list of available models on the model path.

    *kind* can be one of the following:

        * all: all models
        * py: python models only
        * c: compiled models only
        * single: models which support single precision
        * double: models which require double precision
        * opencl: controls if OpenCL is supperessed
        * 1d: models which are 1D only, or 2D using abs(q)
        * 2d: models which can be 2D
        * magnetic: models with an sld
        * nommagnetic: models without an sld

    For multiple conditions, combine with plus.  For example, *c+single+2d*
    would return all oriented models implemented in C which can be computed
    accurately with single precision arithmetic.
    """
    if kind and any(k not in KINDS for k in kind.split('+')):
        raise ValueError("kind not in " + ", ".join(KINDS))
    files = sorted(glob(joinpath(generate.MODEL_PATH, "[a-zA-Z]*.py")))
    available_models = [basename(f)[:-3] for f in files]
    if kind and '+' in kind:
        all_kinds = kind.split('+')
        condition = lambda name: all(_matches(name, k) for k in all_kinds)
    else:
        condition = lambda name: _matches(name, kind)
    selected = [name for name in available_models if condition(name)]

    return selected

def _matches(name, kind):
    if kind is None or kind == "all":
        return True
    info = load_model_info(name)
    pars = info.parameters.kernel_parameters
    if kind == "py" and callable(info.Iq):
        return True
    elif kind == "c" and not callable(info.Iq):
        return True
    elif kind == "double" and not info.single:
        return True
    elif kind == "single" and info.single:
        return True
    elif kind == "opencl" and info.opencl:
        return True
    elif kind == "2d" and any(p.type == 'orientation' for p in pars):
        return True
    elif kind == "1d" and all(p.type != 'orientation' for p in pars):
        return True
    elif kind == "magnetic" and any(p.type == 'sld' for p in pars):
        return True
    elif kind == "nonmagnetic" and any(p.type != 'sld' for p in pars):
        return True
    return False

def load_model(model_name, dtype=None, platform='ocl'):
    # type: (str, str, str) -> KernelModel
    """
    Load model info and build model.

    *model_name* is the name of the model, or perhaps a model expression
    such as sphere*hardsphere or sphere+cylinder.

    *dtype* and *platform* are given by :func:`build_model`.
    """
    return build_model(load_model_info(model_name),
                       dtype=dtype, platform=platform)

def load_model_info(model_string):
    # type: (str) -> modelinfo.ModelInfo
    """
    Load a model definition given the model name.

    *model_string* is the name of the model, or perhaps a model expression
    such as sphere*cylinder or sphere+cylinder. Use '@' for a structure
    factor product, e.g. sphere@hardsphere. Custom models can be specified by
    prefixing the model name with 'custom.', e.g. 'custom.MyModel+sphere'.

    This returns a handle to the module defining the model.  This can be
    used with functions in generate to build the docs or extract model info.
    """
    if '@' in model_string:
        parts = model_string.split('@')
        if len(parts) != 2:
            raise ValueError("Use P@S to apply a structure factor S to model P")
        P_info, Q_info = [load_model_info(part) for part in parts]
        return product.make_product_info(P_info, Q_info)

    product_parts = []
    addition_parts = []

    addition_parts_names = model_string.split('+')
    if len(addition_parts_names) >= 2:
        addition_parts = [load_model_info(part) for part in addition_parts_names]
    elif len(addition_parts_names) == 1:
        product_parts_names = model_string.split('*')
        if len(product_parts_names) >= 2:
            product_parts = [load_model_info(part) for part in product_parts_names]
        elif len(product_parts_names) == 1:
            if "custom." in product_parts_names[0]:
                # Extract ModelName from "custom.ModelName"
                pattern = "custom.([A-Za-z0-9_-]+)"
                result = re.match(pattern, product_parts_names[0])
                if result is None:
                    raise ValueError("Model name in invalid format: " + product_parts_names[0])
                model_name = result.group(1)
                # Use ModelName to find the path to the custom model file
                model_path = joinpath(CUSTOM_MODEL_PATH, model_name + ".py")
                if not os.path.isfile(model_path):
                    raise ValueError("The model file {} doesn't exist".format(model_path))
                kernel_module = custom.load_custom_kernel_module(model_path)
                return modelinfo.make_model_info(kernel_module)
            # Model is a core model
            kernel_module = generate.load_kernel_module(product_parts_names[0])
            return modelinfo.make_model_info(kernel_module)

    model = None
    if len(product_parts) > 1:
        model = mixture.make_mixture_info(product_parts, operation='*')
    if len(addition_parts) > 1:
        if model is not None:
            addition_parts.append(model)
        model = mixture.make_mixture_info(addition_parts, operation='+')
    return model


def build_model(model_info, dtype=None, platform="ocl"):
    # type: (modelinfo.ModelInfo, str, str) -> KernelModel
    """
    Prepare the model for the default execution platform.

    This will return an OpenCL model, a DLL model or a python model depending
    on the model and the computing platform.

    *model_info* is the model definition structure returned from
    :func:`load_model_info`.

    *dtype* indicates whether the model should use single or double precision
    for the calculation.  Choices are 'single', 'double', 'quad', 'half',
    or 'fast'.  If *dtype* ends with '!', then force the use of the DLL rather
    than OpenCL for the calculation.

    *platform* should be "dll" to force the dll to be used for C models,
    otherwise it uses the default "ocl".
    """
    composition = model_info.composition
    if composition is not None:
        composition_type, parts = composition
        models = [build_model(p, dtype=dtype, platform=platform) for p in parts]
        if composition_type == 'mixture':
            return mixture.MixtureModel(model_info, models)
        elif composition_type == 'product':
            P, S = models
            return product.ProductModel(model_info, P, S)
        else:
            raise ValueError('unknown mixture type %s'%composition_type)

    # If it is a python model, return it immediately
    if callable(model_info.Iq):
        return kernelpy.PyModel(model_info)

    numpy_dtype, fast, platform = parse_dtype(model_info, dtype, platform)

    source = generate.make_source(model_info)
    if platform == "dll":
        #print("building dll", numpy_dtype)
        return kerneldll.load_dll(source['dll'], model_info, numpy_dtype)
    else:
        #print("building ocl", numpy_dtype)
        return kernelcl.GpuModel(source, model_info, numpy_dtype, fast=fast)

def precompile_dlls(path, dtype="double"):
    # type: (str, str) -> List[str]
    """
    Precompile the dlls for all builtin models, returning a list of dll paths.

    *path* is the directory in which to save the dlls.  It will be created if
    it does not already exist.

    This can be used when build the windows distribution of sasmodels
    which may be missing the OpenCL driver and the dll compiler.
    """
    numpy_dtype = np.dtype(dtype)
    if not os.path.exists(path):
        os.makedirs(path)
    compiled_dlls = []
    for model_name in list_models():
        model_info = load_model_info(model_name)
        if not callable(model_info.Iq):
            source = generate.make_source(model_info)['dll']
            old_path = kerneldll.DLL_PATH
            try:
                kerneldll.DLL_PATH = path
                dll = kerneldll.make_dll(source, model_info, dtype=numpy_dtype)
            finally:
                kerneldll.DLL_PATH = old_path
            compiled_dlls.append(dll)
    return compiled_dlls

def parse_dtype(model_info, dtype=None, platform=None):
    # type: (ModelInfo, str, str) -> (np.dtype, bool, str)
    """
    Interpret dtype string, returning np.dtype and fast flag.

    Possible types include 'half', 'single', 'double' and 'quad'.  If the
    type is 'fast', then this is equivalent to dtype 'single' but using
    fast native functions rather than those with the precision level
    guaranteed by the OpenCL standard.  'default' will choose the appropriate
    default for the model and platform.

    Platform preference can be specfied ("ocl" vs "dll"), with the default
    being OpenCL if it is availabe.  If the dtype name ends with '!' then
    platform is forced to be DLL rather than OpenCL.

    This routine ignores the preferences within the model definition.  This
    is by design.  It allows us to test models in single precision even when
    we have flagged them as requiring double precision so we can easily check
    the performance on different platforms without having to change the model
    definition.
    """
    # Assign default platform, overriding ocl with dll if OpenCL is unavailable
    # If opencl=False OpenCL is switched off

    if platform is None:
        platform = "ocl"
    if platform == "ocl" and not HAVE_OPENCL or not model_info.opencl:
        platform = "dll"

    # Check if type indicates dll regardless of which platform is given
    if dtype is not None and dtype.endswith('!'):
        platform = "dll"
        dtype = dtype[:-1]

    # Convert special type names "half", "fast", and "quad"
    fast = (dtype == "fast")
    if fast:
        dtype = "single"
    elif dtype == "quad":
        dtype = "longdouble"
    elif dtype == "half":
        dtype = "float16"

    # Convert dtype string to numpy dtype.
    if dtype is None or dtype == "default":
        numpy_dtype = (generate.F32 if platform == "ocl" and model_info.single
                       else generate.F64)
    else:
        numpy_dtype = np.dtype(dtype)

    # Make sure that the type is supported by opencl, otherwise use dll
    if platform == "ocl":
        env = kernelcl.environment()
        if not env.has_type(numpy_dtype):
            platform = "dll"
            if dtype is None:
                numpy_dtype = generate.F64

    return numpy_dtype, fast, platform

def list_models_main():
    # type: () -> None
    """
    Run list_models as a main program.  See :func:`list_models` for the
    kinds of models that can be requested on the command line.
    """
    import sys
    kind = sys.argv[1] if len(sys.argv) > 1 else "all"
    print("\n".join(list_models(kind)))

if __name__ == "__main__":
    list_models_main()
