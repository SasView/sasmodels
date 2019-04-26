"""
Core model handling routines.
"""
from __future__ import print_function

__all__ = [
    "list_models", "load_model", "load_model_info",
    "build_model", "precompile_dlls",
    ]

import os
from os.path import basename, join as joinpath
from glob import glob
import re

import numpy as np # type: ignore

from . import generate
from . import modelinfo
from . import product
from . import mixture
from . import kernelpy
from . import kernelcuda
from . import kernelcl
from . import kerneldll
from . import custom

# pylint: disable=unused-import
try:
    from typing import List, Union, Optional, Any
    from .kernel import KernelModel
    from .modelinfo import ModelInfo
except ImportError:
    pass
# pylint: enable=unused-import

CUSTOM_MODEL_PATH = os.environ.get('SAS_MODELPATH', "")
if CUSTOM_MODEL_PATH == "":
    CUSTOM_MODEL_PATH = joinpath(os.path.expanduser("~"), ".sasmodels", "custom_models")
    #if not os.path.isdir(CUSTOM_MODEL_PATH):
    #    os.makedirs(CUSTOM_MODEL_PATH)

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
        * c: c models only
        * single: c models which support single precision
        * double: c models which require double precision
        * opencl: c models which run in opencl
        * dll: c models which do not run in opencl
        * 1d: models without orientation
        * 2d: models with orientation
        * magnetic: models supporting magnetic sld
        * nommagnetic: models without magnetic parameter

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
    # TODO: may be adding Fq to the list at some point
    is_pure_py = callable(info.Iq)
    if kind == "py":
        return is_pure_py
    elif kind == "c":
        return not is_pure_py
    elif kind == "double":
        return not info.single and not is_pure_py
    elif kind == "single":
        return info.single and not is_pure_py
    elif kind == "opencl":
        return info.opencl
    elif kind == "dll":
        return not info.opencl and not is_pure_py
    elif kind == "2d":
        return any(p.type == 'orientation' for p in pars)
    elif kind == "1d":
        return all(p.type != 'orientation' for p in pars)
    elif kind == "magnetic":
        return any(p.type == 'sld' for p in pars)
    elif kind == "nonmagnetic":
        return not any(p.type == 'sld' for p in pars)
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
    if "+" in model_string:
        parts = [load_model_info(part)
                 for part in model_string.split("+")]
        return mixture.make_mixture_info(parts, operation='+')
    elif "*" in model_string:
        parts = [load_model_info(part)
                 for part in model_string.split("*")]
        return mixture.make_mixture_info(parts, operation='*')
    elif "@" in model_string:
        p_info, q_info = [load_model_info(part)
                          for part in model_string.split("@")]
        return product.make_product_info(p_info, q_info)
    # We are now dealing with a pure model
    elif "custom." in model_string:
        pattern = "custom.([A-Za-z0-9_-]+)"
        result = re.match(pattern, model_string)
        if result is None:
            raise ValueError("Model name in invalid format: " + model_string)
        model_name = result.group(1)
        # Use ModelName to find the path to the custom model file
        model_path = joinpath(CUSTOM_MODEL_PATH, model_name + ".py")
        if not os.path.isfile(model_path):
            raise ValueError("The model file {} doesn't exist".format(model_path))
        kernel_module = custom.load_custom_kernel_module(model_path)
        return modelinfo.make_model_info(kernel_module)
    kernel_module = generate.load_kernel_module(model_string)
    return modelinfo.make_model_info(kernel_module)


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
    elif platform == "cuda":
        return kernelcuda.GpuModel(source, model_info, numpy_dtype, fast=fast)
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
            old_path = kerneldll.SAS_DLL_PATH
            try:
                kerneldll.SAS_DLL_PATH = path
                dll = kerneldll.make_dll(source, model_info, dtype=numpy_dtype)
            finally:
                kerneldll.SAS_DLL_PATH = old_path
            compiled_dlls.append(dll)
    return compiled_dlls

def parse_dtype(model_info, dtype=None, platform=None):
    # type: (ModelInfo, str, str) -> (np.dtype, bool, str)
    """
    Interpret dtype string, returning np.dtype, fast flag and platform.

    Possible types include 'half', 'single', 'double' and 'quad'.  If the
    type is 'fast', then this is equivalent to dtype 'single' but using
    fast native functions rather than those with the precision level
    guaranteed by the OpenCL standard.  'default' will choose the appropriate
    default for the model and platform.

    Platform preference can be specfied ("ocl", "cuda", "dll"), with the
    default being OpenCL or CUDA if available, otherwise DLL.  If the dtype
    name ends with '!' then platform is forced to be DLL rather than GPU.
    The default platform is set by the environment variable SAS_OPENCL,
    SAS_OPENCL=driver:device for OpenCL, SAS_OPENCL=cuda:device for CUDA
    or SAS_OPENCL=none for DLL.

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

    # Check if type indicates dll regardless of which platform is given
    if dtype is not None and dtype.endswith('!'):
        platform = "dll"
        dtype = dtype[:-1]

    # Make sure model allows opencl/gpu
    if not model_info.opencl:
        platform = "dll"

    # Make sure opencl is available, or fallback to cuda then to dll
    if platform == "ocl" and not kernelcl.use_opencl():
        platform = "cuda" if kernelcuda.use_cuda() else "dll"

    # Convert special type names "half", "fast", and "quad"
    fast = (dtype == "fast")
    if fast:
        dtype = "single"
    elif dtype == "quad":
        dtype = "longdouble"
    elif dtype == "half":
        dtype = "float16"

    # Convert dtype string to numpy dtype.  Use single precision for GPU
    # if model allows it, otherwise use double precision.
    if dtype is None or dtype == "default":
        numpy_dtype = (generate.F32 if model_info.single and platform in ("ocl", "cuda")
                       else generate.F64)
    else:
        numpy_dtype = np.dtype(dtype)

    # Make sure that the type is supported by GPU, otherwise use dll
    if platform == "ocl":
        env = kernelcl.environment()
    elif platform == "cuda":
        env = kernelcuda.environment()
    else:
        env = None
    if env is not None and not env.has_type(numpy_dtype):
        platform = "dll"
        if dtype is None:
            numpy_dtype = generate.F64

    return numpy_dtype, fast, platform

def test_composite_order():
    """
    Check that mixture models produce the same result independent of ordder.
    """
    def test_models(fst, snd):
        """Confirm that two models produce the same parameters"""
        fst = load_model(fst)
        snd = load_model(snd)
        # Un-disambiguate parameter names so that we can check if the same
        # parameters are in a pair of composite models. Since each parameter in
        # the mixture model is tagged as e.g., A_sld, we ought to use a
        # regex subsitution s/^[A-Z]+_/_/, but removing all uppercase letters
        # is good enough.
        fst = [[x for x in p.name if x == x.lower()] for p in fst.info.parameters.kernel_parameters]
        snd = [[x for x in p.name if x == x.lower()] for p in snd.info.parameters.kernel_parameters]
        assert sorted(fst) == sorted(snd), "{} != {}".format(fst, snd)

    test_models(
        "cylinder+sphere",
        "sphere+cylinder")
    test_models(
        "cylinder*sphere",
        "sphere*cylinder")
    test_models(
        "cylinder@hardsphere*sphere",
        "sphere*cylinder@hardsphere")
    test_models(
        "barbell+sphere*cylinder@hardsphere",
        "sphere*cylinder@hardsphere+barbell")
    test_models(
        "barbell+cylinder@hardsphere*sphere",
        "cylinder@hardsphere*sphere+barbell")
    test_models(
        "barbell+sphere*cylinder@hardsphere",
        "barbell+cylinder@hardsphere*sphere")
    test_models(
        "sphere*cylinder@hardsphere+barbell",
        "cylinder@hardsphere*sphere+barbell")
    test_models(
        "barbell+sphere*cylinder@hardsphere",
        "cylinder@hardsphere*sphere+barbell")
    test_models(
        "barbell+cylinder@hardsphere*sphere",
        "sphere*cylinder@hardsphere+barbell")

def test_composite():
    # type: () -> None
    """Check that model load works"""
    from .product import RADIUS_ID, VOLFRAC_ID, STRUCTURE_MODE_ID, RADIUS_MODE_ID
    #Test the the model produces the parameters that we would expect
    model = load_model("cylinder@hardsphere*sphere")
    actual = [p.name for p in model.info.parameters.kernel_parameters]
    target = ["sld", "sld_solvent", "radius", "length", "theta", "phi",
              RADIUS_ID, VOLFRAC_ID, STRUCTURE_MODE_ID, RADIUS_MODE_ID,
              "A_sld", "A_sld_solvent", "A_radius"]
    assert target == actual, "%s != %s"%(target, actual)

def list_models_main():
    # type: () -> int
    """
    Run list_models as a main program.  See :func:`list_models` for the
    kinds of models that can be requested on the command line.
    """
    import sys
    kind = sys.argv[1] if len(sys.argv) > 1 else "all"
    try:
        models = list_models(kind)
        print("\n".join(models))
    except Exception:
        print(list_models.__doc__)
        return 1
    return 0

if __name__ == "__main__":
    list_models_main()
