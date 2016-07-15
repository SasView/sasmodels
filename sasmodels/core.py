"""
Core model handling routines.
"""
from __future__ import print_function

__all__ = [
    "list_models", "load_model", "load_model_info",
    "build_model", "precompile_dll",
    ]

import os
from os.path import basename, dirname, join as joinpath, splitext
from glob import glob

import numpy as np # type: ignore

from . import generate
from . import modelinfo
from . import product
from . import mixture
from . import kernelpy
from . import kerneldll
try:
    from . import kernelcl
    HAVE_OPENCL = True
except Exception:
    HAVE_OPENCL = False

try:
    from typing import List, Union, Optional, Any
    from .kernel import KernelModel
    from .modelinfo import ModelInfo
except ImportError:
    pass

try:
    np.meshgrid([])
    meshgrid = np.meshgrid
except Exception:
    # CRUFT: np.meshgrid requires multiple vectors
    def meshgrid(*args):
        if len(args) > 1:
            return np.meshgrid(*args)
        else:
            return [np.asarray(v) for v in args]

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

def list_models():
    # type: () -> List[str]
    """
    Return the list of available models on the model path.
    """
    root = dirname(__file__)
    files = sorted(glob(joinpath(root, 'models', "[a-zA-Z]*.py")))
    available_models = [basename(f)[:-3] for f in files]
    return available_models

def load_model(model_name, dtype=None, platform='ocl'):
    # type: (str, str, str) -> KernelModel
    """
    Load model info and build model.

    *model_name* is the name of the model as used by :func:`load_model_info`.
    Additional keyword arguments are passed directly to :func:`build_model`.
    """
    return build_model(load_model_info(model_name),
                       dtype=dtype, platform=platform)


def load_model_info(model_name):
    # type: (str) -> modelinfo.ModelInfo
    """
    Load a model definition given the model name.

    This returns a handle to the module defining the model.  This can be
    used with functions in generate to build the docs or extract model info.
    """
    parts = model_name.split('+')
    if len(parts) > 1:
        model_info_list = [load_model_info(p) for p in parts]
        return mixture.make_mixture_info(model_info_list)

    parts = model_name.split('*')
    if len(parts) > 1:
        if len(parts) > 2:
            raise ValueError("use P*S to apply structure factor S to model P")
        P_info, Q_info = [load_model_info(p) for p in parts]
        return product.make_product_info(P_info, Q_info)

    kernel_module = generate.load_kernel_module(model_name)
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
            from . import product
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
        return kerneldll.load_dll(source, model_info, numpy_dtype)
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
        source = generate.make_source(model_info)
        if source:
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
    type is 'fast', then this is equivalent to dtype 'single' with the
    fast flag set to True.
    """
    # Assign default platform, overriding ocl with dll if OpenCL is unavailable
    if platform is None:
        platform = "ocl"
    if platform=="ocl" and not HAVE_OPENCL:
        platform = "dll"

    # Check if type indicates dll regardless of which platform is given
    if dtype is not None and dtype.endswith('!'):
        platform = "dll"
        dtype = dtype[:-1]

    # Convert special type names "half", "fast", and "quad"
    fast = (dtype=="fast")
    if fast:
        dtype = "single"
    elif dtype=="quad":
        dtype = "longdouble"
    elif dtype=="half":
        dtype = "f16"

    # Convert dtype string to numpy dtype.
    if dtype is None:
        numpy_dtype = generate.F32 if platform=="ocl" and model_info.single else generate.F64
    else:
        numpy_dtype = np.dtype(dtype)

    # Make sure that the type is supported by opencl, otherwise use dll
    if platform=="ocl":
        env = kernelcl.environment()
        if not env.has_type(numpy_dtype):
            platform = "dll"
            if dtype is None:
                numpy_dtype = generate.F64

    return numpy_dtype, fast, platform
