"""
Core model handling routines.
"""
from __future__ import print_function

__all__ = [
    "list_models", "load_model", "load_model_info",
    "build_model", "precompile_dll",
    ]

from os.path import basename, dirname, join as joinpath
from glob import glob

import numpy as np

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
    DType = Union[None, str, np.dtype]
    from .kernel import KernelModel
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

def list_models():
    # type: () -> List[str]
    """
    Return the list of available models on the model path.
    """
    root = dirname(__file__)
    files = sorted(glob(joinpath(root, 'models', "[a-zA-Z]*.py")))
    available_models = [basename(f)[:-3] for f in files]
    return available_models

def isstr(s):
    # type: (Any) -> bool
    """
    Return True if *s* is a string-like object.
    """
    try: s + ''
    except Exception: return False
    return True

def load_model(model_name, dtype=None, platform='ocl'):
    # type: (str, DType, str) -> KernelModel
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
    # type: (modelinfo.ModelInfo, DType, str) -> KernelModel
    """
    Prepare the model for the default execution platform.

    This will return an OpenCL model, a DLL model or a python model depending
    on the model and the computing platform.

    *model_info* is the model definition structure returned from
    :func:`load_model_info`.

    *dtype* indicates whether the model should use single or double precision
    for the calculation. Any valid numpy single or double precision identifier
    is valid, such as 'single', 'f', 'f32', or np.float32 for single, or
    'double', 'd', 'f64'  and np.float64 for double.  If *None*, then use
    'single' unless the model defines single=False.

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

    ## for debugging:
    ##  1. uncomment open().write so that the source will be saved next time
    ##  2. run "python -m sasmodels.direct_model $MODELNAME" to save the source
    ##  3. recomment the open.write() and uncomment open().read()
    ##  4. rerun "python -m sasmodels.direct_model $MODELNAME"
    ##  5. uncomment open().read() so that source will be regenerated from model
    # open(model_info.name+'.c','w').write(source)
    # source = open(model_info.name+'.cl','r').read()
    source = generate.make_source(model_info)
    if dtype is None:
        dtype = 'single' if model_info.single else 'double'
    if callable(model_info.Iq):
        return kernelpy.PyModel(model_info)
    if (platform == "dll"
            or not HAVE_OPENCL
            or not kernelcl.environment().has_type(dtype)):
        return kerneldll.load_dll(source, model_info, dtype)
    else:
        return kernelcl.GpuModel(source, model_info, dtype)

def precompile_dll(model_name, dtype="double"):
    # type: (str, DType) -> Optional[str]
    """
    Precompile the dll for a model.

    Returns the path to the compiled model, or None if the model is a pure
    python model.

    This can be used when build the windows distribution of sasmodels
    (which may be missing the OpenCL driver and the dll compiler), or
    otherwise sharing models with windows users who do not have a compiler.

    See :func:`sasmodels.kerneldll.make_dll` for details on controlling the
    dll path and the allowed floating point precision.
    """
    model_info = load_model_info(model_name)
    source = generate.make_source(model_info)
    return kerneldll.make_dll(source, model_info, dtype=dtype) if source else None
