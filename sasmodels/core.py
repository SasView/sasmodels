"""
Core model handling routines.
"""

from os.path import basename, dirname, join as joinpath, splitext
from glob import glob
import imp

import numpy as np

from . import models
from . import weights
from . import generate
# TODO: remove circular references between product and core
# product uses call_ER/call_VR, core uses make_product_info/ProductModel
#from . import product
from . import mixture
from . import kernelpy
from . import kerneldll
try:
    from . import kernelcl
    HAVE_OPENCL = True
except:
    HAVE_OPENCL = False

__all__ = [
    "list_models", "load_model_info", "precompile_dll",
    "build_model", "call_kernel", "call_ER_VR",
]

try:
    # Python 3.5 and up
    from importlib.util import spec_from_file_location, module_from_spec
    def load_module(fullname, path):
        spec = spec_from_file_location(fullname, path)
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
except ImportError:
    # CRUFT: python 2
    import imp
    def load_module(fullname, path):
        module = imp.load_source(fullname, path)
        return module

def list_models():
    """
    Return the list of available models on the model path.
    """
    root = dirname(__file__)
    files = sorted(glob(joinpath(root, 'models', "[a-zA-Z]*.py")))
    available_models = [basename(f)[:-3] for f in files]
    return available_models

def isstr(s):
    """
    Return True if *s* is a string-like object.
    """
    try: s + ''
    except: return False
    return True

def load_model(model_name, **kw):
    """
    Load model info and build model.
    """
    if model_name.endswith('.py'):
        model_info = load_model_info_from_path(model_name)
    else:
        model_info = load_model_info(model_name)
    return build_model(model_info, **kw)

def load_model_info_from_path(path):
    # Pull off the last .ext if it exists; there may be others
    name = basename(splitext(path)[0])

    # Not cleaning name since don't need to be able to reload this
    # model later
    # Should probably turf the model from sys.modules after we are done...

    # Placing the model in the 'sasmodels.custom' name space, even
    # though it doesn't actually exist.  imp.load_source doesn't seem
    # to care.
    import sasmodels.custom
    kernel_module = load_module('sasmodels.custom.'+name, path)

    # Now that we have the module, we can load it as usual
    model_info = generate.make_model_info(kernel_module)
    return model_info

def load_model_info(model_name):
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
        from . import product
        # Note: currently have circular reference
        if len(parts) > 2:
            raise ValueError("use P*S to apply structure factor S to model P")
        P_info, Q_info = [load_model_info(p) for p in parts]
        return product.make_product_info(P_info, Q_info)

    #import sys; print "\n".join(sys.path)
    __import__('sasmodels.models.'+model_name)
    kernel_module = getattr(models, model_name, None)
    return generate.make_model_info(kernel_module)


def build_model(model_info, dtype=None, platform="ocl"):
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
    composition = model_info.get('composition', None)
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
    # open(model_info['name']+'.c','w').write(source)
    # source = open(model_info['name']+'.cl','r').read()
    source = generate.make_source(model_info)
    if dtype is None:
        dtype = 'single' if model_info['single'] else 'double'
    if callable(model_info.get('Iq', None)):
        return kernelpy.PyModel(model_info)
    if (platform == "dll"
            or not HAVE_OPENCL
            or not kernelcl.environment().has_type(dtype)):
        return kerneldll.load_dll(source, model_info, dtype)
    else:
        return kernelcl.GpuModel(source, model_info, dtype)

def precompile_dll(model_name, dtype="double"):
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


def get_weights(parameter, values):
    """
    Generate the distribution for parameter *name* given the parameter values
    in *pars*.

    Uses "name", "name_pd", "name_pd_type", "name_pd_n", "name_pd_sigma"
    from the *pars* dictionary for parameter value and parameter dispersion.
    """
    value = values.get(parameter.name, parameter.default)
    if parameter.type not in ('volume', 'orientation'):
        return [value], []
    relative = parameter.type == 'volume'
    limits = parameter.limits
    disperser = values.get(parameter.name+'_pd_type', 'gaussian')
    npts = values.get(parameter.name+'_pd_n', 0)
    width = values.get(parameter.name+'_pd', 0.0)
    nsigma = values.get(parameter.name+'_pd_nsigma', 3.0)
    value, weight = weights.get_weights(
        disperser, npts, width, nsigma, value, limits, relative)
    return value, weight / np.sum(weight)

def dispersion_mesh(pars):
    """
    Create a mesh grid of dispersion parameters and weights.

    Returns [p1,p2,...],w where pj is a vector of values for parameter j
    and w is a vector containing the products for weights for each
    parameter set in the vector.
    """
    value, weight = zip(*pars)
    value = [v.flatten() for v in np.meshgrid(*value)]
    weight = np.vstack([v.flatten() for v in np.meshgrid(*weight)])
    weight = np.prod(weight, axis=0)
    return value, weight

def call_kernel(kernel, pars, cutoff=0, mono=False):
    """
    Call *kernel* returned from *model.make_kernel* with parameters *pars*.

    *cutoff* is the limiting value for the product of dispersion weights used
    to perform the multidimensional dispersion calculation more quickly at a
    slight cost to accuracy. The default value of *cutoff=0* integrates over
    the entire dispersion cube.  Using *cutoff=1e-5* can be 50% faster, but
    with an error of about 1%, which is usually less than the measurement
    uncertainty.

    *mono* is True if polydispersity should be set to none on all parameters.
    """
    if mono:
        active = lambda name: False
    elif kernel.dim == '1d':
        pars_1d = set(p.name for p in kernel.info['parameters'].type['1d'])
        active = lambda name: name in pars_1d
    elif kernel.dim == '2d':
        pars_2d = set(p.name for p in kernel.info['parameters'].type['2d'])
        active = lambda name: name in pars_2d
    else:
        active = lambda name: True

    vw_pairs = [(get_weights(p, pars) if active(p.name) else ([p.default], []))
                for p in kernel.info['parameters']]
    values, weights = zip(*vw_pairs)

    if max([len(w) for w in weights]) > 1:
        details = generate.poly_details(kernel.info, weights)
    else:
        details = kernel.info['mono_details']

    weights, values = [np.hstack(v) for v in (weights, values)]

    weights = weights.astype(dtype=kernel.dtype)
    values = values.astype(dtype=kernel.dtype)
    return kernel(details, weights, values, cutoff)

def call_ER_VR(model_info, vol_pars):
    """
    Return effect radius and volume ratio for the model.

    *info* is either *kernel.info* for *kernel=make_kernel(model,q)*
    or *model.info*.

    *pars* are the parameters as expected by :func:`call_kernel`.
    """
    ER = model_info.get('ER', None)
    VR = model_info.get('VR', None)
    value, weight = dispersion_mesh(vol_pars)

    individual_radii = ER(*value) if ER else 1.0
    whole, part = VR(*value) if VR else (1.0, 1.0)

    effect_radius = np.sum(weight*individual_radii) / np.sum(weight)
    volume_ratio = np.sum(weight*part)/np.sum(weight*whole)
    return effect_radius, volume_ratio


def call_ER(model_info, values):
    """
    Call the model ER function using *values*. *model_info* is either
    *model.info* if you have a loaded model, or *kernel.info* if you
    have a model kernel prepared for evaluation.
    """
    ER = model_info.get('ER', None)
    if ER is None:
        return 1.0
    else:
        vol_pars = [get_weights(parameter, values)
                    for parameter in model_info['parameters']
                    if parameter.type == 'volume']
        value, weight = dispersion_mesh(vol_pars)
        individual_radii = ER(*value)
        #print(values[0].shape, weights.shape, fv.shape)
        return np.sum(weight*individual_radii) / np.sum(weight)

def call_VR(model_info, values):
    """
    Call the model VR function using *pars*.
    *info* is either *model.info* if you have a loaded model, or *kernel.info*
    if you have a model kernel prepared for evaluation.
    """
    VR = model_info.get('VR', None)
    if VR is None:
        return 1.0
    else:
        vol_pars = [get_weights(parameter, values)
                    for parameter in model_info['parameters']
                    if parameter.type == 'volume']
        value, weight = dispersion_mesh(vol_pars)
        whole, part = VR(*value)
        return np.sum(weight*part)/np.sum(weight*whole)

# TODO: remove call_ER, call_VR

