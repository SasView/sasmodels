"""
Core model handling routines.
"""
__all__ = ["list_models", "load_model_definition",  "precompile_dll",
           "load_model", "make_kernel", "call_kernel", "call_ER", "call_VR" ]

from os.path import basename, dirname, join as joinpath
from glob import glob

import numpy as np

from . import models
from . import weights
from . import generate

from . import kernelpy
from . import kerneldll
try:
    from . import kernelcl
    HAVE_OPENCL = True
except:
    HAVE_OPENCL = False


def list_models():
    """
    Return the list of available models on the model path.
    """
    root = dirname(__file__)
    files = sorted(glob(joinpath(root, 'models', "[a-zA-Z]*.py")))
    available_models = [basename(f)[:-3] for f in files]
    return available_models


def load_model_definition(model_name):
    """
    Load a model definition given the model name.
    """
    __import__('sasmodels.models.'+model_name)
    model_definition = getattr(models, model_name, None)
    return model_definition


def precompile_dll(model_name, dtype="double"):
    """
    Precompile the dll for a model.

    Returns the path to the compiled model.

    This can be used when build the windows distribution of sasmodels
    (which may be missing the OpenCL driver and the dll compiler), or
    otherwise sharing models with windows users who do not have a compiler.

    See :func:`sasmodels.kerneldll.make_dll` for details on controlling the
    dll path and the allowed floating point precision.
    """
    model_definition = load_model_definition(model_name)
    source, info = generate.make(model_definition)
    return kerneldll.make_dll(source, info, dtype=dtype)


def isstr(s):
    try: s + ''
    except: return False
    return True

def load_model(model_definition, dtype="single", platform="ocl"):
    """
    Prepare the model for the default execution platform.

    This will return an OpenCL model, a DLL model or a python model depending
    on the model and the computing platform.

    *model_definition* is the python module which defines the model.  If the
    model name is given instead, then :func:`load_model_definition` will be
    called with the model name.

    *dtype* indicates whether the model should use single or double precision
    for the calculation. Any valid numpy single or double precision identifier
    is valid, such as 'single', 'f', 'f32', or np.float32 for single, or
    'double', 'd', 'f64'  and np.float64 for double.

    *platform* should be "dll" to force the dll to be used for C models,
    otherwise it uses the default "ocl".
    """
    if isstr(model_definition):
        model_definition = load_model_definition(model_definition)
    source, info = generate.make(model_definition)
    if callable(info.get('Iq', None)):
        return kernelpy.PyModel(info)

    ## for debugging:
    ##  1. uncomment open().write so that the source will be saved next time
    ##  2. run "python -m sasmodels.direct_model $MODELNAME" to save the source
    ##  3. recomment the open.write() and uncomment open().read()
    ##  4. rerun "python -m sasmodels.direct_model $MODELNAME"
    ##  5. uncomment open().read() so that source will be regenerated from model
    # open(info['name']+'.c','w').write(source)
    # source = open(info['name']+'.cl','r').read()

    dtype = np.dtype(dtype)
    if (platform=="dll"
            or not HAVE_OPENCL
            or (dtype == np.float64 and not kernelcl.environment().has_double)):
        return kerneldll.load_dll(source, info, dtype)
    else:
        return kernelcl.GpuModel(source, info, dtype)

def make_kernel(model, q_vectors):
    """
    Return a computation kernel from the model definition and the q input.
    """
    model_input = model.make_input(q_vectors)
    return model(model_input)

def get_weights(info, pars, name):
    """
    Generate the distribution for parameter *name* given the parameter values
    in *pars*.

    Uses "name", "name_pd", "name_pd_type", "name_pd_n", "name_pd_sigma"
    from the *pars* dictionary for parameter value and parameter dispersion.
    """
    relative = name in info['partype']['pd-rel']
    limits = info['limits']
    disperser = pars.get(name+'_pd_type', 'gaussian')
    value = pars.get(name, info['defaults'][name])
    npts = pars.get(name+'_pd_n', 0)
    width = pars.get(name+'_pd', 0.0)
    nsigma = pars.get(name+'_pd_nsigma', 3.0)
    value,weight = weights.get_weights(
        disperser, npts, width, nsigma,
        value, limits[name], relative)
    return value,weight/np.sum(weight)

def dispersion_mesh(pars):
    """
    Create a mesh grid of dispersion parameters and weights.

    Returns [p1,p2,...],w where pj is a vector of values for parameter j
    and w is a vector containing the products for weights for each
    parameter set in the vector.
    """
    value, weight = zip(*pars)
    if len(value) > 1:
        value = [v.flatten() for v in np.meshgrid(*value)]
        weight = np.vstack([v.flatten() for v in np.meshgrid(*weight)])
        weight = np.prod(weight, axis=0)
    return value, weight

def call_kernel(kernel, pars, cutoff=0):
    """
    Call *kernel* returned from :func:`make_kernel` with parameters *pars*.

    *cutoff* is the limiting value for the product of dispersion weights used
    to perform the multidimensional dispersion calculation more quickly at a
    slight cost to accuracy. The default value of *cutoff=0* integrates over
    the entire dispersion cube.  Using *cutoff=1e-5* can be 50% faster, but
    with an error of about 1%, which is usually less than the measurement
    uncertainty.
    """
    fixed_pars = [pars.get(name, kernel.info['defaults'][name])
                  for name in kernel.fixed_pars]
    pd_pars = [get_weights(kernel.info, pars, name) for name in kernel.pd_pars]
    return kernel(fixed_pars, pd_pars, cutoff=cutoff)

def call_ER(info, pars):
    """
    Call the model ER function using *pars*.

    *info* is either *model.info* if you have a loaded model, or *kernel.info*
    if you have a model kernel prepared for evaluation.
    """
    ER = info.get('ER', None)
    if ER is None:
        return 1.0
    else:
        vol_pars = [get_weights(info, pars, name)
                    for name in info['partype']['volume']]
        value, weight = dispersion_mesh(vol_pars)
        individual_radii = ER(*value)
        #print values[0].shape, weights.shape, fv.shape
        return np.sum(weight*individual_radii) / np.sum(weight)

def call_VR(info, pars):
    """
    Call the model VR function using *pars*.

    *info* is either *model.info* if you have a loaded model, or *kernel.info*
    if you have a model kernel prepared for evaluation.
    """
    VR = info.get('VR', None)
    if VR is None:
        return 1.0
    else:
        vol_pars = [get_weights(info, pars, name)
                    for name in info['partype']['volume']]
        value, weight = dispersion_mesh(vol_pars)
        whole,part = VR(*value)
        return np.sum(weight*part)/np.sum(weight*whole)

