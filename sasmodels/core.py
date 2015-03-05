__all__ = ["list_models", "load_model_cl", "load_model_dll",
           "load_model_definition", ]

from os.path import basename, dirname, join as joinpath
from glob import glob

import numpy as np

from . import models
from . import weights

try:
    from .kernelcl import load_model as load_model_cl
except:
    # pylint: disable=invalid-name
    load_model_cl = None
from .kerneldll import load_model as load_model_dll

def list_models():
    root = dirname(__file__)
    files = sorted(glob(joinpath(root, 'models', "[a-zA-Z]*.py")))
    available_models = [basename(f)[:-3] for f in files]
    return available_models

def load_model_definition(model_name):
    __import__('sasmodels.models.'+model_name)
    model_definition = getattr(models, model_name, None)
    return model_definition

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

    Searches for "name", "name_pd", "name_pd_type", "name_pd_n", "name_pd_sigma"
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

def call_kernel(kernel, pars, cutoff=1e-5):
    fixed_pars = [pars.get(name, kernel.info['defaults'][name])
                  for name in kernel.fixed_pars]
    pd_pars = [get_weights(kernel.info, pars, name) for name in kernel.pd_pars]
    return kernel(fixed_pars, pd_pars, cutoff=cutoff)

def call_ER(info, pars):
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
    VR = info.get('VR', None)
    if VR is None:
        return 1.0
    else:
        vol_pars = [get_weights(info, pars, name)
                    for name in info['partype']['volume']]
        value, weight = dispersion_mesh(vol_pars)
        whole,part = VR(*value)
        return np.sum(weight*part)/np.sum(weight*whole)

