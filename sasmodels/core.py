__all__ = ["list_models", "load_model_cl", "load_model_dll",
           "load_model_definition", ]

from os.path import basename, dirname, join as joinpath
from glob import glob

import numpy as np

from . import models
from . import weights

try:
    from .kernelcl import load_model as load_model_cl
except Exception,exc:
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
    input = model.make_input(q_vectors)
    return model(input)

def get_weights(kernel, pars, name):
    """
    Generate the distribution for parameter *name* given the parameter values
    in *pars*.

    Searches for "name", "name_pd", "name_pd_type", "name_pd_n", "name_pd_sigma"
    """
    relative = name in kernel.info['partype']['pd-rel']
    limits = kernel.info['limits']
    disperser = pars.get(name+'_pd_type', 'gaussian')
    value = pars.get(name, kernel.info['defaults'][name])
    npts = pars.get(name+'_pd_n', 0)
    width = pars.get(name+'_pd', 0.0)
    nsigma = pars.get(name+'_pd_nsigma', 3.0)
    v,w = weights.get_weights(
        disperser, npts, width, nsigma,
        value, limits[name], relative)
    return v,w/np.sum(w)

def dispersion_mesh(pars):
    """
    Create a mesh grid of dispersion parameters and weights.

    Returns [p1,p2,...],w where pj is a vector of values for parameter j
    and w is a vector containing the products for weights for each
    parameter set in the vector.
    """
    values, weights = zip(*pars)
    if len(values) > 1:
        values = [v.flatten() for v in np.meshgrid(*values)]
        weights = np.vstack([v.flatten() for v in np.meshgrid(*weights)])
        weights = np.prod(weights, axis=0)
    return values, weights

def call_kernel(kernel, pars, cutoff=1e-5):
    fixed_pars = [pars.get(name, kernel.info['defaults'][name])
                  for name in kernel.fixed_pars]
    pd_pars = [get_weights(kernel, pars, name) for name in kernel.pd_pars]
    return kernel(fixed_pars, pd_pars, cutoff=cutoff)

def call_ER(kernel, pars):
    ER = kernel.info.get('ER', None)
    if ER is None:
        return 1.0
    else:
        vol_pars = [get_weights(kernel, pars, name)
                    for name in kernel.info['partype']['volume']]
        values, weights = dispersion_mesh(vol_pars)
        fv = ER(*values)
        #print values[0].shape, weights.shape, fv.shape
        return np.sum(weights*fv) / np.sum(weights)

def call_VR(kernel, pars):
    VR = kernel.info.get('VR', None)
    if VR is None:
        return 1.0
    else:
        vol_pars = [get_weights(kernel, pars, name)
                    for name in kernel.info['partype']['volume']]
        values, weights = dispersion_mesh(vol_pars)
        whole,part = VR(*values)
        return np.sum(weights*part)/np.sum(weights*whole)

