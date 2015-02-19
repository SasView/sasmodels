import warnings

import numpy as np

from . import models
from . import weights

try:
    from .kernelcl import load_model
except ImportError,exc:
    warnings.warn(str(exc))
    warnings.warn("using ctypes instead")
    from .kerneldll import load_model

def load_model_definition(model_name):
    __import__('sasmodels.models.'+model_name)
    model_definition = getattr(models, model_name, None)
    return model_definition

# load_model is imported above.  It looks like the following
#def load_model(model_definition, dtype='single):
#    if kerneldll:
#        if source is newer than compiled: compile
#        load dll
#        return kernel
#    elif kernelcl:
#        compile source on context
#        return kernel


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

def call_kernel(kernel, pars):
    fixed_pars = [pars.get(name, kernel.info['defaults'][name])
                  for name in kernel.fixed_pars]
    pd_pars = [get_weights(kernel, pars, name) for name in kernel.pd_pars]
    return kernel(fixed_pars, pd_pars)

class DirectModel:
    def __init__(self, name, q_vectors, dtype='single'):
        self.model_definition = load_model_definition(name)
        self.model = load_model(self.model_definition, dtype=dtype)
        q_vectors = [np.ascontiguousarray(q,dtype=dtype) for q in q_vectors]
        self.kernel = make_kernel(self.model, q_vectors)
    def __call__(self, **pars):
        return call_kernel(self.kernel, pars)

def demo():
    import sys
    if len(sys.argv) < 3:
        print "usage: python -m sasmodels.direct_model modelname (q|qx,qy) par=val ..."
        sys.exit(1)
    model_name = sys.argv[1]
    values = [float(v) for v in sys.argv[2].split(',')]
    if len(values) == 1:
        q = values[0]
        q_vectors = [[q]]
    elif len(values) == 2:
        qx,qy = values
        q_vectors = [[qx],[qy]]
    else:
        print "use q or qx,qy"
        sys.exit(1)
    model = DirectModel(model_name, q_vectors)
    pars = dict((k,float(v))
                for pair in sys.argv[3:]
                for k,v in [pair.split('=')])
    Iq = model(**pars)
    print Iq[0]

if __name__ == "__main__":
    demo()
