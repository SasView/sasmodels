import warnings

import numpy as np

from .core import load_model_definition, make_kernel, call_kernel
from .core import load_model_cl as load_model
if load_model is None:
    warnings.warn("unable to load opencl; using ctypes instead")
    from .core import load_model_dll as load_model

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
