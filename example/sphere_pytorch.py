"""
Minimal example of calling a kernel for a specific set of q values.

        npts = values.pop(parameter.name+'_pd_n', 0)
        width = values.pop(parameter.name+'_pd', 0.0)
        nsigma = values.pop(parameter.name+'_pd_nsigma', 3.0)
        distribution = values.pop(parameter.name+'_pd_type', 'gaussian')

"""
import time

import torch

from numpy import logspace, sqrt
from matplotlib import pyplot as plt
from sasmodels.core import load_model
from sasmodels.direct_model import call_kernel, get_mesh
from sasmodels.details import make_kernel_args, dispersion_mesh

import sasmodels.kerneltorch as kt

device = torch.device('mps')

def make_kernel(model, q_vectors):
    """Instantiate the python kernel with input *q_vectors*"""
    q_input = kt.PyInput(q_vectors, dtype=torch.float32)
    return kt.PyKernel(model.info, q_input)


model = load_model('_spherepy')

q = torch.logspace(-3, -1, 200).to(device)


#qq = logspace(-3, -1, 200)

kernel = make_kernel(model, [q])



pars = {'radius': 200, 'radius_pd': 0.2, 'radius_pd_n':10000, 'scale': 2}

#mesh = get_mesh(kernel.info, pars, dim=kernel.dim)
#print(mesh)

#call_details, values, is_magnetic = make_kernel_args(kernel, mesh)
#print(call_details)
#print(values)

t0 = time.time()
Iq = call_kernel(kernel, pars)
elapsed = time.time() - t0
print('Computation time:', elapsed)

print(Iq)
