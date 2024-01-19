import torch
import time
from numpy import logspace, sqrt
from matplotlib import pyplot as plt
from sasmodels.core import load_model
from sasmodels.direct_model import call_kernel,get_mesh
from sasmodels.details import make_kernel_args

import sasmodels.kerneltorch as kt


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
#device = torch.device('mps')

print("device",device)

def make_kernel(model, q_vectors,device):
    """Instantiate the python kernel with input *q_vectors*"""
    q_input = kt.PyInput(q_vectors, dtype=torch.double)
    return kt.PyKernel(model.info, q_input, device = device)



model = load_model('_spherepy')
q = logspace(-3, -1, 200)
print("q",q[6])
kernel = model.make_kernel([q])

pars = {'radius': 200, 'radius_pd': 0.1, 'radius_pd_n':1000, 'sld':2, 'sld_pd': 0.1, 'sld_pd_n':100, 'scale': 2}

t_before = time.time()
Iq = call_kernel(kernel, pars)
t_after = time.time()
total_np = t_after -t_before
print("Iq",Iq[6])
print("Tota Numpy: ",total_np)

t_before = time.time()
q_t = torch.logspace(start=-3, end=-1, steps=200).to(device)
kernel = make_kernel(model, [q_t],device)
Iq_t = call_kernel(kernel, pars)

# call_kernel unwrap
#calculator = kernel
#cutoff=0.
#mono=False

#mesh = get_mesh(calculator.info, pars, dim=calculator.dim, mono=mono)
#print("in call_kernel: pars:", list(zip(*mesh))[0])
#call_details, values, is_magnetic = make_kernel_args(calculator, mesh)
#print("in call_kernel: values:", values)
#Iq_t = calculator(call_details, values, cutoff, is_magnetic)

t_after = time.time()
total_torch = t_after -t_before

print("pqt",q_t[6])
print("Iq_t",Iq_t[6])

print("Total Pytorch: ",total_torch)


