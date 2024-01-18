import torch
import time
from numpy import logspace, sqrt
from matplotlib import pyplot as plt
from sasmodels.core import load_model
from sasmodels.direct_model import call_kernel,get_mesh
from sasmodels.details import make_kernel_args, dispersion_mesh




#cuda_src =  sas_3j1x_x + sphere_c

#Step 1. Define r and q vectors

model = load_model('_spherepy')
q = logspace(-3, -1, 20)
print("q",q[6])
kernel = model.make_kernel([q])

pars = {'radius': 20, 'radius_pd': 0.1, 'radius_pd_n':100, 'scale': 2}

t_before = time.time()
Iq = call_kernel(kernel, pars)
t_after = time.time()
total_np = t_after -t_before
print("Iq",Iq[6])
print("Tota Numpy: ",total_np)

t_before = time.time()
q_t = torch.logspace(start=-3, end=-1, steps=200)
kernel = model.make_kernel([q_t])

# call_kernel unwrap
calculator = kernel
cutoff=0.
mono=False

mesh = get_mesh(calculator.info, pars, dim=calculator.dim, mono=mono)
#print("in call_kernel: pars:", list(zip(*mesh))[0])
call_details, values, is_magnetic = make_kernel_args(calculator, mesh)
#print("in call_kernel: values:", values)
Iq_t = calculator(call_details, values, cutoff, is_magnetic)

t_after = time.time()
total_torch = t_after -t_before

print("pqt",q_t[6])
print("Iq_t",Iq_t[6])

print("Total Pytorch: ",total_torch)


