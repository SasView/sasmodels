"""
Minimal example of calling a kernel for a specific set of q values.
"""

from numpy import logspace, sqrt
from matplotlib import pyplot as plt
from sasmodels.core import load_model
from sasmodels.direct_model import call_kernel, call_Fq

model = load_model('cylinder')
q = logspace(-3, -1, 200)
kernel = model.make_kernel([q])
pars = {'radius': 200, 'radius_pd': 0.1, 'scale': 2}
Iq = call_kernel(kernel, pars)
F, Fsq, Reff, V, Vratio = call_Fq(kernel, pars)
plt.loglog(q, Iq, label='2 I(q)')
plt.loglog(q, F**2/V, label='<F(q)>^2/V')
plt.loglog(q, Fsq/V, label='<F^2(q)>/V')
plt.xlabel('q (1/A)')
plt.ylabel('I(q) (1/cm)')
plt.title('Cylinder with radius 200.')
plt.legend()
plt.show()
