"""
Minimal example of calling a kernel for a specific set of q values.
"""

from numpy import logspace
from matplotlib import pyplot as plt
from sasmodels.core import load_model
from sasmodels.direct_model import call_kernel

model = load_model('cylinder')
q = logspace(-3, -1, 200)
kernel = model.make_kernel([q])
Iq = call_kernel(kernel, dict(radius=200.))
plt.loglog(q, Iq)
plt.xlabel('q (1/A)')
plt.ylabel('I(q)')
plt.title('Cylinder with radius 200.')
plt.show()
