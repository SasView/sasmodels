r"""
This model fits the Porod function

.. math:: I(q) = C/q^4

to the data directly without any need for linearisation (cf. Log I(q) vs Log q).

Here $C = 2\pi (\Delta\rho)^2 S_v$ is the scale factor where $S_v$ is
the specific surface area (ie, surface area / volume) of the sample, and
$\Delta\rho$ is the contrast factor.

For 2D data: The 2D scattering intensity is calculated in the same way as 1D,
where the q vector is defined as

.. math:: q = \sqrt{q_x^2+q_y^2}

References
----------

G Porod. *Kolloid Zeit*. 124 (1951) 83.

L A Feigin, D I Svergun, G W Taylor. *Structure Analysis by Small-Angle
X-ray and Neutron Scattering*. Springer. (1987)
"""

import numpy as np
from numpy import inf, errstate

name = "porod"
title = "Porod function"
description = """\
          I(q) = scale/q^4 + background
"""

category = "shape-independent"

parameters = []

def Iq(q):
    """
    @param q: Input q-value
    """
    with errstate(divide='ignore'):
        return q**-4

Iq.vectorized = True  # Iq accepts an array of q values

def random():
    sld, solvent = np.random.uniform(-0.5, 12, size=2)
    radius = 10**np.random.uniform(1, 4.7)
    Vf = 10**np.random.uniform(-3, -1)
    scale = 1e-4 * Vf * 2*np.pi*(sld-solvent)**2/(3*radius)
    pars = dict(
        scale=scale,
    )
    return pars

demo = dict(scale=1.5, background=0.5)

tests = [
    [{'scale': 0.00001, 'background':0.01}, 0.04, 3.916250],
    [{}, 0.0, inf],
]
