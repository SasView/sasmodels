r"""
This model fits the Porod function

.. math::

    I(q) = C/q^4
    \\
    C = 2\pi (\Delta\rho)^2 S_v

to the data directly without any need for linearisation (cf. Log I(q) vs Log q).

Here $C$ is the scale factor and $S_v$ is the specific surface area (ie, surface area / volume)
of the sample, and $\Delta\rho$ is the contrast factor.

For 2D data: The 2D scattering intensity is calculated in the same way as 1D,
where the q vector is defined as

.. math::
    q = \sqrt{q_x^2+q_y^2}

References
----------

G Porod. *Kolloid Zeit*. 124 (1951) 83.
L A Feigin, D I Svergun, G W Taylor. *Structure Analysis by Small-Angle X-ray and Neutron Scattering*. Springer. (1987)
"""

from numpy import sqrt, power

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
    return 1.0/power(q, 4)

Iq.vectorized = True  # Iq accepts an array of q values

def Iqxy(qx, qy, *args):
    """
    @param qx:   Input q_x-value
    @param qy:   Input q_y-value
    @param args: Remaining arguments
    """
    return Iq(sqrt(qx ** 2 + qy ** 2), *args)

Iqxy.vectorized = True # Iqxy accepts an array of qx, qy values

demo = dict(scale=1.5, background=0.5)

oldname = "PorodModel"
oldpars = dict(scale='scale', background='background')

tests = [[{'scale': 0.00001, 'background':0.01}, 0.04, 3.916250]]
