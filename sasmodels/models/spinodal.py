r"""
Definition
----------

This model calculates the SAS signal of a phase separating solution
under spinodal decomposition. The scattering intensity $I(q)$ is calculated as

.. math::
    I(q) = I_{max}\frac{(1+\gamma/2)x^2}{\gamma/2+x^{2+\gamma}}+B

where $x=q/q_0$ and $B$ is a flat background. The characteristic structure
length scales with the correlation peak at $q_0$. The exponent $\gamma$ is
equal to $d+1$ with d the dimensionality of the off-critical concentration
mixtures. A transition to $\gamma=2d$ is seen near the percolation threshold
into the critical concentration regime.

References
----------

H. Furukawa. Dynamics-scaling theory for phase-separating unmixing mixtures:
Growth rates of droplets and scaling properties of autocorrelation functions.
Physica A 123,497 (1984).

Authorship and Verification
----------------------------

* **Author:** Dirk Honecker **Date:** Oct 7, 2016
"""

import numpy as np
from numpy import inf, errstate

name = "spinodal"
title = "Spinodal decomposition model"
description = """\
      I(q) = scale ((1+gamma/2)x^2)/(gamma/2+x^(2+gamma))+background

      List of default parameters:
      scale = scaling
      gamma = exponent
      x = q/q_0
      q_0 = correlation peak position [1/A]
      background = Incoherent background"""
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["gamma",      "",    3.0, [-inf, inf], "", "Exponent"],
              ["q_0",  "1/Ang",     0.1, [-inf, inf], "", "Correlation peak position"]
             ]
# pylint: enable=bad-whitespace, line-too-long

def Iq(q,
       gamma=3.0,
       q_0=0.1):
    """
    :param q:              Input q-value
    :param gamma:          Exponent
    :param q_0:            Correlation peak position
    :return:               Calculated intensity
    """

    with errstate(divide='ignore'):
        x = q/q_0
        inten = ((1 + gamma / 2) * x ** 2) / (gamma / 2 + x ** (2 + gamma))
    return inten
Iq.vectorized = True  # Iq accepts an array of q values

def random():
    pars = dict(
        scale=10**np.random.uniform(1, 3),
        gamma=np.random.uniform(0, 6),
        q_0=10**np.random.uniform(-3, -1),
    )
    return pars

demo = dict(scale=1, background=0,
            gamma=1, q_0=0.1)
