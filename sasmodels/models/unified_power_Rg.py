r"""
Definition
----------

The Beaucage model employs the empirical multiple level unified
Exponential/Power-law fit method developed by G. Beaucage. Four functions
are included so that 1, 2, 3, or 4 levels can be used. In addition a 0 level
has been added which simply calculates

.. math::

    I(q) = \text{scale} / q + \text{background}

The Beaucage method is able to reasonably approximate the scattering from
many different types of particles, including fractal clusters, random coils
(Debye equation), ellipsoidal particles, etc.

The empirical fit function is:

.. math::

    I(q) = \text{background}
    + \sum_{i=1}^N \Bigl[
        G_i \exp\Bigl(-\frac{q^2R_{gi}^2}{3}\Bigr)
       + B_i \exp\Bigl(-\frac{q^2R_{g(i+1)}^2}{3}\Bigr)
             \Bigl(\frac{1}{q_i^*}\Bigr)^{P_i} \Bigr]

where

.. math::

    q_i^* = \frac{q}{\operatorname{erf}^3(q R_{gi}/\sqrt{6}}


For each level, the four parameters $G_i$, $R_{gi}$, $B_i$ and $P_i$ must
be chosen.  Beaucage has an additional factor $k$ in the definition of
$q_i^*$ which is ignored here.

For example, to approximate the scattering from random coils (Debye equation),
set $R_{gi}$ as the Guinier radius, $P_i = 2$, and $B_i = 2 G_i / R_{gi}$

See the references for further information on choosing the parameters.

For 2D data: The 2D scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


References
----------

G Beaucage, *J. Appl. Cryst.*, 28 (1995) 717-728

G Beaucage, *J. Appl. Cryst.*, 29 (1996) 134-146

"""

from __future__ import division

import numpy as np
from numpy import inf, exp, sqrt, errstate
from scipy.special import erf

category = "shape-independent"
name = "unified_power_Rg"
title = "Unified Power Rg"
description = """
        The Beaucage model employs the empirical multiple level unified
        Exponential/Power-law fit method developed by G. Beaucage. Four functions
        are included so that 1, 2, 3, or 4 levels can be used.
        """

# pylint: disable=bad-whitespace, line-too-long
parameters = [
    ["level",     "",     1,      [0, 6], "", "Level number"],
    ["rg[level]", "Ang",  15.8,   [0, inf], "", "Radius of gyration"],
    ["power[level]", "",  4,      [-inf, inf], "", "Power"],
    ["B[level]",  "1/cm", 4.5e-6, [-inf, inf], "", ""],
    ["G[level]",  "1/cm", 400,    [0, inf], "", ""],
    ]
# pylint: enable=bad-whitespace, line-too-long

def Iq(q, level, rg, power, B, G):
    level = int(level + 0.5)
    if level == 0:
        with errstate(divide='ignore'):
            return 1./q

    with errstate(divide='ignore', invalid='ignore'):
        result = np.zeros(q.shape, 'd')
        for i in range(level):
            exp_now = exp(-(q*rg[i])**2/3.)
            pow_now = (erf(q*rg[i]/sqrt(6.))**3/q)**power[i]
            if i < level-1:
                exp_next = exp(-(q*rg[i+1])**2/3.)
            else:
                exp_next = 1
            result += G[i]*exp_now + B[i]*exp_next*pow_now

    result[q == 0] = np.sum(G[:level])
    return result

Iq.vectorized = True

demo = dict(
    level=2,
    rg=[15.8, 21],
    power=[4, 2],
    B=[4.5e-6, 0.0006],
    G=[400, 3],
    scale=1.,
    background=0.,
)
