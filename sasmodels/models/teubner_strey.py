r"""
Definition
----------

This model calculates the scattered intensity of a two-component system
using the Teubner-Strey model. Unlike :ref:`dab` this function generates
a peak.

.. math::

    I(q) = \frac{1}{a_2 + c_1 q^2 + c_2 q^4} + \text{background}

The parameters $a_2$, $c_1$ and $c_2$ can be used to determine the
characteristic domain size $d$,

.. math::

    d = 2\pi\left[\frac12\left(\frac{a_2}{c_2}\right)^{1/2}
                  + \frac14\frac{c_1}{c_2}\right]^{-1/2}


and the correlation length $\xi$,

.. math::

    \xi = \left[\frac12\left(\frac{a_2}{c_2}\right)^{1/2}
                  - \frac14\frac{c_1}{c_2}\right]^{-1/2}


For 2D data, scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


.. image:: img/image185.jpg

    1D plot using the default values (w/200 data point).

Reference
---------

M Teubner, R Strey, *J. Chem. Phys.*, 87 (1987) 3195

K V Schubert, R Strey, S R Kline and E W Kaler,
*J. Chem. Phys.*, 101 (1994) 5343

"""

import numpy as np
from numpy import inf, sqrt

name = "teubner_strey"
title = "Teubner-Strey model of microemulsions"
description = """\
   Scattering model class for the Teubner-Strey model given by
    Provide F(x) = 1/( a2 + c1 q^2+  c2 q^4 ) + background
    a2>0, c1<0, c2>0, 4 a2 c2 - c1^2 > 0
"""
category = "shape-independent"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["a2", "", 0.1, [0, inf], "",
               "a2"],
              ["c1", "1e-6/Ang^2", -30., [-inf, 0], "",
               "c1"],
              ["c2", "Ang", 5000., [0, inf], "volume",
               "c2"],
             ]


def form_volume(radius):
    return 1.0

def Iq(q, a2, c1, c2):
    return 1. / np.polyval([c2, c1, a2], q**2)
Iq.vectorized = True  # Iq accepts an array of Q values

def Iqxy(qx, qy, a2, c1, c2):
    return Iq(sqrt(qx**2 + qy**2), a2, c1, c2)
Iqxy.vectorized = True  # Iqxy accepts arrays of Qx, Qy values

# ER defaults to 0.0

# VR defaults to 1.0

demo = dict(scale=1, background=0, a2=0.1, c1=-30.0, c2=5000.0)
oldname = "TeubnerStreyModel"
oldpars = dict(a2='scale')

tests = [[{}, 0.2, 0.144927536232]]
