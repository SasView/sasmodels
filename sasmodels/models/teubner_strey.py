r"""
Definition
----------

This model calculates the scattered intensity of a two-component system
using the Teubner-Strey model. Unlike :ref:`dab` this function generates
a peak. A two-phase material can be characterised by two length scales -
a correlation length and a domain size (periodicity).

The original paper by Teubner and Strey defined the function as:

.. math::

    I(q) \propto \frac{1}{a_2 + c_1 q^2 + c_2 q^4} + \text{background}

where the parameters $a_2$, $c_1$ and $c_2$ are defined in terms of the
periodicity, $d$, and correlation length $\xi$ as:

.. math::

    a_2 &= \biggl[1+\bigl(\frac{2\pi\xi}{d}\bigr)^2\biggr]\\
    c_1 &= -2\xi^2\bigl(\frac{2\pi\xi}{d}\bigr)^2+2\xi^2\\
    c_2 &= \xi^4

and thus, the periodicity, $d$ is given by

.. math::

    d = 2\pi\left[\frac12\left(\frac{a_2}{c_2}\right)^{1/2}
                  - \frac14\frac{c_1}{c_2}\right]^{-1/2}

and the correlation length, $\xi$, is given by

.. math::

    \xi = \left[\frac12\left(\frac{a_2}{c_2}\right)^{1/2}
                  + \frac14\frac{c_1}{c_2}\right]^{-1/2}

Here the model is parameterised in terms of  $d$ and $\xi$ and with an explicit
volume fraction for one phase, $\phi_a$, and contrast,
$\delta\rho^2 = (\rho_a - \rho_b)^2$ :

.. math::

    I(q) = \frac{8\pi\phi_a(1-\phi_a)(\Delta\rho)^2c_2/\xi}
        {a_2 + c_1q^2 + c_2q^4}

where :math:`8\pi\phi_a(1-\phi_a)(\Delta\rho)^2c_2/\xi` is the constant of
proportionality from the first equation above.

In the case of a microemulsion, $a_2 > 0$, $c_1 < 0$, and $c_2 >0$.

For 2D data, scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}

References
----------

M Teubner, R Strey, *J. Chem. Phys.*, 87 (1987) 3195

K V Schubert, R Strey, S R Kline and E W Kaler,
*J. Chem. Phys.*, 101 (1994) 5343

H Endo, M Mihailescu, M. Monkenbusch, J Allgaier, G Gompper, D Richter,
B Jakobs, T Sottmann, R Strey, and I Grillo, *J. Chem. Phys.*, 115 (2001), 580
"""

import numpy as np
from numpy import inf,power,pi

name = "teubner_strey"
title = "Teubner-Strey model of microemulsions"
description = """\
    Calculates scattering according to the Teubner-Strey model
"""
category = "shape-independent"

#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["volfraction_a", "", 0.5, [0, 1.0], "", "Volume fraction of phase a"],
    ["sld_a", "1e-6/Ang^2", 0.3, [-inf, inf], "", "SLD of phase a"],
    ["sld_b", "1e-6/Ang^2", 6.3, [-inf, inf], "", "SLD of phase b"],
    ["d", "Ang", 100.0, [0, inf], "", "Domain size (periodicity)"],
    ["xi", "Ang", 30.0, [0, inf], "", "Correlation length"],
    ]

def Iq(q, volfraction, sld, sld_solvent,d,xi):
    """SAS form"""
    drho2 = (sld-sld_solvent)*(sld-sld_solvent)
    a2 = power(1.0+power(2.0*pi*xi/d,2.0),2.0)
    c1 = -2.0*xi*xi*power(2.0*pi*xi/d,2.0)+2*xi*xi
    c2 = power(xi,4.0)
    prefactor = 8.0*pi*volfraction*(1.0-volfraction)*drho2*c2/xi
    #k2 = (2.0*pi/d)*(2.0*pi/d)
    #xi2 = 1/(xi*xi)
    #q2 = q*q
    #result = prefactor/((xi2+k2)*(xi2+k2)+2.0*(xi2-k2)*q2+q2*q2)
    return 1.0e-4*prefactor / np.polyval([c2, c1, a2], q**2)

Iq.vectorized = True  # Iq accepts an array of q values

demo = dict(scale=1, background=0, volfraction_a=0.5,
                     sld_a=0.3, sld_b=6.3,
                     d=100.0, xi=30.0)
tests = [[{}, 0.06, 41.5918888453]]
