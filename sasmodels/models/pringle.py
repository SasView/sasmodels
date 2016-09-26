r"""
Definition
----------

The form factor for this bent disc is essentially that of a hyperbolic
paraboloid and calculated as

.. math::

    P(q) = (\Delta \rho )^2 V \int^{\pi/2}_0 d\psi \sin{\psi} sinc^2
    \left( \frac{qd\cos{\psi}}{2} \right)
    \left[ \left( S^2_0+C^2_0\right) + 2\sum_{n=1}^{\infty}
     \left( S^2_n+C^2_n\right) \right]

where

.. math::

    C_n = \int^{R}_{0} r dr\cos(qr^2\alpha \cos{\psi})
    J_n\left( qr^2\beta \cos{\psi}\right)
    J_{2n}\left( qr \sin{\psi}\right)

.. math::

    S_n = \int^{R}_{0} r dr\sin(qr^2\alpha \cos{\psi})
    J_n\left( qr^2\beta \cos{\psi}\right)
    J_{2n}\left( qr \sin{\psi}\right)

and $\Delta \rho \text{ is } \rho_{pringle}-\rho_{solvent}$, $V$ is the volume of
the disc, $\psi$ is the angle between the normal to the disc and the q vector,
$d$ and $R$ are the "pringle" thickness and radius respectively, $\alpha$ and
$\beta$ are the two curvature parameters, and $J_n$ is the n\ :sup:`th` order
Bessel function of the first kind.

.. figure:: img/pringles_fig1.png

    Schematic of model shape (Graphic from Matt Henderson, matt@matthen.com)

Reference
---------

Karen Edler, Universtiy of Bath, Private Communication. 2012.
Derivation by Stefan Alexandru Rautu.

**Author:** Andrew Jackson **on:** 2008

**Last Modified by:** Wojciech Wpotrzebowski **on:** March 20, 2016

**Last Reviewed by:** Paul Butler **on:** March 21, 2016

"""

from numpy import inf, pi

name = "pringle"
title = "The Pringle model provides the form factor, $P(q)$, for a 'pringle' \
or 'saddle-shaped' disc that is bent in two directions."
description = """\

"""
category = "shape:cylinder"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["radius",      "Ang",         60.0,   [0, inf],    "volume", "Pringle radius"],
    ["thickness",   "Ang",         10.0,   [0, inf],    "volume", "Thickness of pringle"],
    ["alpha",       "",            0.001,  [-inf, inf], "volume", "Curvature parameter alpha"],
    ["beta",        "",            0.02,   [-inf, inf], "volume", "Curvature paramter beta"],
    ["sld_pringle", "1e-6/Ang^2",  1.0,    [-inf, inf], "sld", "Pringle sld"],
    ["sld_solvent", "1e-6/Ang^2",  6.3,    [-inf, inf], "sld", "Solvent sld"]
    ]
# pylint: enable=bad-whitespace, line-too-long


source = ["lib/polevl.c", "lib/sas_J0.c", "lib/sas_J1.c", \
          "lib/sas_JN.c", "lib/gauss76.c", "pringle.c"]

def ER(radius, thickness, alpha, beta):
    """
    Return equivalent radius (ER)
    """
    ddd = 0.75 * radius * (2 * radius * thickness + (thickness + radius) \
                           * (thickness + pi * radius))
    return 0.5 * (ddd) ** (1. / 3.)

demo = dict(background=0.0,
            scale=1.0,
            radius=60.0,
            thickness=10.0,
            alpha=0.001,
            beta=0.02,
            sld_pringle=1.0,
            sld_solvent=6.35)

tests = [
    [{'scale' : 1.0,
      'radius': 60.0,
      'thickness': 10.0,
      'alpha': 0.001,
      'beta': 0.02,
      'sld_pringle': 1.0,
      'sld_solvent': 6.3,
      'background': 6.3,
     }, 0.1, 16.185532],

    [{'scale' : 1.0,
      'radius': 60.0,
      'thickness': 10.0,
      'alpha': 0.001,
      'beta': 0.02,
      'sld_pringle': 1.0,
      'sld_solvent': 6.3,
      'background': 6.3,
     }, 0.01, 297.153496],

    [{'scale' : 1.0,
      'radius': 60.0,
      'thickness': 10.0,
      'alpha': 0.001,
      'beta': 0.02,
      'sld_pringle': 1.0,
      'sld_solvent': 6.3,
      'background': 6.3,
     }, 0.001, 324.021256415],

    [{'scale' : 1.0,
      'radius': 60.0,
      'thickness': 10.0,
      'alpha': 0.001,
      'beta': 0.02,
      'sld_pringle': 1.0,
      'sld_solvent': 6.3,
      'background': 6.3,
     }, (0.001, 90.0), 6.30000026876],
]
