r"""
This model provides the form factor, $P(q)$, for a 'pringle' or 'saddle-shaped'
object (a hyperbolic paraboloid).

.. figure:: img/pringles_fig1.png

    (Graphic from Matt Henderson, matt@matthen.com)

The returned value is in units of cm-1, on absolute scale.

Definition
----------

The form factor calculated is

.. math::

    I(q) = (\Delta \rho )^2 V \int^{\pi/2}_0 d\psi \sin{\psi} sinc^2
    \left( \frac{qd\cos{\psi}}{2} \right)
    \left[ \left( S^2_0+C^2_0\right) + 2\sum_{n=1}^{\infty}
     \left( S^2_0+C^2_0\right) \right]

where

.. math::

    C_n = \int^{R}_{0} r dr\cos{qr^2\alpha \cos{\psi}}
    J_n\left( qr^2\beta \cos{\psi}\right)
    J_{2n}\left( qr \sin{\psi}\right)

.. math::

    S_n = \int^{R}_{0} r dr\sin{qr^2\alpha \cos{\psi}}
    J_n\left( qr^2\beta \cos{\psi}\right)
    J_{2n}\left( qr \sin{\psi}\right)


.. figure:: img/pringle-vs-cylinder.png

    1D plot using the default values (with 150 data points).

Reference
---------

S Alexandru Rautu, Private Communication.

"""

from numpy import inf, pi

name = "pringles"
title = "Pringles model for K Edler. Represents a disc that is bent in two directions."
description = """\

"""
category = "shape:cylinder"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["radius",      "Ang",         60.0,   [0, inf],    "volume", "Pringle radius"],
    ["thickness",   "Ang",         10.0,   [0, inf],    "volume", "Thickness of pringle"],
    ["alpha",       "",            0.001,  [-inf, inf], "", "Curvature parameter"],
    ["beta",        "",            0.02,   [-inf, inf], "", "Curvature paramter"],
    ["sld_pringle", "1e-6/Ang^2",  1.0,    [-inf, inf], "", "Fractal exponent"],
    ["sld_solvent", "1e-6/Ang^2",  6.3,    [-inf, inf], "", "Correlation length"]
    ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/polevl.c", "lib/sas_J0.c", "lib/sas_J1.c", "lib/sas_JN.c", "lib/gauss76.c", "pringles.c"]

def ER(radius, thickness):
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

oldname = 'PringlesModel'
oldpars = dict(background='background',
                scale='scale',
                radius = 'radius',
                thickness = 'thickness',
                alpha = 'alpha',
                beta='beta',
                pringle_sld='sld_pringle',
                solvent_sld='sld_solvent')

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
