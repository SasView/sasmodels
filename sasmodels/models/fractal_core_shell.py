r"""
Definition
----------
Calculates the scattering from a fractal structure with a primary building
block of core-shell spheres, as opposed to just homogeneous spheres in
the fractal model. It is an extension of the well known Teixeira\ [#teixeira]_
fractal model replacing the $P(q)$ of a solid sphere with that of a core-shell
sphere. This model could find use for aggregates of coated particles, or
aggregates of vesicles for example.

.. math::

    I(q) = P(q)S(q) + \text{background}

Where $P(q)$ is the core-shell form factor and $S(q)$ is the
Teixeira\ [#teixeira]_ fractal structure factor both of which are given again
below:

.. math::

    P(q) &= \frac{\phi}{V_s}\left[3V_c(\rho_c-\rho_s)
    \frac{\sin(qr_c)-qr_c\cos(qr_c)}{(qr_c)^3}+
    3V_s(\rho_s-\rho_{solv})
    \frac{\sin(qr_s)-qr_s\cos(qr_s)}{(qr_s)^3}\right]^2 \\
    S(q) &= 1 + \frac{D_f\ \Gamma\!(D_f-1)}{[1+1/(q\xi)^2]^{(D_f-1)/2}}
    \frac{\sin[(D_f-1)\tan^{-1}(q\xi)]}{(qr_s)^{D_f}}

where $\phi$ is the volume fraction of particles, $V_s$ is the volume of the
whole particle, $V_c$ is the volume of the core, $\rho_c$, $\rho_s$, and
$\rho_{solv}$ are the scattering length densities of the core, shell, and
solvent respectively, $r_c$ and $r_s$ are the radius of the core and the radius
of the whole particle respectively, $D_f$ is the fractal dimension, and $\xi$ the
correlation length.

Polydispersity of radius and thickness are also provided for.

This model does not allow for anisotropy and thus the 2D scattering intensity
is calculated in the same way as 1D, where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}

Our model is derived from the form factor calculations implemented in IGOR
macros by the NIST Center for Neutron Research\ [#Kline]_

References
----------

.. [#teixeira] J Teixeira, *J. Appl. Cryst.*, 21 (1988) 781-785
.. [#Kline]  S R Kline, *J Appl. Cryst.*, 39 (2006) 895

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Paul Butler and Paul Kienzle **Date:** November 27, 2016
* **Last Reviewed by:** Paul Butler and Paul Kienzle **Date:** November 27, 2016
"""

import numpy as np
from numpy import pi, inf

name = "fractal_core_shell"
title = "Scattering from a fractal structure formed from core shell spheres"
description = """\
    Model for fractal aggregates of core-shell primary particles. It is based on
    the Teixeira model for the S(q) of a fractal * P(q) for a core-shell sphere

    radius =  the radius of the core
    thickness = thickness of the shell
    thick_layer = thickness of a layer
    sld_core = the SLD of the core
    sld_shell = the SLD of the shell
    sld_solvent = the SLD of the solvent
    volfraction = volume fraction of core-shell particles
    fractal_dim = fractal dimension
    cor_length = correlation length of the fractal like aggretates
    """
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["radius",      "Ang",        60.0, [0.0, inf],  "volume", "Sphere core radius"],
    ["thickness",   "Ang",        10.0, [0.0, inf],  "volume", "Sphere shell thickness"],
    ["sld_core",    "1e-6/Ang^2", 1.0,  [-inf, inf], "sld",    "Sphere core scattering length density"],
    ["sld_shell",   "1e-6/Ang^2", 2.0,  [-inf, inf], "sld",    "Sphere shell scattering length density"],
    ["sld_solvent", "1e-6/Ang^2", 3.0,  [-inf, inf], "sld",    "Solvent scattering length density"],
    ["volfraction", "",           1.0,  [0.0, inf],  "",       "Volume fraction of building block spheres"],
    ["fractal_dim",    "",        2.0,  [0.0, 6.0],  "",       "Fractal dimension"],
    ["cor_length",  "Ang",      100.0,  [0.0, inf],  "",       "Correlation length of fractal-like aggregates"],
]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sas_3j1x_x.c", "lib/sas_gamma.c", "lib/core_shell.c",
          "lib/fractal_sq.c", "fractal_core_shell.c"]

def random():
    outer_radius = 10**np.random.uniform(0.7, 4)
    # Use a distribution with a preference for thin shell or thin core
    # Avoid core,shell radii < 1
    thickness = np.random.beta(0.5, 0.5)*(outer_radius-2) + 1
    radius = outer_radius - thickness
    cor_length = 10**np.random.uniform(0.7, 2)*outer_radius
    volfraction = 10**np.random.uniform(-3, -1)
    fractal_dim = 2*np.random.beta(3, 4) + 1
    pars = dict(
        #background=0, sld_block=1, sld_solvent=0,
        volfraction=volfraction,
        radius=radius,
        cor_length=cor_length,
        fractal_dim=fractal_dim,
    )
    return pars

demo = dict(scale=0.05,
            background=0,
            radius=20,
            thickness=5,
            sld_core=3.5,
            sld_shell=1.0,
            sld_solvent=6.35,
            volfraction=0.05,
            fractal_dim=2.0,
            cor_length=100.0)

def ER(radius, thickness):
    """
        Equivalent radius
        @param radius: core radius
        @param thickness: shell thickness
    """
    return radius + thickness

def VR(radius, thickness):
    """
        Volume ratio
        @param radius: core radius
        @param thickness: shell thickness
    """
    whole = 4.0/3.0 * pi * (radius + thickness)**3
    core = 4.0/3.0 * pi * radius**3
    return whole, whole-core

tests = [[{'radius': 20.0, 'thickness': 10.0}, 'ER', 30.0],
         [{'radius': 20.0, 'thickness': 10.0}, 'VR', 0.703703704]]

#         # The SasView test result was 0.00169, with a background of 0.001
#         # They are however wrong as we now know.  IGOR might be a more
#         # appropriate source. Otherwise will just have to assume this is now
#         # correct and self generate a correct answer for the future. Until we
#         # figure it out leave the tests commented out
#         [{'radius': 60.0,
#           'thickness': 10.0,
#           'sld_core': 1.0,
#           'sld_shell': 2.0,
#           'sld_solvent': 3.0,
#           'background': 0.0
#          }, 0.015211, 692.84]]
