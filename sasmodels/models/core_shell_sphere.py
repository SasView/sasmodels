r"""
This model provides the form factor, $P(q)$, for a spherical particle with a core-shell structure.
The form factor is normalized by the particle volume.

Definition
----------

The 1D scattering intensity is calculated in the following way (Guinier, 1955)

.. math::

    P(q) = \frac{\text{scale}}{V} F^2(q) + \text{background}

where

.. math::
    F^2(q)=\frac{3}{V_s}\left[V_c(\rho_c-\rho_s)\frac{\sin(qr_c)-qr_c\cos(qr_c)}{(qr_c)^3}+
    V_s(\rho_s-\rho_{solv})\frac{\sin(qr_s)-qr_s\cos(qr_s)}{(qr_s)^3}\right]

where $V_s$ is the volume of the outer shell, $V_c$ is
the volume of the core, $r_s$ is the radius of the shell, $r_c$ is the radius of the
core, $\rho_c$ is the scattering length density of the core, $\rho_s$ is the scattering length
density of the shell, $\rho_{solv}$ is the scattering length density of the solvent.

The 2D scattering intensity is the same as $P(q)$ above, regardless of the
orientation of the $q$ vector.

NB: The outer most radius (ie, = radius + thickness) is used as the
effective radius for $S(Q)$ when $P(Q) \cdot S(Q)$ is applied.

Reference
---------

A Guinier and G Fournet, *Small-Angle Scattering of X-Rays*, John Wiley and Sons, New York, (1955)

Validation
----------

Validation of our code was done by comparing the output of the 1D model to the output of
the software provided by NIST (Kline, 2006). Figure 1 shows a comparison of the output of
our model and the output of the NIST software.

.. image:: img/core_shell_sphere_1d.jpg

    Figure 1: Comparison of the SasView scattering intensity for a core-shell sphere with
    the output of the NIST SANS analysis software. The parameters were set to:
    *scale* = 1.0, *radius* = 60 , *contrast* = 1e-6 |Ang^-2|, and
    *background* = 0.001 |cm^-1|.
"""

from numpy import pi, inf

name = "core_shell_sphere"
title = "Form factor for a monodisperse spherical particle with particle with a core-shell structure."
description = """
    F^2(q) = 3/V_s [V_c (core_sld-shell_sld) (sin(q*radius)-q*radius*cos(q*radius))/(q*radius)^3 
                   + V_s (shell_sld-solvent_sld) (sin(q*r_s)-q*r_s*cos(q*r_s))/(q*r_s)^3]

            V_s: Volume of the sphere shell
            V_c: Volume of the sphere core
            r_s: Shell radius = radius + thickness
"""
category = "shape:sphere"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["radius",      "Ang",        60.0, [0, inf],    "volume", "Sphere core radius"],
              ["thickness",   "Ang",        10.0, [0, inf],    "volume", "Sphere shell thickness"],
              ["core_sld",    "1e-6/Ang^2", 1.0,  [-inf, inf], "",       "Sphere core scattering length density"],
              ["shell_sld",   "1e-6/Ang^2", 2.0,  [-inf, inf], "",       "Sphere shell scattering length density"],
              ["solvent_sld", "1e-6/Ang^2", 3.0,  [-inf, inf],  "",      "Solvent scattering length density"]]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sph_j1c.c", "core_shell_sphere.c"]

demo = dict(scale=1, background=0, radius=60, thickness=10,
            core_sld=1.0, shell_sld=2.0, solvent_sld=0.0)

oldname = 'CoreShellModel'
oldpars = {}

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
    whole = 4.0 * pi / 3.0 * pow((radius + thickness), 3)
    core = 4.0 * pi / 3.0 * radius * radius * radius
    return whole, whole - core

tests = [[{'radius': 20.0, 'thickness': 10.0}, 'ER', 30.0],
         [{'radius': 20.0, 'thickness': 10.0}, 'VR', 0.703703704],

         # The SasView test result was 0.00169, with a background of 0.001
         [{'radius': 60.0,
           'thickness': 10.0,
           'core_sld': 1.0,
           'shell_sld':2.0,
           'solvent_sld':3.0,
           'background':0.0
          }, 0.4, 0.000698838]]
