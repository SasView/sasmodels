r"""
Calculates the scattering from a fractal structure with a primary building
block of core-shell spheres, as opposed to just homogeneous spheres in
the fractal model.
This model could find use for aggregates of coated particles, or aggregates
of vesicles.

Definition
----------

.. math::

    I(q) = \text{background} + P(q)S(q)

The form factor $P(q)$ is that from core_shell model with $bkg$ = 0


.. math::

    P(q)=\frac{scale}{V_s}\left[3V_c(\rho_c-\rho_s)
    \frac{\sin(qr_c)-qr_c\cos(qr_c)}{(qr_c)^3}+
    3V_s(\rho_s-\rho_{solv})
    \frac{\sin(qr_s)-qr_s\cos(qr_s)}{(qr_s)^3}\right]^2


while the fractal structure factor $S(q)$ is

.. math::

    S(q) = \frac{D_f\Gamma(D_f-1)\sin((D_f-1)\tan^{-1}(q\xi))}
    {(qr_c)^{D_f}\left(1+\frac{1}{q^2\xi ^2} \right)^{\frac{D_f-1}{2}}}

where $D_f$ = frac_dim, |xi| = cor_length, $r_c$ = (core) radius, and
$scale$ = volume fraction.

The fractal structure is as documented in the fractal model.
Polydispersity of radius and thickness is provided for.

For 2D data: The 2D scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}

References
----------

See the core_shell and fractal model descriptions

"""

from numpy import pi, inf

name = "fractal_core_shell"
title = ""
description = """

"""
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["radius",      "Ang",        60.0, [0, inf],    "volume", "Sphere core radius"],
    ["thickness",   "Ang",        10.0, [0, inf],    "volume", "Sphere shell thickness"],
    ["sld_core",    "1e-6/Ang^2", 1.0,  [-inf, inf], "",       "Sphere core scattering length density"],
    ["sld_shell",   "1e-6/Ang^2", 2.0,  [-inf, inf], "",       "Sphere shell scattering length density"],
    ["sld_solvent", "1e-6/Ang^2", 3.0,  [-inf, inf], "",       "Solvent scattering length density"],
    ["volfraction", "",           1.0,  [0, inf],    "",       "Volume fraction of building block spheres"],
    ["frac_dim",    "",           2.0,  [-inf, inf], "",       "Fractal dimension"],
    ["cor_length",  "Ang",      100.0,  [0, inf],    "",       "Correlation length of fractal-like aggregates"]]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sph_j1c.c", "lib/sas_gamma.c", "lib/core_shell.c", "fractal_core_shell.c"]

demo = dict(scale=0.05,
            background=0,
            radius=20,
            thickness=5,
            sld_core=3.5,
            sld_shell=1.0,
            sld_solvent=6.35,
            volfraction=0.05,
            frac_dim=2.0,
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
    whole = 4.0 * pi / 3.0 * pow((radius + thickness), 3)
    core = 4.0 * pi / 3.0 * radius * radius * radius
    return whole, whole-core

tests = [[{'radius': 20.0, 'thickness': 10.0}, 'ER', 30.0],
         [{'radius': 20.0, 'thickness': 10.0}, 'VR', 0.703703704],

         # The SasView test result was 0.00169, with a background of 0.001
         [{'radius': 60.0,
           'thickness': 10.0,
           'sld_core': 1.0,
           'sld_shell': 2.0,
           'sld_solvent': 3.0,
           'background': 0.0
          }, 0.4, 0.00070126]]
