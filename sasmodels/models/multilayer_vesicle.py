r"""
Definition
----------

This model is a trivial extension of the core_shell_sphere function to include
*N* shells where the core is filled with solvent and the shells are interleaved
with layers of solvent. For *N = 1*, this returns the same as the vesicle model,
except for the normalisation, which here is to outermost volume.
The shell thicknessess and SLD are constant for all shells as expected for
a multilayer vesicle.

.. figure:: img/multi_shell_geometry.jpg

    Geometry of the multilayer_vesicle model.

See the :ref:`core-shell-sphere` model for more documentation.

The 1D scattering intensity is calculated in the following way (Guinier, 1955)

.. math::
    P(q) = \text{scale} \cdot \frac{V_f}{V_t} F^2(q) + \text{background}

for

.. math::
    F(q) = (\rho_\text{shell}-\rho_\text{solv}) \sum_{i=1}^{n_\text{pairs}}
        \left[
          3V(R_i)\frac{\sin(qR_i)-qR_i\cos(qR_i)}{(qR_i)^3} \\
          - 3V(R_i+t_s)\frac{\sin(q(R_i+t_s))-q(R_i+t_s)\cos(q(R_i+t_s))}{(q(R_i+t_s))^3}
        \right]

and

.. math::
     R_i = r_c + (i-1)(t_s + t_w)

where $V_f$ is the volume fraction of particles, $V_t$ is the volume of the
whole particle, $V(r)$ is the volume of a sphere of radius $r$, $r_c$ is the
radius of the core, $\rho_\text{shell}$ is the scattering length density of a
shell, $\rho_\text{solv}$ is the scattering length density of the solvent.

The outer most radius, $r_o = R_n + t_s$, is used for both the volume fraction
normalization and for the effective radius for *S(Q)* when $P(Q) * S(Q)$
is applied.

The 2D scattering intensity is the same as 1D, regardless of the orientation
of the q vector which is defined as:

.. math::

    q = \sqrt{q_x^2 + q_y^2}

For information about polarised and magnetic scattering, see
the :ref:`magnetism` documentation.

This code is based on the form factor calculations implemented in the NIST
Center for Neutron Research provided c-library (Kline, 2006).

References
----------

B Cabane, *Small Angle Scattering Methods*,
in *Surfactant Solutions: New Methods of Investigation*,
Ch.2, Surfactant Science Series Vol. 22, Ed. R Zana and M Dekker,
New York, (1987).

**Author:** NIST IGOR/DANSE **on:** pre 2010

**Last Modified by:** Piotr Rozyczko **on:** Feb 24, 2016

**Last Reviewed by:** Paul Butler **on:** March 20, 2016

"""

from numpy import inf

name = "multilayer_vesicle"
title = "P(Q) for a Multi-lamellar vesicle"
description = """
    multilayer_vesicle model parameters;
    scale : scale factor for abs intensity if needed else 1.0
    volfraction: volume fraction
    radius : Core radius of the multishell
    thick_shell: shell thickness
    thick_solvent: water thickness
    sld_solvent: solvent scattering length density
    sld: shell scattering length density
    n_pairs:number of "shell plus solvent" layer pairs
    background: incoherent background
        """
category = "shape:sphere"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["volfraction", "",  0.05, [0.0, 1],  "", "volume fraction of vesicles"],
    ["radius", "Ang", 60.0, [0.0, inf],  "", "radius of solvent filled core"],
    ["thick_shell", "Ang",        10.0, [0.0, inf],  "", "thickness of one shell"],
    ["thick_solvent", "Ang",        10.0, [0.0, inf],  "", "solvent thickness between shells"],
    ["sld_solvent",    "1e-6/Ang^2",  6.4, [-inf, inf], "sld", "solvent scattering length density"],
    ["sld",   "1e-6/Ang^2",  0.4, [-inf, inf], "sld", "Shell scattering length density"],
    ["n_pairs",     "",            2.0, [1.0, inf],  "", "Number of shell plus solvent layer pairs"],
    ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sph_j1c.c", "multilayer_vesicle.c"]

# TODO: the following line does nothing
polydispersity = ["radius", "n_pairs"]

demo = dict(scale=1, background=0,
            volfraction=0.05,
            radius=60.0,
            thick_shell=10.0,
            thick_solvent=10.0,
            sld_solvent=6.4,
            sld=0.4,
            n_pairs=2.0)

tests = [
    # Accuracy tests based on content in test/utest_other_models.py
    [{'radius': 60.0,
      'thick_shell': 10.0,
      'thick_solvent': 10.0,
      'sld_solvent': 6.4,
      'sld': 0.4,
      'n_pairs': 2.0,
      'scale': 1.0,
      'background': 0.001,
     }, 0.001, 122.1405],

    [{'volfraction': 1.0,
      'radius': 60.0,
      'thick_shell': 10.0,
      'thick_solvent': 10.0,
      'sld_solvent': 6.4,
      'sld': 0.4,
      'n_pairs': 2.0,
      'scale': 1.0,
      'background': 0.001,
     }, (0.001, 0.30903), 1.61873],
    ]
