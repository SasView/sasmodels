r"""
Definition
----------

This model is a trivial extension of the core_shell_sphere function to include
*N* shells where the core is filled with solvent and the shells are interleaved
with layers of solvent. For *N = 1*, this returns the same as the vesicle model.
The shell thicknessess and SLD are constant across all shells as expected for
a multilayer vesicle.

.. figure:: img/multi_shell_geometry.jpg

    Geometry of the multilayer_vesicle model.

See the :ref:`core-shell` model for more documentation.

The 2D scattering intensity is the same as 1D, regardless of the orientation
of the q vector which is defined as:

.. math::

    q = \sqrt{q_x^2 + q_y^2}

.. note:
    The outer most radius
    $radius + n_pairs * thicn_shell + (n_pairs - 1) * thick_solvent$
    is used as the effective radius for *S(Q)* when $P(Q) * S(Q)$ is applied.

For information about polarised and magnetic scattering, see
the :doc:`magnetic help <../sasgui/perspectives/fitting/mag_help>` documentation.

This code is based on the form factor calculations implemented in the NIST
Center for Neutron Research provided c-library (Kline, 2006).

References
----------

B Cabane, *Small Angle Scattering Methods*,
in *Surfactant Solutions: New Methods of Investigation*,
Ch.2, Surfactant Science Series Vol. 22, Ed. R Zana and M Dekker,
New York, (1987).

**Author:** NIST IGOR/DANSE **on:** pre 2010

**Last Modified by:** Piotr Rozyczko**on:** Feb 24, 2016

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
    n_pairs:number of pairs of water/shell
    background: incoherent background
        """
category = "shape:sphere"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["volfraction", "",  0.05, [0.0, 1],  "", "volume fraction of vesicles"],
    ["radius", "Ang", 60.0, [0.0, inf],  "", "Core radius of the multishell"],
    ["thick_shell", "Ang",        10.0, [0.0, inf],  "", "Shell thickness"],
    ["thick_solvent", "Ang",        10.0, [0.0, inf],  "", "Water thickness"],
    ["sld_solvent",    "1e-6/Ang^2",  6.4, [-inf, inf], "sld", "Core scattering length density"],
    ["sld",   "1e-6/Ang^2",  0.4, [-inf, inf], "sld", "Shell scattering length density"],
    ["n_pairs",     "",            2.0, [1.0, inf],  "", "Number of pairs of water and shell"],
    ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sph_j1c.c", "multilayer_vesicle.c"]

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
