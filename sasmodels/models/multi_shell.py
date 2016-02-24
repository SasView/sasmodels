r"""
This model provides the form factor, *P(q)*, for a multi-lamellar vesicle
with *N* shells where the core is filled with solvent and the shells are
interleaved with layers of solvent. For *N = 1*, this returns the VesicleModel.

Definition
----------

.. figure:: img/multi_shell_fig1.jpg

The 2D scattering intensity is the same as 1D, regardless of the orientation
of the q vector which is defined as:

.. math::

    q = \sqrt{q_x^2 + q_y^2}

.. note:
    The outer most radius
    $core_radius + n_pairs * s_thickness + (n_pairs - 1) * w_thickness$
    is used as the effective radius for *S(Q)* when $P(Q) * S(Q)$ is applied.


.. figure:: img/multi_shell_1d.jpg

    1D plot using the default values (with 200 data point).

Our model uses the form factor calculations implemented in a c-library provided
by the NIST Center for Neutron Research (Kline, 2006).

Reference
---------
B Cabane, *Small Angle Scattering Methods*,
in *Surfactant Solutions: New Methods of Investigation*,
Ch.2, Surfactant Science Series Vol. 22, Ed. R Zana and M Dekker,
New York, (1987).

"""

from numpy import inf

name = "multi_shell"
title = "Multi shell model"
description = """
    MultiShell (Sphere) Model (or Multilamellar Vesicles): Model parameters;
    scale : scale factor
    core_radius : Core radius of the multishell
    s_thickness: shell thickness
    w_thickness: water thickness
    core_sld: core scattering length density
    shell_sld: shell scattering length density
    n_pairs:number of pairs of water/shell
    background: incoherent background
        """
category = "shape:sphere"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["core_radius", "Ang",        60.0, [0.0, inf],  "", "Core radius of the multishell"],
    ["s_thickness", "Ang",        10.0, [0.0, inf],  "", "Shell thickness"],
    ["w_thickness", "Ang",        10.0, [0.0, inf],  "", "Water thickness"],
    ["core_sld",    "1e-6/Ang^2",  6.4, [-inf, inf], "", "Core scattering length density"],
    ["shell_sld",   "1e-6/Ang^2",  0.4, [-inf, inf], "", "Shell scattering length density"],
    ["n_pairs",     "",            2.0, [1.0, inf],  "", "Number of pairs of water and shell"],
    ]
# pylint: enable=bad-whitespace, line-too-long

source = ["multi_shell.c"]

demo = dict(scale=1, background=0,
            core_radius=60.0,
            s_thickness=10.0,
            w_thickness=10.0,
            core_sld=6.4,
            shell_sld=0.4,
            n_pairs=2.0)

oldname = 'MultiShellModel'
oldpars = dict()

tests = [
    # Accuracy tests based on content in test/utest_other_models.py
    [{'core_radius': 60.0,
      's_thickness': 10.0,
      'w_thickness': 10.0,
      'core_sld': 6.4,
      'shell_sld': 0.4,
      'n_pairs': 2.0,
      'scale': 1.0,
      'background': 0.001,
     }, 0.001, 2442.81],

    [{'core_radius': 60.0,
      's_thickness': 10.0,
      'w_thickness': 10.0,
      'core_sld': 6.4,
      'shell_sld': 0.4,
      'n_pairs': 2.0,
      'scale': 1.0,
      'background': 0.001,
     }, (0.001, 0.30903), 1.61873],
    ]
