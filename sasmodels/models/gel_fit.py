r"""
*This model was implemented by an interested user!*

Unlike a concentrated polymer solution, the fine-scale polymer distribution
in a gel involves at least two characteristic length scales,
a shorter correlation length ( $a1$ ) to describe the rapid fluctuations
in the position of the polymer chains that ensure thermodynamic equilibrium,
and a longer distance (denoted here as $a2$ ) needed to account for the static
accumulations of polymer pinned down by junction points or clusters of such
points. The latter is derived from a simple Guinier function. Compare also the 
gauss_lorentz_gel model.


Definition
----------

The scattered intensity $I(q)$ is calculated as

.. math::

    I(Q) = I(0)_L \frac{1}{\left( 1+\left[ ((D+1/3)Q^2a_{1}^2
    \right]\right)^{D/2}} + I(0)_G exp\left( -Q^2a_{2}^2\right) + B

where

.. math::

    a_{2}^2 \approx \frac{R_{g}^2}{3}

Note that the first term reduces to the Ornstein-Zernicke equation
when $D = 2$; ie, when the Flory exponent is 0.5 (theta conditions).
In gels with significant hydrogen bonding $D$ has been reported to be
~2.6 to 2.8.


References
----------

Mitsuhiro Shibayama, Toyoichi Tanaka, Charles C Han,
*J. Chem. Phys.* 1992, 97 (9), 6829-6841

Simon Mallam, Ferenc Horkay, Anne-Marie Hecht, Adrian R Rennie, Erik Geissler,
*Macromolecules* 1991, 24, 543-548

"""

from numpy import inf

name = "gel_fit"
title = "Fitting using fine-scale polymer distribution in a gel."
description = """\
    Structure factor for interacting particles:

    Shibayama-Geissler Two-Length Scale Fit for Gels (GelFit)

    Shibayama; Tanaka; Han J Chem Phys (1992), 97(9), 6829-6841
    Mallam; Horkay; Hecht; Rennie; Geissler, Macromol (1991), 24, 543
"""
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["guinier_scale",    "cm^-1",   1.7, [-inf, inf], "", "Guinier length scale"],
              ["lorentzian_scale", "cm^-1",   3.5, [-inf, inf], "", "Lorentzian length scale"],
              ["gyration_radius",  "Ang",     104.0, [2, inf],    "", "Radius of gyration"],
              ["fractal_exp",      "",          2.0, [0, inf],    "", "Fractal exponent"],
              ["cor_length",       "Ang",      16.0, [0, inf],    "", "Correlation length"]
             ]
# pylint: enable=bad-whitespace, line-too-long

source = ["gel_fit.c"]

demo = dict(background=0.01,
            guinier_scale=1.7,
            lorentzian_scale=3.5,
            gyration_radius=104,
            fractal_exp=2.0,
            cor_length=16.0)

tests = [[{'guinier_scale': 1.0,
           'lorentzian_scale': 1.0,
           'gyration_radius': 10.0,
           'fractal_exp': 10.0,
           'cor_length': 20.0,
           'background': 0.0,
          }, 0.1, 0.716532],

         [{'guinier_scale': 4.0,
           'lorentzian_scale': 10.0,
           'gyration_radius': 500.0,
           'fractal_exp': 1.0,
           'cor_length': 20.0,
           'background': 20.0,
          }, 5.0, 20.1224653026],
        ]
