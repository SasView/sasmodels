r"""
The Benoit model for a simple star polymer, with Gaussian coils arms from 
a common point.

Definition
----------

For a star with $f$ arms the scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = \frac{2}{fv^2}\left[ v-1+exp(-v)+\frac{f-1}{2}
           \left[ 1-exp(-v)\right]^2\right]

where

.. math::

    v=\frac{u^2f}{(3f-2)}

and

.. math::

    u = \left\langle R_{g}^2\right\rangle q^2

contains the square of the ensemble average radius-of-gyration of an arm.  Note that 
when there is only one arm, $f$ = 1, the Debye Gaussian coil equation is recovered.
Star polymers in solutions tend to have strong interparticle and osmotic effects, so 
the Benoit equation may not work well. At small q the Guinier term and hence I(q=0)
is the same as for $f$ arms of radius of gyration $R_g$, as described for the :ref:`mono_gauss_coil <mono-gauss-coil>`.


References
----------

H Benoit *J. Polymer Science*, 11, 596-599 (1953)


"""

from numpy import inf

name = "star_polymer"
title = "Star polymer model with Gaussian statistics"
description = """
        Benoit 'Star polymer with Gaussian statistics'
        with
        P(q) = 2/{fv^2} * (v - (1-exp(-v)) + {f-1}/2 * (1-exp(-v))^2)
        where
        - v = u^2f/(3f-2)
        - u = <R_g^2>q^2, where <R_g^2> is the ensemble average radius of
        gyration squared of an arm
        - f is the number of arms on the star
        """
category = "shape-independent"
single = False
# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["rg_squared", "Ang^2", 100.0, [0.0, inf], "", "Ensemble radius of gyration SQUARED of an arm"],
              ["arms",    "",      3,   [1.0, 6.0], "", "Number of arms in the model"],
             ]
# pylint: enable=bad-whitespace, line-too-long

source = ["star_polymer.c"]

demo = dict(scale=1, background=0,
            rg_squared=100.0,
            arms=3.0)

oldname = 'StarPolymer'

oldpars = dict(rg_squared='R2',
               arms='arms')

tests = [[{'rg_squared': 2.0,
           'arms':    3.3,
          }, 0.5, 0.851646091108],

         [{'rg_squared':    1.0,
           'arms':       2.0,
           'background': 1.8,
          }, 1.0, 2.53575888234],
        ]
# 23Mar2016  RKH edited docs, would this better use rg not rg^2 ? Numerical noise at extremely small q.rg