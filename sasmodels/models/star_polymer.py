r"""
Calculates the scattering from fractal-like aggregates based on
the Mildner reference. This model is also known as the Benoit Star model.

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

is the square of the ensemble average radius-of-gyration of an arm.

.. figure:: img/star_polymer_1d.jpg

    1D plot using the default values.


Reference
---------

H Benoit *J. Polymer Science*, 11, 596-599 (1953)


"""

from numpy import inf

name = "star_polymer"
title = "Star polymer model with Gaussian statistics"
description = """
        Scattering model class for 'Star polymer with Gaussian statistics'
        with
        P(q) = 2/{fv^2} * (v - (1-exp(-v)) + {f-1}/2 * (1-exp(-v))^2)
        where
        - v = u^2f/(3f-2)
        - u = <R_g^2>q^2, where <R_g^2> is the ensemble average radius of
        giration squared of an arm
        - f is the number of arms on the star
        """
category = "shape-independent"
single = False
# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["radius2", "Ang", 100.0, [0.0, inf], "", "Ensemble radius of gyration squared of an arm"],
              ["arms",    "",      3,   [1.0, 6.0], "", "Number of arms in the model"],
             ]
# pylint: enable=bad-whitespace, line-too-long

source = ["star_polymer.c"]

demo = dict(scale=1, background=0,
            radius2=100.0,
            arms=3.0)

oldname = 'StarPolymer'

oldpars = dict(radius2='R2',
               arms='arms')

tests = [[{'radius2': 2.0,
           'arms':    3.3,
          }, 0.5, 0.850646091108],

         [{'radius2':    1.0,
           'arms':       2.0,
           'background': 1.8,
          }, 1.0, 2.53575888234],
        ]
