r"""
Fractal Structure Factor S(q)
=============================

Structure factor for fractal-like aggregates composed of spherical primary
particles, using a discrete-chain high-q baseline and a low-q crossover.

WARNING: This model assumes isotropic aggregates and should not be used in
combination with strongly anisotropic particle shapes.

The model is based on the mass-fractal structure factor originally derived by
Teixeira (1988) and describes aggregates characterized by a fractal
dimension :math:`D_f` and a correlation length :math:`\xi`. The fractal
dimension governs the compactness of the aggregate, with larger values
corresponding to denser structures, while :math:`\xi` represents the aggregate
size.

The original expression is

.. math::

    S(q)
    =
    1
    +
    \frac{D_f\Gamma(D_f-1)}
         {\left[1+(q\xi)^{-2}\right]^{(D_f-1)/2}}
    \frac{\sin\left[(D_f-1)\tan^{-1}(q\xi)\right]}
         {(qR_0)^{D_f}},

where :math:`R_0` is the radius of the primary particle.

To allow better comparisons with the other aggregate structure factors, the model is reparameterized so that the
weight-average aggregation number :math:`N_{\mathrm{agg}}`
becomes a fitting parameter instead of the correlation length
:math:`\xi`.

Assuming the aggregation number is related to the cutoff length by

.. math::

    N
    =
    \left(\frac{\xi}{R_0}\right)^{D_f},

the correlation length may be written as

.. math::

    \xi
    =
    R_0 N^{1/D_f}.

Defining

.. math::

    y
    =
    qR_0N^{1/D_f},

the structure factor may be expressed as

.. math::

    S(q)
    =
    1+A\,K(q),

with

.. math::

    K(q)
    =
    \sin[(D_f-1)\arctan(y)]
    (qR_0)^{-D_f}
    \left(1+\frac{1}{y^2}\right)^{-(D_f-1)/2}.

In the limit :math:`q\rightarrow0`,

.. math::

    K(q)
    \rightarrow
    (D_f-1)N.

The normalization constant :math:`A` is chosen so that

.. math::

    S(0)=N,

yielding

.. math::

    A
    =
    \frac{N-1}
         {N(D_f-1)}.

The final expression implemented in the model is therefore

.. math::

    S(q)
    =
    1
    +
    \frac{N-1}{N(D_f-1)}
    \sin[(D_f-1)\arctan(qR_0N^{1/D_f})]
    (qR_0)^{-D_f}
    \left(
        1+\frac{1}
                 {(qR_0N^{1/D_f})^2}
    \right)^{-(D_f-1)/2},

where :math:`N=N_{\mathrm{agg}}`.

The model therefore fits the physically intuitive aggregation number directly
while preserving the low- :math:`q` limit

.. math::

    S(0)=N_{\mathrm{agg}}.

When combined with a polydisperse form factor, the simple
:math:`P(q)S(q)` approximation evaluated at an effective particle size is not
generally equivalent to averaging :math:`P(q,r)S(q,r)` over the particle-size
distribution. For systems with substantial primary-particle polydispersity,
the :math:`\beta` approximation (*structure_factor_mode = 1*) is often more
appropriate.

High-q Baseline Modification
----------------------------

The high-q limit is modified to simulate local point correlations, as in a
random-flight model.
The term ``1`` in the expression for :math:`S(q)` is replaced by:

.. math::

   S(q) = 1 +
   ( 2\frac{\sin(q d)}{q d}
   + 2(\frac{\sin(q d)}{q d})^2 )
   \left( 1 - e^{-(q d / C)^4} \right)

The low-q crossover constant is :math:`C = 5.2`.

Parameters
----------
radius_effective : effective scatterer radius (Å), half the center-to-center distance :math:`d`; see *radius_effective_mode* when combining with a form factor
volfraction : unused in this :math:`S(q)`; required for :math:`P@S` products
fractal_dim : fractal dimension :math:`D_f`
n_aggreg : number of particles in the fractal cluster

See also
Larsen, A. H., Pedersen, J. S., & Arleth, L. (2020). Assessment of
structure factors for analysis of small-angle scattering data from
desired or undesired aggregates. Applied Crystallography, 53(4), 991-1005.


Validation
----------

Translated from original FORTRAN code.

References
----------

* Teixeira, J. (1988). *Small-angle scattering by fractal systems*.
  **Journal of Applied Crystallography**, 21(6), 781–785.

Authorship and Verification
---------------------------

* **Author:** Jan Skov Pedersen
* **Last Modified by:** Jan Skov Pedersen, April 12, 2026
* **Last Reviewed by:** Reviewer Name Here (Date)
"""

# Sasmodels plugin wrapper for fractal S(q) with discrete-chain baseline.

import numpy as np

name = "fractal_aggregate_discrete_chain"
title = "Fractal aggregate with discrete-chain baseline"
description = """
Fractal Structure Factor S(q) with local correlation between points
"""

category = "structure-factor"
structure_factor = True

# Must match C: Iq(q, radius_effective, volfraction, fractal_dim, n_aggreg)
parameters = [
    ["radius_effective", "Ang", 10.0, [0.0, np.inf], "",
     "effective scatterer radius (half center-to-center distance)"],
    ["volfraction", "", 0.2, [0.0, 1.0], "",
     "unused in S(q); required for P@S products"],
    ["fractal_dim", "", 2.0, [1.0, 3.0], "", "Fractal dimension"],
    ["n_aggreg", "", 50.0, [1.0, np.inf], "", "Number of particles in cluster"],
]

source = ["fractal_aggregate_discrete_chain.c"]

def random():
    import random

    return {
        "radius_effective": random.uniform(5.0, 40.0),
        "volfraction": random.uniform(0.01, 0.3),
        "fractal_dim": random.uniform(1.1, 2.9),
        "n_aggreg": random.uniform(5.0, 300.0),
    }

def test():
    print("fractal_aggregate_discrete_chain plugin loaded correctly.")
