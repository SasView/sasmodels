r"""
Fractal Structure Factor S(q)
=============================

Fractal structure factor :math:`S(q)` using a discrete-chain high-q baseline
and a low-q crossover.

WARNING: Should not be used in combination with very anisotrpic particle shapes. 

This model calculates the structure factor of fractal-like aggregates
according to the following equation:

.. math::

   S(q) = 1 +
   \frac{D_f\,\Gamma(D_f - 1)}
        {\left[1 + 1/(q \xi)^2\right]^{(D_f - 1)/2}}
   \frac{\sin\!\left[(D_f - 1)\tan^{-1}(q \xi)\right]}
        {(q R_0)^{D_f}}

Here, :math:`\xi` is the correlation length representing the cluster size,
and :math:`D_f` is the fractal dimension characterizing the internal
self-similarity of the structure.

The expression has been reformulated so that the aggregation number
:math:`N_{\mathrm{agg}}` is used as a fit parameter instead of :math:`\xi`.

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

See also
Larsen, A. H., Pedersen, J. S., & Arleth, L. (2020). Assessment of 
structure factors for analysis of small-angle scattering data from 
desired or undesired aggregates. Applied Crystallography, 53(4), 991-1005.


Parameters
----------

``d``  
    Distance between adjacent scatterers.

``D_fract``  
    Fractal dimension :math:`D_f`.

``N_agg``  
    Number of particles in the fractal cluster.

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

# fractal_aggregate_discrete_chain.py
# Sasmodels plugin wrapper for fractal S(q) using inter-scatterer distance d.

import numpy as np

name = "fractal_aggregate_discrete_chain"
title = "Fractal aggregate with discrete-chain baseline"
description = """
Fractal Structure Factor S(q) with local correlation between points
"""

category = "structure-factor"

# Must match C signature: Iq(q, d, D, N)
parameters = [
    ["dist_points", "Ang", 20.0, [0.0, np.inf], "", "Inter-scatterer distance"],
    ["D_fract", "",      2.0, [1.0, 3.0],     "", "Fractal dimension"],
    ["N_agg", "",     50.0, [1.0, np.inf], "", "Number of particles in cluster"]
]

source = ["fractal_aggregate_discrete_chain.c"]

def random():
    import random
    return {
        "d_points": random.uniform(10.0, 80.0),
        "D_fract": random.uniform(1.1, 2.9),
        "N_agg": random.uniform(5.0, 300.0)
    }

def test():
    print("fractal_aggregate_discrete_chain plugin loaded correctly.")