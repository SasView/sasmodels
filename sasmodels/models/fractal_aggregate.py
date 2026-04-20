# fractal_aggregate.py
# Sasmodels plugin wrapper for the fractal structure factor
#
# Matches the format and syntax of free_rotating_chain.py
r"""
Definition
----------
This model calculates the structure factor of a fractal-like aggregates.

WARNING: Should not be used in combination with very anisotrpic particle shapes.

It uses the following equation:

.. math::

       S(q) &= 1 + \frac{D_f\  \Gamma\!(D_f-1)}{[1+1/(q \xi)^2\  ]^{(D_f -1)/2}}
    \frac{\sin[(D_f-1) \tan^{-1}(q \xi) ]}{(q R_0)^{D_f}}

where $\xi$ is the correlation length representing the cluster size and $D_f$
is the fractal dimension, representing the self similarity of the structure.

The expresion has been reformulated so that the aggregation number $N_{agg}$ is a fit paramter in stead of $\xi$
Parameters
----------
d : distance between adjacent scatterers
D_fract : fractal dimension
N_agg : number of particles in the fractal cluster

See also
Larsen, A. H., Pedersen, J. S., & Arleth, L. (2020). Assessment of
structure factors for analysis of small-angle scattering data from
desired or undesired aggregates. Applied Crystallography, 53(4), 991-1005.



Validation
----------

Translated FORTRAN code

References
----------

# Teixeira, J. (1988). Small-angle scattering by fractal systems. Applied Crystallography, 21(6), 781-785.

Authorship and Verification
----------------------------

* **Author:** Jan Skov Pedersen
* **Last Modified by:** Jan Skov Pedersen, April 12, 2026
* **Last Reviewed by:** Reviewer Name Here **Date:**
"""



import numpy as np

name = "fractal_aggregate"
title = "Fractal aggregate structure factor"
description = """
Fractal structure factor S(q)
"""

category = "structure-factor"

# Must match C kernel signature: Iq(q, dist_points, D_fract, N_agg)
parameters = [
    ["dist_points",  "Ang", 20.0, [0.0, np.inf], "", "distance between scatterers"],
    ["D_fract",  "",     2.0, [1.0, 3.0],     "", "Fractal dimension"],
    ["N_agg",  "",    50.0, [1.0, np.inf], "", "Number of particles"],
]

# Kernel source
source = ["fractal_aggregate.c"]

def random():
    import random
    return {
        "dist_popints": random.uniform(5.0, 50.0),
        "D_fract": random.uniform(1.1, 2.9),
        "N_agg": random.uniform(5.0, 300.0),
    }

def test():
    print("fractal_aggregate plugin loaded correctly.")
