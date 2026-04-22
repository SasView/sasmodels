# fractal_aggregate.py
# Sasmodels plugin wrapper for the fractal structure factor
#
# Matches the format and syntax of free_rotating_chain.py
r"""
Definition
----------
This model calculates the structure factor of a fractal-like aggregates.

WARNING: Should not be used in combination with very anisotropic particle shapes.

It uses the following equation:

.. math::

    S(q) = 1 + \frac{D_f\, \Gamma(D_f-1)}{\left[1+1/(q\xi)^2\right]^{(D_f-1)/2}}
    \frac{\sin\left[(D_f-1)\tan^{-1}(q\xi)\right]}{(q R_0)^{D_f}}

where $\xi$ is the correlation length representing the cluster size and $D_f$
is the fractal dimension, representing the self similarity of the structure.

The expression has been reformulated so that the aggregation number $N_{agg}$
is a fit parameter instead of $\xi$.

When combined with a polydisperse form factor, simple $P(q)\,S(q)$ with
$S(q)$ evaluated at an effective radius differs from averaging
$P(q,r)\,S(q,r)$ over $r$. For polydisperse primary particles, the
$\beta$-/decoupling approximation (*structure_factor_mode = 1*) is usually
more appropriate.

Parameters
----------
radius_effective : distance between adjacent scatterers (Å); use *unconstrained*
    *radius_effective_mode* to fit independently of the form-factor radius.
volfraction : required for $P@S$ products; **not used** in this $S(q)$
    (see ``sasmodels.product`` until optional structure-factor volfraction exists).
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
structure_factor = True

# Must match C kernel: Iq(q, radius_effective, volfraction, D_fract, N_agg)
parameters = [
    ["radius_effective", "Ang", 20.0, [0.0, np.inf], "",
     "distance between adjacent scatterers"],
    ["volfraction", "", 0.2, [0.0, 1.0], "",
     "unused in S(q); required as second parameter for P@S products"],
    ["D_fract", "", 2.0, [1.0, 3.0], "", "Fractal dimension"],
    ["N_agg", "", 50.0, [1.0, np.inf], "", "Number of particles"],
]

# Kernel source
source = ["fractal_aggregate.c"]

def random():
    import random

    return {
        "radius_effective": random.uniform(5.0, 50.0),
        "volfraction": random.uniform(0.01, 0.3),
        "D_fract": random.uniform(1.1, 2.9),
        "N_agg": random.uniform(5.0, 300.0),
    }

def test():
    print("fractal_aggregate plugin loaded correctly.")
