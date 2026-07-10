r"""
Definition
----------

Linear aggregate structure factor.

WARNING: Should not be used in combination with very anisotropic particle shapes.

Equations are given in
Larsen, A. H., Pedersen, J. S., & Arleth, L. (2020). Assessment of
structure factors for analysis of small-angle scattering data from
desired or undesired aggregates. Applied Crystallography, 53(4), 991-1005.

The model is more specifically eq. (21,22) in the paper.

A non‑integer effective number of points is handled by linear
interpolation between integer chain lengths.

The model evaluates:

.. math::

    S(Q) = (1-w)\,S_N(Q) + w\,S_{N+1}(Q)

where:

- :math:`N = \lfloor RN \rfloor`
- :math:`w = RN - N`

and :math:`S_N(Q)` is the Debye sum for the aggregate.

Parameters
----------
radius_effective : effective scatterer radius (Å), half the center-to-center distance between scatterers on the chain.
volfraction : unused in this :math:`S(q)`; required for :math:`P@S` products.
n_aggreg : aggregation number.

Validation
----------

Translated FORTRAN code

References
----------

# Larsen, A. H., Pedersen, J. S., & Arleth, L. (2020). Assessment of
structure factors for analysis of small-angle scattering data from
desired or undesired aggregates. Applied Crystallography, 53(4), 991-1005.

Authorship and Verification
----------------------------

* **Author:** Jan Skov Pedersen
* **Last Modified by:** Jan Skov Pedersen, April 12, 2026
* **Last Reviewed by:** Reviewer Name Here **Date:**
"""

import numpy as np

name = "linear_aggregate"
title = "Linear aggregate structure factor"
description = """\
Linear aggregate structure factor
"""
category = "structure-factor"
structure_factor = True

parameters = [
    ["radius_effective", "Ang", 10.0, [0.0, np.inf], "",
     "effective scatterer radius (half center-to-center distance)"],
    ["volfraction", "", 0.2, [0.0, 1.0], "",
     "unused in S(q); required for P@S products"],
    ["n_aggreg", "", 50.0, [0.0, 100.0], "", "Aggregation number"],
]

valid = "n_aggreg <= 100"

source = ["linear_aggregate.c"]

def random():
    import random

    return {
        "radius_effective": random.uniform(5, 25),
        "volfraction": random.uniform(0.01, 0.3),
        "n_aggreg": random.uniform(1, 100),
    }

def test():
    print("Structure factor of linear aggregate")
