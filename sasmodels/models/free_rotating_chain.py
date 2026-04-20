r"""
Definition
----------

Free‑rotating chain structure factor: Eandom flight.

WARNING: Should not be used in combination with very anisotrpic particle shapes. 

This model calculates the structure factor of a free‑rotating chain
of scattering points, following the treatment of Burchard and Kajiwara.

Equations are given in 
Larsen, A. H., Pedersen, J. S., & Arleth, L. (2020). Assessment of 
structure factors for analysis of small-angle scattering data from 
desired or undesired aggregates. Applied Crystallography, 53(4), 991-1005.

The model is more specifically eq. (23,24) in the paper.

A non‑integer effective number of rotating points is handled by linear
interpolation between integer chain lengths.

The model evaluates:

.. math::

    S(Q) = (1-w)\,S_N(Q) + w\,S_{N+1}(Q)

where:

- :math:`N = \lfloor RN \rfloor`
- :math:`w = RN - N`

and :math:`S_N(Q)` is the Debye sum for a free‑rotating chain.

References
----------

Burchard, W., & Kajiwara, K. (1970).
The statistics of stiff chain molecules I.
The particle scattering factor.
*Proceedings of the Royal Society of London A* **316**, 185–199.

Authorship and Verification
----------------------------

* **Author:** Jan Skov Pedersen
* **Converted to sasmodels kernel by:** —
* **Last Reviewed by:** — **Date:** April 10, 2026
"""

import numpy as np

# Model identification
name = "free_rotating_chain"
title = "Free-rotating chain structure factor"
description = """
Structure factor of a free-rotating chain (random flight)
"""

category = "structure-factor"

# Parameters must match the kernel signature: Iq(Q, RN, d)
parameters = [
    ["N_agg", "", 20.0, [1.0, np.inf], "", "Aggregation number"],
    ["dist_points", "Ang", 20.0, [0.0, np.inf], "", "Separation between points"],
]

# Kernel source
source = ["free_rotating_chain.c"]

def random():
    import random
    return {
        "N_agg": random.uniform(1.0, 100.0),
        "dist_points": random.uniform(10.0, 50.0),
    }

def test():
    print("Structure factor of free‑rotating chain (free_rotating_chain)")