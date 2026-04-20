r"""
Definition
----------

Linear aggregate structure factor.

WARNING: Should not be used in combination with very anisotrpic particle shapes.

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

parameters = [
    ['N_agg', '', 50.0, [0, np.inf], '', 'Aggregation number'],
    ['dist_points', 'Ang', 20.0, [0, np.inf], '', 'distance between scatterers']
]

source = ["linear_aggregate.c"]

def random():
    import random
    return {
        'N_agg': random.uniform(1, 100),
        'dist_points': random.uniform(10, 50)
    }

def test():
    print('Steucture factor of linear aggregate')
