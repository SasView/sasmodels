r"""
Definition
----------
Structure factor for a polydisperse spherical cluster with internal correlations.

WARNING: Should not be used in combination with very anisotropic particle shapes.

Polydispersity described by Schulz distribution
and expressions taken from Pedersen, J. S., Møller, T. L., Raak, N., & Corredig,
M. (2022). A model on an absolute scale for the small-angle X-ray scattering
from bovine casein micelles. Soft Matter, 18(45), 8613-8625.

The implemented expressions are described in
Pedersen, J. S., Møller, T. L., & Corredig, M. (2026). Scattering from 'Babinet 'particles (or not…):
spherical particles made up of spheres and spherical particles with s
pherical voids. Applied Crystallography, 59(1).

Modifications have been done for adapting to monodisperse distance between
points in cluster and for allowing weight-average aggregation number to be a fit parameter.

See also
Larsen, A. H., Pedersen, J. S., & Arleth, L. (2020). Assessment of
structure factors for analysis of small-angle scattering data from
desired or undesired aggregates. Applied Crystallography, 53(4), 991-1005.


Parameters
----------
radius_effective : effective scatterer radius (Å), half the minimum center-to-center distance; first parameter for :math:`P@S` wiring (see sasmodels structure-factor conventions)
volfraction : unused in this :math:`S(q)`; required as second parameter for :math:`P@S` products
R_clust : mean radius of large clusters
sig_rel_R : relative polydispersity of R_cluster
N_agg : weight-average aggregation number

References
----------

# Pedersen, J. S., Møller, T. L., & Corredig, M. (2026). Scattering from 'Babinet 'particles (or not…):
spherical particles made up of spheres and spherical particles with s
pherical voids. Applied Crystallography, 59(1).

Authorship and Verification
----------------------------

* **Author:** Jan Skov Pedersen
* **Last Modified by:** Jan Skov Pedersen, April 12, 2026
* **Last Reviewed by:**  **Date:**
"""


import numpy as np

name = "compact_polydisperse_cluster"
title = "Compact Polydisperse Cluster Structure Factor"
description = """
Structure factor for a compact polydisperse cluster
with internal correlations.
"""

category = "structure-factor"
structure_factor = True

# Must match C: Iq(Q, radius_effective, volfraction, R_clust, sig_rel_R, N_agg)
parameters = [
    ["radius_effective", "Ang", 10.0, [0.0, np.inf], "",
     "effective scatterer radius (half center-to-center distance)"],
    ["volfraction", "", 0.2, [0.0, 1.0], "",
     "unused in S(q); required for P@S products"],
    ["R_clust", "Ang", 40.0, [0.0, np.inf], "", "Average cluster radius"],
    ["sig_rel_R", "", 0.4, [0.0, 1.0], "", "Relative size polydispersity"],
    ["N_agg", "", 50.0, [10.0, 100.0], "", "Weight-average aggregation number"],
]

source = ["compact_polydisperse_cluster.c"]

def random():
    import random

    return {
        "radius_effective": random.uniform(2.5, 15.0),
        "volfraction": random.uniform(0.01, 0.3),
        "R_clust": random.uniform(50, 200),
        "sig_rel_R": random.uniform(0.01, 0.3),
        "N_agg": random.uniform(10, 100),
    }

def test():
    print("compact_polydisperse_cluster model loaded successfully.")
