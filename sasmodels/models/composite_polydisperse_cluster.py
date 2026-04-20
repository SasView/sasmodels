r"""
Definition
----------
Structure factor for a polydisperse spherical cluster with internal correlations..

WARNING: Should not be used in combination with very anisotrpic particle shapes.

Polydispersity described by Schulz distribution
and expressions taken from Pedersen, J. S., Møller, T. L., Raak, N., & Corredig,
M. (2022). A model on an absolute scale for the small-angle X-ray scattering
from bovine casein micelles. Soft Matter, 18(45), 8613-8625.

The implemented expressions are described in
Pedersen, J. S., Møller, T. L., & Corredig, M. (2026). Scattering from 'Babinet 'particles (or not…):
spherical particles made up of spheres and spherical particles with s
pherical voids. Applied Crystallography, 59(1).

Modifications have been done for adapting to monodisperse distance between
points in cluster  and for allowing weight averrage aggregation number to  be a fit paramter.

See also
Larsen, A. H., Pedersen, J. S., & Arleth, L. (2020). Assessment of
structure factors for analysis of small-angle scattering data from
desired or undesired aggregates. Applied Crystallography, 53(4), 991-1005.


Parameters:
R_clust   : mean radius of large clusters
sig_rel_R : relative polydispersity of R_cluster
dist_points   : minimum distance between scatterers
N_agg  : weight average aggregation number

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

name = "composite_polydisperse_cluster"
title = "Composite Polydisperse Cluster Structure Factor"
description = """
Structure factor for a polydisperse spherical polydisperse cluster
with internal correlations.
"""

category = "structure-factor"

# Must match C kernel signature: Iq(Q, RL, SIGL, RS, ETA)
parameters = [
    ["R_clust",   "Ang", 40.0, [0.0, np.inf], "", "Average cluster radius"],
    ["sig_rel_R", "",      0.4, [0.0, 1.0],    "", "Relative size polydispersity"],
    ["dist_points",   "Ang",  20.0, [0.0, np.inf], "", "Minimum distance between scatterers"],
    ["N_agg",  "",      50, [10, 100],   "", "Weight average aggregation number"]
]

source = ["composite_polydisperse_cluster.c"]

def random():
    import random
    return {
        "R_clust":   random.uniform(50, 200),
        "sig_rel_R": random.uniform(0.01, 0.3),
        "dist_points":   random.uniform(5, 30),
        "N_agg":  random.uniform(10, 1000)
    }

def test():
    print("composite_polydisperse_cluster model loaded successfully.")
