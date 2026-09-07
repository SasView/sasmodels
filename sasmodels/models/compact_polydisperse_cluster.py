r"""
Definition
----------
Structure factor for a polydisperse spherical cluster with internal correlations.

WARNING: Should not be used in combination with very anisotropic particle shapes.

Structure factor for polydisperse compact spherical clusters with internal
correlations. The model is intended for clusters composed of spherical building
blocks and should not be used together with highly anisotropic form factors.

The cluster size distribution is described by a Schulz distribution. The model
is based on the treatment developed by Pedersen *et al.* (2022) and the
expressions presented in Pedersen *et al.* (2026). Modifications have been
introduced to allow a monodisperse nearest-neighbor distance within the cluster
and to use the weight-average aggregation number as a fitting parameter.

The model has four adjustable parameters,

.. math::

    R_L,\quad \sigma_L,\quad R_S,\quad N_{\mathrm{agg}}

where :math:`R_L` is the mean cluster radius (radius_cluster), :math:`\sigma_L` is the relative
polydispersity of the cluster radius (radius_cluster_sigma_relative), :math:`R_S` is the radius of the
elementary spheres (radius_effective), and

.. math::

    N_{\mathrm{agg}}
    = \frac{\langle N^2\rangle}{\langle N\rangle}

is the weight-average aggregation number (n_aggreg).

Let us first consider the monodisperse case. The elementary spheres are distributed with in a volume :math:`V(R_L-R_S)` in order to have the elementary spheres completely embedded. For a cluster containing :math:`N` elementary spheres,

.. math::

    N = \eta
    \frac{V(R_L-R_S)}{V_S},

where

.. math::

    V_S = \frac{4\pi R_S^3}{3}

is the volume of an elementary sphere and :math:`\eta` is the volume fraction
of spheres within the cluster.

Now, consider also polydispersity. The compact-cluster structure factor is then given by

.. math::

    S_{\mathrm{comp}}(q)
    = \tilde S_{\mathrm{HS}}(q)
    + \alpha
      \frac{\langle P_{R_L-R_S}(q)\rangle}
           {\langle V(R_L-R_S)\rangle},

where :math:`\langle P_{R_L-R_S}(q)\rangle` is the polydisperse form factor of a sphere for a Schulz distribution normaælized to :math:`\langle V(R_L-R_S)^2\rangle` at  :math:`q=0` and

.. math::

    \alpha = \frac{\eta}{V_S}.

Rather than fitting :math:`\eta` directly, it is calculated from the
aggregation number through

.. math::

    N_{\mathrm{agg}}
    =
    \alpha
    \frac{\langle V(R_L-R_S)^2\rangle}
         {\langle V(R_L-R_S)\rangle},

which yields

.. math::

    \eta
    =
    N_{\mathrm{agg}} V_S
    \frac{\langle V(R_L-R_S)\rangle}
         {\langle V(R_L-R_S)^2\rangle}.

For a Schulz distribution,

.. math::

    \eta
    =
    N_{\mathrm{agg}}
    \frac{R_S^3}{(R_L-R_S)^3}
    \frac{(Z_L+3)(Z_L+2)}
         {(Z_L+6)(Z_L+5)},

where

.. math::

    Z_L = \frac{1}{\sigma_L^2}-1.

The model employs an effective hard-sphere volume fraction

.. math::

    \eta_{\mathrm{eff}}
    = \eta\left(1-\frac{f}{2}\right),

with

.. math::

    f =
    \frac{(R_L-R_S)^3-(R_L-2R_S)^3}
         {(R_L-R_S)^3}.

Thus, :math:`R_L`, :math:`\sigma_L`, :math:`R_S`, and
:math:`N_{\mathrm{agg}}` are fitted parameters, whereas
:math:`\eta` and :math:`\eta_{\mathrm{eff}}` are calculated internally.



Parameters
----------
radius_effective : effective scatterer radius (Å), half the minimum center-to-center distance; first parameter for :math:`P@S` wiring (see sasmodels structure-factor conventions)
volfraction : unused in this :math:`S(q)`; required as second parameter for :math:`P@S` products
radius_cluster : mean radius of large clusters
radius_cluster_sigma_relative : relative polydispersity of radius_cluster
n_aggreg : weight-average aggregation number

References
----------

# Pedersen, J. S., Møller, T. L., & Corredig, M. (2026). Scattering from 'Babinet 'particles (or not…):
spherical particles made up of spheres and spherical particles with s
pherical voids. Applied Crystallography, 59(1).

# Pedersen, J. S., Møller, T. L., Raak, N., & Corredig,
M. (2022). A model on an absolute scale for the small-angle X-ray scattering
from bovine casein micelles. Soft Matter, 18(45), 8613-8625.

# Larsen, A. H., Pedersen, J. S., & Arleth, L. (2020). Assessment of
structure factors for analysis of small-angle scattering data from
desired or undesired aggregates. Applied Crystallography, 53(4), 991-1005.

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
single = False

# Must match C: Iq(Q, radius_effective, volfraction, radius_cluster, radius_cluster_sigma_relative, n_aggreg)
parameters = [
    ["radius_effective", "Ang", 10.0, [0.0, np.inf], "",
     "effective scatterer radius (half center-to-center distance)"],
    ["volfraction", "", 0.2, [0.0, 1.0], "",
     "unused in S(q); required for P@S products"],
    ["radius_cluster", "Ang", 40.0, [0.0, np.inf], "", "Average cluster radius"],
    ["radius_cluster_sigma_relative", "", 0.4, [0.0, 1.0], "", "Relative size polydispersity"],
    ["n_aggreg", "", 50.0, [10.0, 100.0], "", "Weight-average aggregation number"],
]

source = ["lib/sas_gammainc.c", "lib/sas_3j1x_x.c", "compact_polydisperse_cluster.c"]

def random():
    import random

    return {
        "radius_effective": random.uniform(2.5, 15.0),
        "volfraction": random.uniform(0.01, 0.3),
        "radius_cluster": random.uniform(50, 200),
        "radius_cluster_sigma_relative": random.uniform(0.01, 0.3),
        "n_aggreg": random.uniform(10, 100),
    }

def test():
    print("compact_polydisperse_cluster model loaded successfully.")
