r"""
This model provides the form factor P(q) for nanoprims with different cross-sections.
It contains both the form factor for a reference orientation and the 1D form factor after orientation average
using the Fibonacci quadrature. This quadrature provides a quasi-uniform distribution of points on the unit sphere
using the golden ratio. The only new parameter is the number of points to generate on the unit sphere.
The default value (500 points) provides a good balance between accuracy and computational efficiency.

Definition
----------

expliquer les paramètres, donner les formules etc (voir article)


.. math::

    formule

.. figure:: img/nanoprisms_intensity_plot.png

Validation
----------

Validation of the code is made using numerical checks.
Comparisons with Debye formula calculations were made using DebyeCalculator library (https://github.com/FrederikLizakJohansen/DebyeCalculator).
Good agreement was found at q < 0.1 1/Angstrom.

References
----------

1. Jules Marcone, Jaime Gabriel Trazo, Rahul Nag, Claire Goldmann, Nicolas Ratel-Ramond, Cyrille Hamon and Marianne Imperor-Clerc (2025) J. Marcone et al. ,
 Form factor of prismatic particles for small-angle scattering analysis, JAC 2025

1. Wei-Ren Chen et al. "Scattering functions of Platonic solids".
   In: Journal of Applied Crystallography - J APPL CRYST 44 (June 2011).
   DOI: 10.1107/S0021889811011691

2. Croset, Bernard, "Form factor of any polyhedron: a general compact
   formula and its singularities" In: J. Appl. Cryst. (2017). 50, 1245–1255
   https://doi.org/10.1107/S1600576717010147

3. Wuttke, J. Numerically stable form factor of any polygon and polyhedron
   J Appl Cryst 54, 580-587 (2021)
   https://doi.org/10.1107/S160057672100171

Authorship and Verification
----------------------------

* **Authors:** Marianne Imperor-Clerc (marianne.imperor@cnrs.fr)
             Helen Ibrahim (helenibrahim1@outlook.com) ?
             Jules Marcone ()
             Sara Mokhtari (smokhtari@insa-toulouse.fr)

* **Last Modified by:** MIC **Date:** 11 December 2025

* **Last Reviewed by:** SM **Date:** 10 December 2025

"""


name = "nanoprisms"
title = "nanoprisms model with Fibonacci quadrature"
description = """
        nanoprisms with Fibonacci quadrature"""

category = "plugin"

import numpy as np
from numpy import inf
from prismformfactors import *

from sasmodels.quadratures.fibonacci import fibonacci_sphere

#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld", "1e-6/Ang^2", 126., [-inf, inf], "sld",
               "nanoprism scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 9.4, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["n_sides", "", 5, [3, 50], "",
               "n_sides"],
              ["R_ave", "Ang", 500, [0., inf], "volume",
               "Average radius"],
              ["L", "Ang", 5000, [0., inf], "volume",
               "length"],
               ]


def form_volume(nsides,Rave,L):
    nsides=int(nsides)
    edge = functions.edge_from_gyration_radius(nsides,Rave)
    return functions.volume_nanoprism(nsides,edge,L)

def Iqabc(qa,qb,qc,nsides,Rave,L): # proportionnal to the volume**2
    nsides=int(nsides)
    edge = functions.edge_from_gyration_radius(nsides,Rave)
    intensity=formfactor.I_nanoprism([qa,qb,qc],nsides,edge,L)
    return intensity


def Iq(q, sld, sld_solvent, nsides:int, Rave, L, npoints_fibonacci:int= 500):
    """
    Integration over the unit sphere with Fibonacci points.
    Each point has an equal weight = 1/npoints_fibonacci

    Parameters
    ----------
    q : float ou array
        Norm of the scattering vector
    sld, sld_solvent : 
        Contrast of scattering length density
    nsides, Rave, L : 
        Geometrical parameters of the prism
    npoints_fibonacci : int
        Number of Fibonacci points on the sphere
    Returns
    -------
    Iq : ndarray
        Scattering intensity averaged over all orientations
    time_fib : float
        Execution time of the function in seconds
    n_points_total : int
        Total number of points used in the quadrature
    """

    nsides = int(nsides)
    q = np.atleast_1d(q)
    q_unit,w = fibonacci_sphere(npoints_fibonacci)   # shape (npoints, 3)

    # Projections
    qa = q[:, np.newaxis] * q_unit[:, 0][np.newaxis, :]
    qb = q[:, np.newaxis] * q_unit[:, 1][np.newaxis, :]
    qc = q[:, np.newaxis] * q_unit[:, 2][np.newaxis, :]

    # Compute intensity
    intensity = Iqabc(qa, qb, qc, nsides, Rave, L)  # shape (nq, npoints)

    # Uniform average over the sphere
    integral = np.sum(w[np.newaxis, :] * intensity, axis=1)
    integral = np.mean(intensity, axis=1)


    return integral * (sld - sld_solvent)**2


Iq.vectorized = True




