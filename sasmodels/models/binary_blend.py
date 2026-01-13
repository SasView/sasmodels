r"""

Two-polymer RPA model with a flat background.

Definition
----------
This model calculates the scattering from a two polymer
blend using the Random Phase Approximation (RPA).

This is a revision of the *binary_blend* model posted to the
Marketplace in May 2020 following User feedback.

.. note:: The two polymers are assumed to be monodisperse and incompressible.

.. note:: There is an equivalent model in sasmodels called *rpa_binary_mixture_homopolymers*.
          While they do the same calculation (although *rpa_binary_mixture_homopolymers* is
          coded in C), they use different input parameters (see the validation section in
          the documentation of *rpa_binary_mixture_homopolymers* for details), so you can
          choose one or another depending on your preferences.

The scattered intensity $I(q)$ is calculated as [1,2]

.. math::

     \frac{(\rho_A - \rho_B)^2}{N_{Av} I(q)} = \text{scale} \cdot
     \left[\frac{1}{\phi_A n_A v_A P_A(q)} + \frac{1}{\phi_B n_B v_B P_B(q)} -
     \frac {2 \chi_{AB}}{v_0}\right] + \text{background}

where $\rho_i$ is the scattering length density, $\phi_i$ is the volume
fraction (such that $\phi_A$ + $\phi_B$ = 1), $n_i$ is the degree of
polymerization (number of repeat units), $v_i$ is the molar specific volume
of one *monomer*, $P_i(q)$ is Debye's Gaussian coil form factor
for polymer $i$, $N_{Av}$ is Avogadro's number, $\chi_{AB}$ is the Flory-Huggins
interaction parameter and $v_0$ is the reference volume,

with

.. math::

     v_i = m_i / \delta_i , \;
     v_0 = \sqrt{v_A \cdot v_B} , \;
     Z = (q \cdot Rg_i)^2 , \;
     P_i(q) = 2 \cdot [exp(-Z) + Z - 1] / Z^2 ,

where $m_i$ is the molecular weight of a *repeat unit* (not of the polymer),
$\delta_i$ is the mass density, and $Rg_i$ is the radius of gyration of polymer $i$.

.. note:: This model works best when as few parameters as possible are allowed
          to optimize. Indeed, most parameters should be known *a priori* anyhow!
          The calculation should also be exact, meaning the *scale* parameter
          should be left at 1.

.. note:: Tip: Try alternately optimizing $n_A$ & $Rg_A$ and $n_B$ & $Rg_B$ and
          only then optimizing $\chi_{AB}$.

Acknowledgments
----------------
The author would like to thank James Cresswell and Richard Thompson for
highlighting some issues with the model as it was originally coded.
The molar specific volumes were being computed from the polymer molecular
weight, not the weight of the repeat unit, meaning in most instances the
values were grossly over-estimated, whilst the reference volume was fixed at
a value which in most instances would have been too small. Both issues are
corrected in this version of the model.

References
----------

#.  Lapp, Picot & Benoit, *Macromolecules*, (1985), 18, 2437-2441 (Appendix) 
#.  Hammouda, *The SANS Toolbox*, Chapters 28, 29 & 34 and Section H

Authorship and Verification
---------------------------

* **Author:** Steve King **Date:** 07/05/2020
* **Last Modified by:** Steve King **Date:** 12/09/2024
* **Last Reviewed by:** Miguel A. Gonzalez **Date:** 16/11/2025

"""
import numpy as np
from numpy import inf

name = "binary_blend"
title = "Two-polymer RPA model"
description = """
    Evaluates a two-
    polymer RPA model
"""
category = "Polymers"
structure_factor = False
single = False

#   ["name", "units", default, [lower, upper], "type","description"],
#   lower limit for m_A/m_B set for CH repeat unit
#   lower limit for n_A/n_B set for PS (see Macromolecules, 37(1), 2004, 161-166);
#   could be reduced for more flexible polymers and vice versa
#   lower limit for rg_A/rg_B set for C-C bond!
parameters = [
    ["chi_AB", "cm^3/mol", 0.001, [0, 100], "", "Flory interaction parameter"],
    ["density_A", "g/cm^3", 1.05, [0.1, 3], "", "Mass density of A"],
    ["density_B", "g/cm^3", 0.90, [0.1, 3], "", "Mass density of B"],
    ["m_A", "g/mol", 112, [13, inf], "", "Mol weight of REPEAT A"],
    ["m_B", "g/mol", 104, [13, inf], "", "Mol weight of REPEAT B"],
    ["n_A", "", 465, [50, inf], "", "No. repeats of A"],
    ["n_B", "", 501, [50, inf], "", "No. repeats of B"],
    ["rg_A", "Ang", 59.3, [1.5, inf], "", "Radius of gyration of A"],
    ["rg_B", "Ang", 59.3, [1.5, inf], "", "Radius of gyration of B"],
    ["sld_A", "1e-6/Ang^2", 6.55, [-1, 7], "sld", "SLD of A"],
    ["sld_B", "1e-6/Ang^2", 1.44, [-1, 7], "sld", "SLD of B"],
    ["volfrac_A", "", 0.48, [0, 1], "", "Volume fraction of A"],
]

def Iq(q, chi_AB, density_A, density_B, m_A, m_B,
    n_A, n_B, rg_A, rg_B, sld_A, sld_B, volfrac_A):

    N_Av = 6.023E+23             # Avogadro number (mol^-1)
    v_A = m_A / density_A        # Molar specific volume of A (cm^3 mol^-1)
    v_B = m_B / density_B        # Molar specific volume of B (cm^3 mol^-1)
    v_0 = np.sqrt(v_A * v_B)     # Reference volume (cm^3 mol^-1)

    Z_A = (q * rg_A) * (q * rg_A)
    Z_B = (q * rg_B) * (q * rg_B)

    contrast = (sld_A - sld_B)   # Contrast (Ang^-2)

    Pq_A = 2.0 * (np.exp(-1.0 * Z_A) - 1.0 + Z_A) / (Z_A * Z_A)
    Pq_B = 2.0 * (np.exp(-1.0 * Z_B) - 1.0 + Z_B) / (Z_B * Z_B)

    U_A = volfrac_A * n_A * v_A * Pq_A
    U_B = (1.0 - volfrac_A) * n_B * v_B * Pq_B
    chiterm = (2.0 * chi_AB) / v_0
    prefactor = 1.0E+20 * (contrast * contrast) / N_Av
    Inverse_Iq = (1.0 / prefactor) * ((1.0 / U_A) + (1.0 / U_B) - chiterm)
    result = 1.0 / Inverse_Iq
    return result

Iq.vectorized = True  # Iq accepts an array of q values

tests = [
   [{'scale': 1.0, 'background' : 0.08, 'chi_AB': 0.0005, 'density_A': 1.05,
     'density_B': 0.9, 'm_A': 112, 'm_B': 104, 'n_A' : 400, 'n_B' : 351,
     'rg_A': 60, 'rg_B': 48, 'sld_A': 6.55, 'sld_B': 1.44, 'volfrac_A': 0.5},
     [0.002009078240301904, 0.24943453429220586], [49.594719935513766, 0.5713619727360871]],
   ]
