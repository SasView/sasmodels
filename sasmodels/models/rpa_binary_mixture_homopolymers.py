r"""

RPA model for a binary mixture of homopolymers

Definition
----------

Calculates the macroscopic scattering intensity for a
two polymer blend using the Random Phase Approximation (RPA).

.. note:: There is an equivalent model in sasmodels called *binary_blend*.
          While they do the same calculation (although *binary_blend* is
          coded purely in python), they use different input parameters (see
          the validation section below for details), so you can choose
          one or another depending on your preferences.

The model is based on the papers by Akcasu *et al.* [1] and by
Hammouda [2], assuming the polymer follows Gaussian statistics such
that $R_g^2 = n b^2/6$, where $b$ is the statistical segment length and $n$ is
the number of statistical segment lengths. A nice tutorial on how these are
constructed and implemented can be found in chapters 28, 31 and 34, and Part H,
of Hammouda's 'SANS Toolbox' [3].

The scattered intensity $I(q)$ is calculated as[2,3]

.. math::

     \frac{\left(l_A/v_A - l_B/v_B\right)^2 N_{Av}}{I(q)} = \text{scale} \cdot
     \left[\frac{1}{\phi_A n_A v_A P_A(q)} + \frac{1}{\phi_B n_B v_B P_B(q)} -
     \frac{2 \chi_{AB}}{v_0}\right] + \text{background}

where $l_i$, $v_i$, and $b_i$ are the scattering length,
molar specific volume and segment length of a single *monomer*,
$\phi_i$ is the volume fraction (noting that $\phi_A$ + $\phi_B$ = 1),
$n_i$ is the degree of polymerization (number of repeat units),
$N_{Av}$ is Avogadro's number, $\chi_{AB}$ is the Flory-Huggins
interaction parameter, $v_0 = \sqrt{v_A \cdot v_B}$ is the
reference volume, and $P_i(q)$ is Debye's Gaussian coil form
factor for an isolated chain of polymer $i$:

.. math::

     P_i(q) = \frac{2 \cdot \left[exp(-q^2 \cdot Rg_i^2) + q^2 \cdot Rg_i^2 - 1\right]}
     {\left(q \cdot Rg_i\right)^4}

.. note::
    * In general, the degrees of polymerization, the volume
      fractions, the molar volumes and the scattering lengths for each
      component are obtained from other methods and held fixed, while the *scale*
      parameter should be held equal to unity.
    * The variables to fit are normally the segment lengths $b_a$ and $b_b$,
      and the Flory-Huggins interaction parameter $\chi_{AB}$.

Validation
----------

The default parameters correspond to the example case of a binary blend
of deuterated polysterene (PSD) and poly (vinyl methyl ether) (PVME) at 60 C,
reported in section 8.2 of [2]. Thus, the generated scattering curve
can be compared with the one shown in Fig. 11 of [2].

The same result is also obtained with the *binary_blend* model, which
uses alternative input entries:

  * Mass density and molecular weight instead of molar volume (use 1.12 $g/cm^3$
    and 112 g/mol for PSD and 1.05 $g/cm^3$ and 58 g/mol for PVME).

  * Radius of gyration instead of the segment length of the monomer (use 136.3 |Ang| for
    PSD and 128.2 |Ang| for PVME).

  * Scattering length densities instead of scattering lengths (use 6.42 and 0.594
    (in units of $10^{-6}/A^2$) for PSD and PVME, respectively).

.. warning:: The original *rpa* model applied to case 0 did not produce the same
             result. However, at present is not completely clear if this is due to
             a problem with the code or to the problems with the interface to select
             and pass the correct parameters.

References
----------

#. A Z Akcasu, R Klein and B Hammouda, *Macromolecules*, 26 (1993) 4136
#. B. Hammouda, *Advances in Polymer Science* 106 (1993) 87
#. B. Hammouda, *SANS Toolbox*
   https://www.nist.gov/system/files/documents/2023/04/14/the_sans_toolbox.pdf
   (accessed 15 November 2025).

Authorship and Verification
----------------------------

* **Author:** Boualem Hammouda - NIST IGOR/DANSE **Date:** pre 2010
* **Converted to sasmodels by:** Paul Kienzle **Date:** July 18, 2016
* **Single case rewritten by:** Miguel A. Gonzalez **Date:** November 16, 2025
"""

from numpy import inf

name = "rpa_binary_mixture_homopolymers"
title = "Random Phase Approximation for a binary mixture of homopolymers"
description = """
Intensity of a binary mixture of homopolymers calculated within the RPA
"""
category = "Polymers"

#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["phiA", "", 0.484, [0, 1], "", "volume fraction of polymer A"],
    ["nA", "", 1741.0, [1, inf], "", "Degree of polymerization of polymer A"],
    ["vA", "cm^3/mol", 100.0, [0, inf], "", "molar volume of polymer A"],
    ["lA", "fm", 106.5, [-inf, inf], "", "scattering length of polymer A"],
    ["bA", "Ang", 8.0, [0, inf], "", "segment length of polymer A"],
    ["nB", "", 2741.0, [1, inf], "", "Degree of polymerization of polymer B"],
    ["vB", "cm^3/mol", 55.4, [0, inf], "", "molar volume of polymer B"],
    ["lB", "fm", 5.5, [-inf, inf], "", "scattering length of polymer B"],
    ["bB", "Ang", 6.0, [0, inf], "", "segment length of polymer B"],
    ["chiAB", "cm^3/mol", -0.021, [-inf, inf], "", "A:B interaction parameter"],
]

source = ["rpa_binary_mixture_homopolymers.c"]
single = False

# TODO: no random parameters generated for RPA
