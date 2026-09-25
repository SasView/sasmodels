r"""
This model adds a $q^{-4}$ Porod regime and hexagonal Bragg peaks with the first
ten peaks present.

It describes objects assembled together with internal 2D hexagonal order, such as
rods packed side by side: the lattice is periodic in the plane perpendicular to the
rod axis and has no periodicity along it, so the reflections are indexed by two
integers $(h,k)$ only. The domains are assumed to be randomly oriented in the beam,
so that the 2D reciprocal lattice gives rings which azimuthally average to the peak
positions below. The Porod term accounts for the contribution at low $q$ due to the
surface of the assemblies.

Definition
----------

Bragg peaks are modelled using the pseudo-Voigt peak function, a weighted linear
summation of a Lorentzian (L) and a Gaussian (G) sharing the same half-width at
half-maximum:

.. math::

    V(q) = w_f \cdot L(q) + (1 - w_f) \cdot G(q)

with

.. math::

    L(q) = \frac{1}{1 + \left( \frac{q - q_0}{HWHM} \right)^2}

    G(q) = \exp\left[ -\frac{1}{2} (q - q_0)^2 / \sigma^2 \right]

    \sigma = HWHM / \sqrt{2 \ln 2}

The peak positions of a 2D hexagonal lattice of cell parameter $a$ all follow from
that single parameter,

.. math::

    q_{hk} = q_{10}\,\sqrt{h^2 + hk + k^2}, \qquad q_{10} = \frac{4\pi}{\sqrt{3}\,a}

so the first ten allowed reflections, for which $h^2+hk+k^2 = 1, 3, 4, 7, 9, 12, 13,
16, 19, 21$, are the (10), (11), (20), (21), (30), (22), (31), (40), (32) and (41)
peaks, at positions in the ratios $1 : \sqrt{3} : 2 : \sqrt{7} : 3 \ldots$ The
integers absent from that list (2, 5, 6, 8, ...) cannot be written as
$h^2+hk+k^2$, so no peak appears at those positions; that pattern of present and
absent reflections is the signature of the hexagonal symmetry.

Here $a$ (``a_cell``) is the centre-to-centre distance between neighbouring
objects, the rows being separated by $d_{10} = a\sqrt{3}/2$.

The scattered intensity is the sum of the Porod term and the ten peaks,

.. math::

    I(q) = scale \cdot \left[ \left(\frac{scale_{Porod}}{q}\right)^4
           + \sum_{hk} A_{hk} V(q, q_{hk}, HWHM_{hk}) \right] + background

where the amplitudes $A_{hk}$ and the widths $HWHM_{hk}$ are free for each
reflection, while the positions are fixed by $a$ alone.

Note on instrumental resolution
-------------------------------

As for the other peak models in sasmodels, the fitted $HWHM$ values are the widths
produced by the sample alone. The widths actually observed are larger, because the
measured intensity is the model convolved with the instrumental resolution function.

References
----------

1. G Porod. *Kolloid Zeit* 124 (1951) 83

2. L A Feigin, D I Svergun, G W Taylor
   Structure Analysis by Small-Angle X-ray and Neutron Scattering
   Springer (1987)

3. Aaron L. Stancik, Eric B. Brauns
   A simple asymmetric lineshape for fitting infrared absorption spectra
   Vibrational Spectroscopy 47 (2008) 66-69


Authorship and Verification
----------------------------

* **Authors:** Jules Marcone (julesmarcone@gmail.com) **Date:** 30 May 2023
               Marianne Imperor-Clerc (marianne.imperor@cnrs.fr)
               Anirban Mandal (mandalanirban2023@gmail.com)
* **Author:**  Steve King **Date:** 24 June 2020

* **Last Modified by:** MIC **Date:** 16 June 2026

* **Last Reviewed by:** **Date:**

"""

import numpy as np
from numpy import errstate, inf

name = "hexagonal_2d_10peaks_porod"
title = "2D hexagonal lattice: ten pseudo-Voigt Bragg peaks plus a Porod term"
description = """\
          I(q) = scale((scale_Porod/q)^4 + sum_of_hexagonal_peaks) + background
"""

category = "shape-independent"

parameters = [["scale_Porod", "", 0.05, [0, inf], "", "Scale factor for Porod"],
              ["a_cell", "Ang", 900, [20, 10000], "", "hexagonal cell parameter"],
              ["w_f", "", 0.8, [0, 1], "", "lorentzian/gaussian weighting factor"],
              ["hwhm_q10", "1/Ang", 1, [0, 1], "", "HWHM of q10 peak"],
              ["hwhm_q11", "1/Ang", 1, [0, 1], "", "HWHM of q11 peak"],
              ["hwhm_q20", "1/Ang", 1, [0, 1], "", "HWHM of q20 peak"],
              ["hwhm_q21", "1/Ang", 1, [0, 1], "", "HWHM of q21 peak"],
              ["hwhm_q30", "1/Ang", 1, [0, 1], "", "HWHM of q30 peak"],
              ["hwhm_q22", "1/Ang", 1, [0, 1], "", "HWHM of q22 peak"],
              ["hwhm_q31", "1/Ang", 1, [0, 1], "", "HWHM of q31 peak"],
              ["hwhm_q40", "1/Ang", 1, [0, 1], "", "HWHM of q40 peak"],
              ["hwhm_q32", "1/Ang", 1, [0, 1], "", "HWHM of q32 peak"],
              ["hwhm_q41", "1/Ang", 1, [0, 1], "", "HWHM of q41 peak"],
              ["scale_q10", "", 1, [0, inf], "", "Scale factor for q10 peak"],
              ["scale_q11", "", 1, [0, inf], "", "Scale factor for q11 peak"],
              ["scale_q20", "", 1, [0, inf], "", "Scale factor for q20 peak"],
              ["scale_q21", "", 0, [0, inf], "", "Scale factor for q21 peak"],
              ["scale_q30", "", 0, [0, inf], "", "Scale factor for q30 peak"],
              ["scale_q22", "", 0, [0, inf], "", "Scale factor for q22 peak"],
              ["scale_q31", "", 0, [0, inf], "", "Scale factor for q31 peak"],
              ["scale_q40", "", 0, [0, inf], "", "Scale factor for q40 peak"],
              ["scale_q32", "", 0, [0, inf], "", "Scale factor for q32 peak"],
              ["scale_q41", "", 0, [0, inf], "", "Scale factor for q41 peak"]]


# Pseudo-Voigt lineshape at an arbitrary peak position, called once per
# reflection below. The same function is used by the single-peak model
# peak_pseudo_voigt; sasmodels.special would be the cleaner home for it.
def Ipeak(q, wf, q0, hwhm):
    # wf = 1 gives a Lorentzian, wf = 0 a Gaussian; the two widths are coupled
    # through sigma = hwhm / sqrt(2 ln 2) so that both have the same HWHM.
    sigma = hwhm / np.sqrt(2 * np.log(2))
    # Protect against zero width peaks (hwhm == 0, sigma == 0), which are zero
    # everywhere except at the centre.
    lorentzian = (1.0 / (1.0 + (q - q0)**2 / hwhm**2) if hwhm > 0
                  else 1.0 * (q == q0))
    gaussian = (np.exp(-0.5 * (q - q0)**2 / sigma**2) if sigma > 0
                else 1.0 * (q == q0))
    return wf * lorentzian + (1.0 - wf) * gaussian


def Iq(q, scale_Porod, a_cell, w_f,
       hwhm_q10, hwhm_q11, hwhm_q20, hwhm_q21, hwhm_q30,
       hwhm_q22, hwhm_q31, hwhm_q40, hwhm_q32, hwhm_q41,
       scale_q10, scale_q11, scale_q20, scale_q21, scale_q30,
       scale_q22, scale_q31, scale_q40, scale_q32, scale_q41):
    """
    scale_Porod: scale coefficient for the Porod term, (scale_Porod/q)^4
    a_cell: 2D hexagonal cell parameter; it fixes all ten peak positions
    w_f: weighting coefficient in the pseudo-Voigt peak function;
    w_f = 1 for a Lorentzian and w_f = 0 for a Gaussian peak.
    hwhm_qhk: HWHM of the (hk) hexagonal peak
    scale_qhk: scale factor for the (hk) hexagonal peak
    """
    # All ten positions follow from the single cell parameter a_cell.
    q10 = 4 * np.pi / (np.sqrt(3) * a_cell)
    q11, q20, q21, q30 = np.sqrt(3) * q10, 2 * q10, np.sqrt(7) * q10, 3 * q10
    q22, q31, q40 = np.sqrt(12) * q10, np.sqrt(13) * q10, 4 * q10
    q32, q41 = np.sqrt(19) * q10, np.sqrt(21) * q10

    # The Porod term diverges at q = 0, which is outside any measured range.
    with errstate(divide='ignore'):
        porod = (scale_Porod / q)**4

    return (porod
            + scale_q10 * Ipeak(q, w_f, q10, hwhm_q10)
            + scale_q11 * Ipeak(q, w_f, q11, hwhm_q11)
            + scale_q20 * Ipeak(q, w_f, q20, hwhm_q20)
            + scale_q21 * Ipeak(q, w_f, q21, hwhm_q21)
            + scale_q30 * Ipeak(q, w_f, q30, hwhm_q30)
            + scale_q22 * Ipeak(q, w_f, q22, hwhm_q22)
            + scale_q31 * Ipeak(q, w_f, q31, hwhm_q31)
            + scale_q40 * Ipeak(q, w_f, q40, hwhm_q40)
            + scale_q32 * Ipeak(q, w_f, q32, hwhm_q32)
            + scale_q41 * Ipeak(q, w_f, q41, hwhm_q41))

Iq.vectorized = True  # Iq accepts an array of q values

tests = [
    # Porod term alone: scale*(scale_Porod/q)^4 + background
    # = 1e-5*(1/0.04)^4 + 0.01 = 3.90625 + 0.01
    [{"scale": 0.00001,
      "background": 0.01,
      "scale_Porod": 1.,
      "scale_q10": 0.,
      "scale_q11": 0.,
      "scale_q20": 0.,
      "scale_q21": 0.,
      "scale_q30": 0.,
      "scale_q22": 0.,
      "scale_q31": 0.,
      "scale_q40": 0.,
      "scale_q32": 0.,
      "scale_q41": 0.,
      },
     0.04, 3.916250],
    # (10) peak alone, at its centre: q10 = 4 pi / (sqrt(3) a_cell) with
    # a_cell = 900 Ang, peak height = scale*scale_q10 = 1.0, plus background
    [{"scale": 1.0,
      "background": 0.0,
      "scale_Porod": 0.,
      "a_cell": 900.,
      "w_f": 0.8,
      "hwhm_q10": 0.001,
      "scale_q10": 1.,
      "scale_q11": 0.,
      "scale_q20": 0.,
      "scale_q21": 0.,
      "scale_q30": 0.,
      "scale_q22": 0.,
      "scale_q31": 0.,
      "scale_q40": 0.,
      "scale_q32": 0.,
      "scale_q41": 0.,
      },
     0.0080613, 1.0],
]
