r"""
This model describes a pseudo-Voigt shaped peak on a flat background.

Definition
----------

This pseudo-Voigt peak function is a weighted linear summation of
Lorentzian (L) and Gaussian (G) peak shapes. 
It is a popular function for modelling peak shape.
It can be tailored to any experimental peak shape. 

The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = scale \cdot \left[ w_f \cdot I(q)_L + (1 - w_f) \cdot I(q)_G \right] + background

where $w_f$ is a weighting factor and

.. math::

    I(q)_L = \frac{1}{1 + \left( \frac{q - q_0}{HWHM} \right)^2}

    I(q)_G = \exp\left[ -\frac{1}{2} (q - q_0)^2 / \sigma^2 \right]

The peak is taken to be centered at $q_0$ with a HWHM (half-width
half-maximum) of $1.17741\,\sigma$, where $\sigma$ is the standard deviation
of the Gaussian. In other words, the widths of the Lorentzian and the
Gaussian have been coupled for convenience of parameterisation:

.. math::

    \sigma = HWHM / \sqrt{2 \ln 2} = HWHM / 1.17741

When $w_f = 1$ a Lorentzian peak is returned, and when $w_f = 0$ a
Gaussian peak is returned.

The pseudo-Voigt is an approximation to the true Voigt profile, which is the
convolution of a Lorentzian with a Gaussian. Because the two widths are coupled
here, the lineshape is controlled by the single weighting factor $w_f$, so the model
has the same number of free parameters as a true Voigt would: peak position, width
and $w_f$, against peak position, $\sigma$ (Gaussian) and $\gamma$ (Lorentzian).
Varying $w_f$ at fixed HWHM changes the weight in the tails rather than the width of
the peak. The advantage over the true Voigt is that the expression is analytic and
inexpensive to evaluate.

For 2D data the scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


Note on instrumental resolution
-------------------------------

As for the other peak models in sasmodels, $HWHM$ is the width of the peak produced
by the sample alone. The width actually observed is larger, because the measured
intensity is the model convolved with the instrumental resolution function. The
fitted $HWHM$ therefore only corresponds to the measured peak width when the
resolution is negligible; otherwise resolution smearing should be applied during the
fit, so that the fitted value remains the physical width.


Validation
----------

The pseudo-Voigt peak reduces exactly to a pure Lorentzian for $w_f = 1$
and to a pure Gaussian for $w_f = 0$; both limits were checked against their
analytic values (see tests section at the end).
The full pseudo-Voigt shape has also been compared, for identical
parameters, against a slightly different SasView implementation (https://marketplace.sasview.org/models/127/) 
of the same function and gives the same result. 


References
----------

1. Aaron L. Stancik, Eric B. Brauns
   A simple asymmetric lineshape for fitting infrared absorption spectra
   Vibrational Spectroscopy 47 (2008) 66-69


Authorship and Verification
----------------------------

* **Author:**  Steve King **Date:** 24 June 2020

* **Last Modified by:** Anirban Mandal (mandalanirban2023@gmail.com) **Date:** 06 July 2026

* **Last Reviewed by:** Marianne Imperor-Clerc (marianne.imperor@cnrs.fr) **Date:** 07 September 2026

"""

import numpy as np
from numpy import inf

name = "peak_pseudo_voigt"
title = "Single pseudo-Voigt peak"
description = """\
          I(q) = scale*peak + background
"""

category = "shape-independent"

parameters = [["w_f", "", 0.8, [0, 1], "", "lorentzian/gaussian weighting factor"],
              ["peak_pos", "1/Ang", 0.05, [0, inf], "", "Position of the peak"],
              ["peak_hwhm", "1/Ang", 0.01, [0, 1], "", "HWHM of the peak"]]


# Kept as a separate function so that models of ordered phases (lamellar,
# 2D hexagonal, ...) can call it once per reflection. Moving it to
# sasmodels.special as sas_pseudovoigt would be the cleaner home for it.
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


def Iq(q, w_f, peak_pos, peak_hwhm):
    """
    w_f: weighting coefficient in the pseudo-Voigt peak function;
    w_f = 1 for a Lorentzian and w_f = 0 for a Gaussian peak.
    peak_pos: position of the peak
    peak_hwhm: HWHM of the peak
    """
    return Ipeak(q, w_f, peak_pos, peak_hwhm)

Iq.vectorized = True  # Iq accepts an array of q values

tests = [
    # pure Lorentzian (w_f = 1): peak centre, half-width, and 2 x HWHM
    [{"scale": 1.0, "background": 0.0, "w_f": 1.0,
      "peak_pos": 0.05, "peak_hwhm": 0.01}, 0.05, 1.0],
    [{"scale": 1.0, "background": 0.0, "w_f": 1.0,
      "peak_pos": 0.05, "peak_hwhm": 0.01}, 0.06, 0.5],
    [{"scale": 1.0, "background": 0.0, "w_f": 1.0,
      "peak_pos": 0.05, "peak_hwhm": 0.01}, 0.07, 0.2],
    # pure Gaussian (w_f = 0): half-width is 0.5 by definition, 2 x HWHM = 1/16
    [{"scale": 1.0, "background": 0.0, "w_f": 0.0,
      "peak_pos": 0.05, "peak_hwhm": 0.01}, 0.06, 0.5],
    [{"scale": 1.0, "background": 0.0, "w_f": 0.0,
      "peak_pos": 0.05, "peak_hwhm": 0.01}, 0.07, 0.0625],
    # mixed pseudo-Voigt (w_f = 0.8) away from the centre
    [{"scale": 1.0, "background": 0.0, "w_f": 0.8,
      "peak_pos": 0.05, "peak_hwhm": 0.01}, 0.07, 0.1725],
    # zero width: unity at the centre, zero elsewhere
    [{"scale": 1.0, "background": 0.0, "w_f": 0.8,
      "peak_pos": 0.05, "peak_hwhm": 0.0}, 0.05, 1.0],
    [{"scale": 1.0, "background": 0.0, "w_f": 0.8,
      "peak_pos": 0.05, "peak_hwhm": 0.0}, 0.06, 0.0],
]
