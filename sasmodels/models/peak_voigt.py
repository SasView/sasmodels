r"""
This model describes a pseudo-Voigt shaped peak on a flat background.

Definition
----------

This pseudo-Voigt peak function is a weighted linear summation of
Lorentzian (L) and Gaussian (G) peak shapes. 
It is a popular function for modelling peak shape.
It can be tailored to any specific peak shape and it can also produce a peak shape with asymmetry. 

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

For 2D data the scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


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

1. L A Feigin, D I Svergun, G W Taylor
   Structure Analysis by Small-Angle X-ray and Neutron Scattering
   Springer (1987)

2. Aaron L. Stancik, Eric B. Brauns
   A simple asymmetric lineshape for fitting infrared absorption spectra
   Vibrational Spectroscopy 47 (2008) 66-69


Authorship and Verification
----------------------------

* **Author:**  Steve King **Date:** 24 June 2020

* **Authors:** Marianne Imperor-Clerc (marianne.imperor@cnrs.fr)
               Anirban Mandal (mandalanirban2023@gmail.com)

* **Last Modified by:** Anirban Mandal **Date:** 06 July 2026

* **Last Reviewed by:** Steve King **Date:**

"""

import numpy as np
from numpy import inf, errstate

name = "peak_voigt"
title = "Single pseudo-Voigt peak"
description = """\
          I(q) = scale*peak + background
"""

category = "shape-independent"

parameters = [["w_f", "", 0.8, [0, 1], "", "lorentzian/gaussian weighting factor"],
              ["peak_pos", "1/Ang", 0.05, [0, inf], "", "Position of the peak"],
              ["peak_hwhm", "1/Ang", 0.01, [0, 1], "", "HWHM of the peak"]]


def Ipeak(q, wf, q0, hwhm):
    """
    When $w_f$ = 1 a Lorentzian peak is returned, and when $w_f$ = 0 a
    Gaussian peak is returned.

    The peak is taken to be centered at $q_0$ with a HWHM (half-width
    half-maximum) for the Lorentzian and sigma = HWHM / 1.17741 for the
    Gaussian, where sigma is the standard deviation of the Gaussian. In
    other words, the widths of the Lorentzian and the Gaussian have been
    coupled for convenience of parameterisation.
    """
    cste = np.sqrt(2 * np.log(2))
    # cste = 1.17741
    sigma = hwhm / cste
    intensity = (wf * (1 / (1 + ((q - q0)**2.0 / hwhm**2.0)))) + \
                ((1.0 - wf) * np.exp((-0.5 * (q - q0)**2.0) / (sigma**2.0)))
    return intensity


def Iq(q, w_f, peak_pos, peak_hwhm):
    """
    w_f: weighting coefficient in the pseudo-Voigt peak function;
    w_f = 1 for a Lorentzian and w_f = 0 for a Gaussian peak.
    peak_pos: position of the peak
    peak_hwhm: HWHM of the peak
    """

    with errstate(divide='ignore'):
        L = Ipeak(q, w_f, peak_pos, peak_hwhm)

        return L

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
]
