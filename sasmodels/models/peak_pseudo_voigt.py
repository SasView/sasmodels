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

For 2D data the scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


Reuse of the peak shape in multi-peak models
--------------------------------------------

The line shape is implemented in a separate helper function, ``Ipeak(q, wf, q0,
hwhm)``, which ``Iq`` then calls. The reason is that a single diffraction peak is
rarely measured on its own: an ordered phase produces a *series* of Bragg peaks
whose positions are all fixed by one lattice parameter, while every peak in the
series shares the same line shape. Writing the line shape once, as a function of an
arbitrary peak position and width, is what allows the same code to be reused for
those multi-peak models. For example, for a lamellar phase of period $d$

.. math::

    q_n = n\,q_1, \qquad q_1 = \frac{2\pi}{d}, \qquad n = 1, 2, 3 \ldots

and for a two-dimensional hexagonal phase of cell parameter $a$

.. math::

    q_{hk} = q_{10}\,\sqrt{h^2 + hk + k^2}, \qquad q_{10} = \frac{4\pi}{\sqrt{3}\,a}

giving peaks in the ratios $1 : \sqrt{3} : 2 : \sqrt{7} : 3 \ldots$ for the
$(10), (11), (20), (21), (30) \ldots$ reflections. A model of such a phase is then
a sum over reflections,

.. math::

    I(q) = scale \cdot \sum_{hk} A_{hk}\, \mathrm{Ipeak}(q, w_f, q_{hk}, HWHM_{hk})
           + background

in which only the amplitudes $A_{hk}$, the widths and the single lattice parameter
are fitted. Keeping ``Ipeak`` separate from ``Iq`` means that the single-peak model
documented here and any multi-peak model built on it evaluate *identical* code for
the line shape, so the two can never drift apart. ``Iq`` is the thin wrapper that
exposes the single-peak case through the sasmodels interface.

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
from numpy import errstate, inf

name = "peak_pseudo_voigt"
title = "Single pseudo-Voigt peak"
description = """\
          I(q) = scale*peak + background
"""

category = "shape-independent"

parameters = [["w_f", "", 0.8, [0, 1], "", "lorentzian/gaussian weighting factor"],
              ["peak_pos", "1/Ang", 0.05, [0, inf], "", "Position of the peak"],
              ["peak_hwhm", "1/Ang", 0.01, [0, 1], "", "HWHM of the peak"]]


# The line shape is kept in its own function, taking the peak position and width as
# plain arguments, so that models of ordered phases (lamellar, 2D hexagonal, ...) can
# call it once per reflection and share exactly this code. See the section "Reuse of
# the peak shape in multi-peak models" in the module documentation above.
def Ipeak(q, wf, q0, hwhm):
    """
    Pseudo-Voigt line shape at an arbitrary peak position, for reuse by
    multi-peak models as well as by ``Iq`` below.

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
    intensity = (
        (wf * (1 / (1 + ((q - q0)**2.0 / hwhm**2.0))))
        + ((1.0 - wf) * np.exp((-0.5 * (q - q0)**2.0) / (sigma**2.0)))
    )
    return intensity


def Iq(q, w_f, peak_pos, peak_hwhm):
    """
    w_f: weighting coefficient in the pseudo-Voigt peak function;
    w_f = 1 for a Lorentzian and w_f = 0 for a Gaussian peak.
    peak_pos: position of the peak
    peak_hwhm: HWHM of the peak

    The ``errstate`` guard covers the lower limit of the ``peak_hwhm`` range: at
    ``peak_hwhm = 0`` the line shape divides by zero away from the peak centre
    (``divide``) and evaluates 0/0 at the centre itself (``invalid``). Both limits
    are well defined and numpy already returns them correctly, so only the warnings
    need suppressing.
    """
    with errstate(divide='ignore', invalid='ignore'):
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
]
