r"""
Definition
----------

This model calculates an empirical functional form for SAS data characterized
by a broad scattering peak. Many SAS spectra are characterized by a broad peak
even though they are from amorphous soft materials. For example, soft systems
that show a SAS peak include copolymers, polyelectrolytes, multiphase systems,
layered structures, etc.

The d-spacing corresponding to the broad peak is a characteristic distance
between the scattering inhomogeneities (such as in lamellar, cylindrical, or
spherical morphologies, or for bicontinuous structures).

The scattering intensity $I(q)$ is calculated as

.. math:: I(q) = \frac{A}{q^n} + \frac{C}{\left[1 + (|q - q_0|\xi)^m\right]^p} + B

Here the peak position is related to the d-spacing as $q_0 = 2\pi / d_0$.

$A$ is the Porod law scale factor, $n$ the Porod exponent, $C$ is the
peak scale factor, $m$ the peak width exponent, $\xi$ the correlation length,
$p$ the peak shape exponent, and $B$ the flat background. The exponent $p$
acts on the peak term only.

With $A = 0$ and $q_0 = 0$, the model reduces to the lorentz model for
$m = 2$, $p = 1$, and to the Debye-Anderson-Brumberger (dab) model for
$m = 2$, $p = 2$ (with $C = \xi^3$).

.. note::
   In sasmodels 1.1.0 and earlier (SasView 6.2.x and earlier), $p$ was
   applied to the sum of the power law and peak terms, contrary to
   the equation above. Fits with $p \ne 1$ from those versions will not
   be reproduced by later versions. Results for $p = 1$ are unchanged.

For 2D data the scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math:: q = \sqrt{q_x^2 + q_y^2}

References
----------

None.

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Adam H. Turner **Date:** November 8, 2026
* **Last Reviewed by:** Richard Heenan **Date:** March 21, 2016
"""

import numpy as np
from numpy import errstate, inf

name = "broad_peak"
title = "Broad peak on top of a power law decay"
description = """\
      I(q) = porod_scale/pow(q,porod_exp) + peak_scale/
      pow(1.0 + pow(fabs(q - peak_pos)*correlation_length, width_exp), shape_exp)
      + background

      List of default parameters:
      porod_scale = Porod term scaling
      porod_exp = Porod exponent
      peak_scale = Lorentzian term scaling
      correlation_length = Correlation length [A]
      peak_pos = peak location [1/A]
      width_exp = peak width exponent
      shape_exp = peak shape exponent
      background = Incoherent background"""
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["porod_scale",    "",  1.0e-05, [-inf, inf], "", "Power law scale factor"],
              ["porod_exp",      "",      3.0, [-inf, inf], "", "Exponent of power law"],
              ["peak_scale",  "",     10.0, [-inf, inf], "", "Scale factor for broad peak"],
              ["correlation_length", "Ang",  50.0, [-inf, inf], "", "screening length"],
              ["peak_pos",       "1/Ang", 0.1, [-inf, inf], "", "Peak position in q"],
              ["width_exp",    "",      2.0, [-inf, inf], "", "Exponent of peak width"],
              ["shape_exp",    "",      1.0, [-inf, inf], "", "Exponent of peak shape"],
             ]
# pylint: enable=bad-whitespace, line-too-long

def Iq(q,
       porod_scale=1.0e-5,
       porod_exp=3.0,
       peak_scale=10.0,
       correlation_length=50.0,
       peak_pos=0.1,
       width_exp=2.0,
       shape_exp=1.0):
    """
    :param q:                    Input q-value
    :param porod_scale:          Power law scale factor
    :param porod_exp:            Exponent of power law
    :param peak_scale:           Scale factor for broad peak
    :param correlation_length:   Correlation length
    :param peak_pos:             Peak position in q
    :param width_exp:            Exponent of peak width
    :param shape_exp:            Exponent of peak shape
    :return:                     Calculated intensity
    """
    z = abs(q - peak_pos) * correlation_length
    with errstate(divide='ignore'):
        inten = (porod_scale / q ** porod_exp
                 + peak_scale / (1 + z ** width_exp) ** shape_exp)
    return inten
Iq.vectorized = True  # Iq accepts an array of q values

def random():
    """Return a random parameter set for the model."""
    pars = dict(
        scale=1,
        porod_scale=10**np.random.uniform(-8, -5),
        porod_exp=np.random.uniform(1, 6),
        peak_scale=10**np.random.uniform(0.3, 6),
        correlation_length=10**np.random.uniform(0, 2),
        peak_pos=10**np.random.uniform(-3, -1),
        width_exp=np.random.uniform(1, 4),
        shape_exp=np.random.uniform(1, 2),
    )
    pars['correlation_length'] /= pars['peak_pos']
    pars['peak_scale'] *= pars['porod_scale'] / pars['peak_pos']**pars['porod_exp']
    #pars['porod_scale'] = 0.
    return pars

tests = [
    # on the peak, z = 0, so shape_exp must not change the result
    [{}, 0.1, 10.011],
    [{'shape_exp': 2.0}, 0.1, 10.011],
    # off the peak: 1e-5/0.12**3 + 10/(1 + 1)**2 + 0.001
    [{'shape_exp': 2.0}, 0.12, 2.506787037],
    # lorentz limit, same point and value as the lorentz model test
    [{'porod_scale': 0.0, 'peak_pos': 0.0, 'peak_scale': 1.0,
      'correlation_length': 250.0, 'width_exp': 2.0, 'shape_exp': 1.0},
     0.01, 0.138931],
    # dab limit with L = 50: L**3/(1 + (q*L)**2)**2 + 0.001
    [{'porod_scale': 0.0, 'peak_pos': 0.0, 'peak_scale': 125000.0,
      'correlation_length': 50.0, 'width_exp': 2.0, 'shape_exp': 2.0},
     0.01, 80000.001],
]
