r"""
This model describes a Gaussian shaped peak on a flat background

.. image:: img/image198.PNG

with the peak having height of *I0* centered at *q0* and having a standard deviation of *B*. The FWHM (full-width
half-maximum) is 2.354 B.

For 2D data: The 2D scattering intensity is calculated in the same way as 1D, where the *q* vector is defined as

.. image:: img/image040.gif

.. image:: img/image199.jpg

*Figure. 1D plot using the default values (w/500 data points).*

REFERENCE
---------

None.
"""

from numpy import inf

name = "gaussian_peak"
title = "Gaussian shaped peak"
description = """
    Model describes a Gaussian shaped peak including a flat background
    Provide F(q) = scale*exp( -1/2 *[(q-q0)/B]^2 )+ background
"""
category = "shape-independent"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["q0", "1/Ang", 0.05, [-inf, inf], "", "Peak position"],
              ["sigma", "1/Ang", 0.005, [-inf, inf], "",
               "Peak width (standard deviation)"],
             ]

# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return 1.0;
    """

Iq = """
    return exp(-0.5*pow((q - q0)/sigma,2.0));
    """


Iqxy = """
    // never called since no orientation or magnetic parameters.
    //return -1.0;
    return Iq(sqrt(qx*qx + qy*qy), q0, sigma);
    """


# VR defaults to 1.0

demo = dict(scale=1, background=0, q0=0.05, sigma=0.005)
oldname = "PeakGaussModel"
oldpars = dict(sigma='B')
