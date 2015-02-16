"""
Conversion of scattering cross section from SANS in absolute
units into SESANS using a Hankel transformation

Everything is in units of metres except specified otherwise

Wim Bouwman (w.g.bouwman@tudelft.nl), June 2013
"""

from __future__ import division

import numpy as np
from numpy import pi, exp
from scipy.special import jv as besselj

def hankel(zz, wavelength, th, q, I):
    """
    Compute the expected SESANS polarization for a given SANS pattern.

    Uses the hankel transform followed by the exponential.  The values
    for zz (or spin echo length, or delta), wavelength and sample thickness
    information should come from the dataset.  *q* should be chosen such
    that the oscillations in *I(q)* are well sampled (e.g., 5*2*pi/d_max).

    *zz* [nm] is the set of z points at which to compute the hankel transform

    *wavelength* [m]  is the wavelength of each individual point *zz*

    *th* [m] is the sample thickness.

    *q* [nm^{-1}] is the set of q points at which the model has been computed.
    These should be equally spaced.

    *qwidth* is the width of the integration

    *I* [m^{-1}] is the value of the SANS model at *q*
    """
    dq=(q[1]-q[0])*1e9   # [m^-1] step size in q, needed for integration
    G = np.zeros(len(zz), 'd')
    for i in range(len(zz)):
        integr = besselj(0,q*zz[i])*I*q
        G[i] = np.sum(integr)
    G *= dq*1e9*2*pi # integr step, conver q into [m**-1] and 2 pi circle integr
    PP = exp(th*wavelength**2/(4*pi**2)*(G-G[0]))

    return PP
