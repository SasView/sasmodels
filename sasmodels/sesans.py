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

def make_q(q_max, Rmax):
    r"""
    Return a $q$ vector suitable for SESANS covering from $2\pi/ (10 R_{\max})$
    to $q_max$.
    """
    q_min = dq = 0.1 * 2*pi / Rmax
    return np.arange(q_min, q_max, dq)

def hankel(SElength, wavelength, thickness, q, Iq):
    r"""
    Compute the expected SESANS polarization for a given SANS pattern.

    Uses the hankel transform followed by the exponential.  The values for *zz*
    (or spin echo length, or delta), wavelength and sample thickness should
    come from the dataset.  $q$ should be chosen such that the oscillations
    in $I(q)$ are well sampled (e.g., $5 \cdot 2 \pi/d_{\max}$).

    *SElength* [A] is the set of $z$ points at which to compute the
    Hankel transform

    *wavelength* [m]  is the wavelength of each individual point *zz*

    *thickness* [cm] is the sample thickness.

    *q* [A$^{-1}$] is the set of $q$ points at which the model has been
    computed. These should be equally spaced.

    *I* [cm$^{-1}$] is the value of the SANS model at *q*
    """
    G = np.zeros(len(SElength), 'd')
    for i in range(len(SElength)):
        integr = besselj(0, q*SElength[i])*Iq*q
        G[i] = np.sum(integr)

    # [m^-1] step size in q, needed for integration
    dq=(q[1]-q[0])*1e10

    # integration step, convert q into [m**-1] and 2 pi circle integration
    G *= dq*1e10*2*pi

    P = exp(thickness*wavelength**2/(4*pi**2)*(G-G[0]))

    return P
