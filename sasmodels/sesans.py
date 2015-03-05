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

def make_q(q_zmax, Rmax):
    q_min = dq = 0.1 * 2*pi / Rmax
    #q_min = 0.00003
    return np.arange(q_min, q_zmax, dq)

# TODO: dead code; for now the call to the hankel transform happens in BumpsModel
class SesansCalculator:
    def __init__(self, kernel, q_zmax, Rmax, SElength, wavelength, thickness):
        self._set_kernel(kernel, q_zmax, Rmax)
        self.SElength = SElength
        self.wavelength = wavelength
        self.thickness = thickness

    def _set_kernel(self, kernel, q_zmax, Rmax):
        kernel_input = kernel.make_input([make_q(q_zmax, Rmax)])
        self.sans_calculator = kernel(kernel_input)

    def __call__(self, pars, pd_pars, cutoff=1e-5):
        Iq = self.sans_calculator(pars, pd_pars, cutoff)
        P = hankel(self.SElength, self.wavelength, self.thickness, self.q, Iq)
        self.Iq = Iq
        return P

def hankel(SElength, wavelength, thickness, q, Iq):
    """
    Compute the expected SESANS polarization for a given SANS pattern.

    Uses the hankel transform followed by the exponential.  The values
    for zz (or spin echo length, or delta), wavelength and sample thickness
    information should come from the dataset.  *q* should be chosen such
    that the oscillations in *I(q)* are well sampled (e.g., 5*2*pi/d_max).

    *SElength* [A] is the set of z points at which to compute the hankel transform

    *wavelength* [m]  is the wavelength of each individual point *zz*

    *thickness* [cm] is the sample thickness.

    *q* [A^{-1}] is the set of q points at which the model has been computed.
    These should be equally spaced.

    *I* [cm^{-1}] is the value of the SANS model at *q*
    """
    G = np.zeros(len(SElength), 'd')
    for i in range(len(SElength)):
        integr = besselj(0,q*SElength[i])*Iq*q
        G[i] = np.sum(integr)
    dq=(q[1]-q[0])*1e10   # [m^-1] step size in q, needed for integration
    G *= dq*1e10*2*pi # integr step, conver q into [m**-1] and 2 pi circle integr
    P = exp(thickness*wavelength**2/(4*pi**2)*(G-G[0]))

    return P
