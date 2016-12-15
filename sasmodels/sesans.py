"""
Conversion of scattering cross section from SANS (I(q), or rather, ds/dO) in absolute
units (cm-1)into SESANS correlation function G using a Hankel transformation, then converting
the SESANS correlation function into polarisation from the SESANS experiment

Everything is in units of metres except specified otherwise (NOT TRUE!!!)
Everything is in conventional units (nm for spin echo length)

Wim Bouwman (w.g.bouwman@tudelft.nl), June 2013
"""

from __future__ import division

import numpy as np  # type: ignore
from numpy import pi, exp  # type: ignore
from scipy.special import j0
#from mpmath import j0 as j0
        
class SesansTransform(object):
    #: Set of spin-echo lengths in the measured data
    SE = None  # type: np.ndarray
    #: Maximum acceptance of scattering vector in the spin echo encoding dimension (for ToF: Q of min(R) and max(lam))
    zaccept = None # type: float
    #: Maximum size sensitivity; larger radius requires more computation
    Rmax = None  # type: float
    #: q values to calculate when computing transform
    q = None  # type: np.ndarray

    # transform arrays
    _H = None  # type: np.ndarray
    _H0 = None # type: np.ndarray

    def set_transform(self, SE, zaccept, Rmax):
        if self.SE is None or len(SE) != len(self.SE) or np.any(SE != self.SE) or zaccept != self.zaccept or Rmax != self.Rmax:
            self.SE, self.zaccept, self.Rmax = SE, zaccept, Rmax
            self._set_q()
            self._set_hankel()

    def apply(self, Iq):
        G0 = np.dot(self._H0, Iq)
        G = np.dot(self._H.T, Iq)
        P = G - G0
        return P

    def _set_q(self):
        #q_min = dq = 0.1 * 2*pi / self.Rmax

        q_max = 2*pi / (self.SE[1]-self.SE[0])
        q_min = dq = 0.1 *2*pi / (np.size(self.SE) * self.SE[-1])

        #q_min = dq = q_max / 100000
        q=np.arange(q_min, q_max, q_min)
        self.q = q
        self.dq = dq

    def _set_hankel(self):
        #Rmax = #value in text box somewhere in FitPage?
        q = self.q
        dq = self.dq
        SElength = self.SE

        H0 = dq / (2 * pi) * q
        q=np.array(q,dtype='float32')
        SElength=np.array(SElength,dtype='float32')

        # Using numpy tile, dtype is conserved
        repq=np.tile(q,(SElength.size,1))
        repSE=np.tile(SElength,(q.size,1))
        H = dq / (2 * pi) * j0(repSE*repq.T)*repq.T

        # Using numpy meshgrid - meshgrid produces float64 from float32 inputs! Problem for 32-bit OS: Memerrors!
        #H0 = dq / (2 * pi) * q
        #repSE, repq = np.meshgrid(SElength, q)
        #repq=np.array(repq,dtype='float32')
        #repSE=np.array(repSE,dtype='float32')
        #H = dq / (2 * pi) * j0(repSE*repq)*repq

        self._H, self._H0 = H, H0

class SESANS1D(SesansTransform):
    def __init__(self, data, _H0, _H, q_calc):
        # x values of the data (Sasview requires these to be named "q")
        self.q = data.x
        self._H0 = _H0
        self._H = _H
        # Pysmear does some checks on the smearer object, these checks require the "data" object...
        self.data=data
        # q values of the SAS model
        self.q_calc = q_calc # These are the MODEL's q values used by the smearer (in this case: the Hankel transform)
    def apply(self, theory):
        return SesansTransform.apply(self,theory)