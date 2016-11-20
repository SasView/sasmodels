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
#from scipy.special import jv as besselj
from scipy.special import j0
#from mpmath import j0 as j0
#from mpmath import besselj
#from mpmath import mpf
from src.sas.sascalc.data_util.nxsunit import Converter
#from sas.sasgui.perspectives.fitting.fitpage import FitPage
#import direct_model.DataMixin as model
        
def make_q(q_max, Rmax):
    r"""
    Return a $q$ vector suitable for SESANS covering from $2\pi/ (10 R_{\max})$
    to $q_max$. This is the integration range of the Hankel transform; bigger range and 
    more points makes a better numerical integration.
    Smaller q_min will increase reliable spin echo length range. 
    Rmax is the "radius" of the largest expected object and can be set elsewhere.
    q_max is determined by the acceptance angle of the SESANS instrument.
    """
    from sas.sascalc.data_util.nxsunit import Converter

    q_min = dq = 0.1 * 2*pi / Rmax
    return np.arange(q_min,
                     Converter(q_max[1])(q_max[0],
                                         units="1/A"),
                     dq)

def Hankelconstructor(data):
    Rmax = 1000000
    #Rmax = #value in text box?
    q_calc = make_q(data.sample.zacceptance, Rmax)
    SElength = Converter(data._xunit)(data.x, "A")
    dq = q_calc[1] - q_calc[0]
    H0 = dq / (2 * pi) * q_calc
    repSE, repq = np.meshgrid(SElength,q_calc)
    repq=np.array(repq,dtype='f')
    repSE=np.array(repSE,dtype='f')
    H = dq / (2 * pi) * j0(repSE*repq)*repq
    return H0, H, q_calc

def hankeltrafo(H0, H, Iq_calc):
    G0 = np.dot(H0, Iq_calc)
    G = np.dot(H.T, Iq_calc)
    P = G - G0
    return P  # This is the logarithmic Polarization, not the linear one as found in Andersson2008!


class SesansTransform(object):
    #: Set of spin-echo lengths in the measured data
    SElength = None  # type: np.ndarray
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
        q_min = dq = 0.1 * 2*pi / self.Rmax
        q_max = self.zaccept
        q=np.arange(q_min,q_max)
        self.q = q
        self.dq = dq

    def _set_hankel(self):
        #Rmax = #value in text box somewhere in FitPage?
        q = self.q
        dq = self.dq
        SElength = self.SE

        H0 = dq / (2 * pi) * q
        repSE, repq = np.meshgrid(SElength,q)
        repq=np.array(repq,dtype='f')
        repSE=np.array(repSE,dtype='f')
        H = dq / (2 * pi) * j0(repSE*repq)*repq

        self._H, self._H0 = H, H0

