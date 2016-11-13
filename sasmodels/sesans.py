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
from scipy.special import jv as besselj
from scipy.special import j0
#from mpmath import j0
from mpmath import besselj
from mpmath import mpf
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
    q_calc = make_q(data.sample.zacceptance, Rmax)
    SElength = Converter(data._xunit)(data.x, "A")
    dq = q_calc[1] - q_calc[0]
    H0 = dq / (2 * pi) * q_calc
    repSE, repq = np.meshgrid(SElength,q_calc)
    H = dq / (2 * pi) * j0(repSE*repq)*repq
    return H0, H, q_calc

def hankeltrafo(H0, H, Iq_calc):
    G0 = np.dot(H0, Iq_calc)
    G = np.dot(H.T, Iq_calc)
    P = G - G0
    return P  # This is the logarithmic Polarization, not the linear one as found in Andersson2008!
