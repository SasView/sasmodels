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
from numpy import pi  # type: ignore
from scipy.special import j0


class SesansTransform(object):
    """
    Spin-Echo SANS transform calculator.  Similar to a resolution function,
    the SesansTransform object takes I(q) for the set of *q_calc* values and
    produces a transformed dataset

    *SElength* (A) is the set of spin-echo lengths in the measured data.

    *zaccept* (1/A) is the maximum acceptance of scattering vector in the spin
    echo encoding dimension (for ToF: Q of min(R) and max(lam)).

    *Rmax* (A) is the maximum size sensitivity; larger radius requires more
    computation time.
    """
    #: SElength from the data in the original data units; not used by transform
    #: but the GUI uses it, so make sure that it is present.
    q = None  # type: np.ndarray

    #: q values to calculate when computing transform
    q_calc = None  # type: np.ndarray

    # transform arrays
    _H = None   # type: np.ndarray
    _H0 = None  # type: np.ndarray

    def __init__(self, z, SElength, lam, zaccept, Rmax, log_spacing=1.0003):
        # type: (np.ndarray, float, float, float, float, float) -> None
        self.q = z
        self.log_spacing = log_spacing
        self._set_hankel(SElength, lam, zaccept, Rmax)

    def apply(self, Iq):
        # type: (np.ndarray) -> np.ndarray
        """
        Apply the SESANS transform to the computed I(q).
        """
        G0 = np.dot(self._H0, Iq)
        G = np.dot(self._H.T, Iq)
        P = G - G0
        return P

    def _set_hankel(self, SElength, lam, zaccept, Rmax):
        # type: (np.ndarray, float, float, float) -> None
        SElength = np.asarray(SElength)
        q_max = 2*pi / (SElength[1] - SElength[0])
        q_min = 0.1 * 2*pi / (np.size(SElength) * SElength[-1])
        q = np.exp(np.arange(np.log(q_min), np.log(q_max),
                             np.log(self.log_spacing)))

        dq = np.diff(q)
        dq = np.insert(dq, 0, dq[0])

        H0 = dq/(2*pi) * q

        H = np.outer(q, SElength)
        j0(H, out=H)
        H *= (dq * q / (2*pi)).reshape((-1, 1))

        reptheta = np.outer(q, lam/(2*pi))
        np.arcsin(reptheta, out=reptheta)
        mask = reptheta > zaccept
        H[mask] = 0

        self.q_calc = q
        self._H, self._H0 = H, H0
