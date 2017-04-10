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

class SesansTransform(object):
    """
    Spin-Echo SANS transform calculator.  Similar to a resolution function,
    the SesansTransform object takes I(q) for the set of *q_calc* values and
    produces a transformed dataset.

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
    _H = None  # type: np.ndarray
    _H0 = None # type: np.ndarray

    def __init__(self, z, SElength, zaccept, Rmax):
        # type: (np.ndarray, float, float) -> None
        #import logging; logging.info("creating SESANS transform")
        self.q = z
        # isoriented flag determines whether data is from an oriented sample or
        # not, should be a selection variable upon entering SESANS data.
        self._set_hankel(SElength, zaccept, Rmax)

    def apply(self, Iq):
        G0 = np.dot(self._H0, Iq)
        G = np.dot(self._H.T, Iq)
        P = G - G0
        return P

    def _set_hankel(self, SElength, zaccept, Rmax):
        # type: (np.ndarray, float, float) -> None
        # Force float32 arrays, otherwise run into memory problems on some machines
        SElength = np.asarray(SElength, dtype='float32')

        #Rmax = #value in text box somewhere in FitPage?
        q_max = 2*pi / (SElength[1] - SElength[0])
        q_min = 0.1 * 2*pi / (np.size(SElength) * SElength[-1])
        q = np.arange(q_min, q_max, q_min, dtype='float32')
        dq = q_min

        H0 = np.float32(dq/(2*pi)) * q

        repq = np.tile(q, (SElength.size, 1)).T
        repSE = np.tile(SElength, (q.size, 1))
        H = np.float32(dq/(2*pi)) * j0(repSE*repq) * repq

        self.q_calc = q
        self._H, self._H0 = H, H0

class OrientedSesansTransform(object):
    """
    Oriented Spin-Echo SANS transform calculator.  Similar to a resolution
    function, the OrientedSesansTransform object takes I(q) for the set
    of *q_calc* values and produces a transformed dataset.

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
    _cosmat = None  # type: np.ndarray
    _cos0 = None # type: np.ndarray
    _Iq_shape = None # type: Tuple[int, int]

    def __init__(self, z, SElength, zaccept, Rmax):
        # type: (np.ndarray, float, float) -> None
        #import logging; logging.info("creating SESANS transform")
        self.q = z
        # isoriented flag determines whether data is from an oriented sample or
        # not, should be a selection variable upon entering SESANS data.
        self._set_cosmat(SElength, zaccept, Rmax)

    def apply(self, Iq):
        dq = self.q_calc[0][0]
        Iq = np.reshape(Iq, self._Iq_shape)
        G0 = self._cos0 * np.sum(Iq) * dq
        G = np.sum(np.dot(Iq, self._cosmat.T), axis=1) * dq
        P = G - G0
        return P

    def _set_cosmat(self, SElength, zaccept, Rmax):
        # type: (np.ndarray, float, float) -> None
        # Force float32 arrays, otherwise run into memory problems on some machines
        SElength = np.asarray(SElength, dtype='float32')

        # Rmax = #value in text box somewhere in FitPage?
        q_max = 2 * pi / (SElength[1] - SElength[0])
        q_min = 0.1 * 2 * pi / (np.size(SElength) * SElength[-1])

        q = np.arange(q_min, q_max, q_min, dtype='float32')
        dq = q_min

        cos0 = np.float32(dq / (2 * pi))
        cosmat = np.float32(dq / (2 * pi)) * np.cos(q[:, None] * SE[None, :])

        qx, qy = np.meshgrid(q, q)
        self._Iq_shape = qx.shape
        self.q_calc = qx.flatten(), qy.flatten()
        self._cosmat, self._cos0 = cosmat, cos0
