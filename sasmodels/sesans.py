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
from scipy.special import j1


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

    def __init__(self, z, SElength, lam, zaccept, Rmax):
        # type: (np.ndarray, float, float) -> None
        self.q = z
        self._set_hankel(SElength, lam, zaccept, Rmax)

    def apply(self, Iq):
        # tye: (np.ndarray) -> np.ndarray
        G0 = np.dot(self._H0, Iq)
        G = np.dot(self._H.T, Iq)
        G0 = G[0]  # FIXME: This is a kludge
        P = G - G0
        return P

    def _set_hankel(self, SElength, lam, zaccept, Rmax):
        # type: (np.ndarray, float, float) -> None
        # Force float32 arrays, otherwise run into memory problems on
        # some machines
        SElength = np.asarray(SElength, dtype='float32')

        # Rmax = #value in text box somewhere in FitPage?
        q_max = 2*pi / (SElength[1] - SElength[0])
        q_min = 0.1 * 2*pi / (np.size(SElength) * SElength[-1])
        # q = np.arange(q_min, q_max, q_min, dtype='float32')
        # q = np.exp(np.arange(np.log(q_min), np.log(q_max), np.log(2),
        #                      dtype=np.float32))
        q = np.exp(np.linspace(np.log(q_min), np.log(q_max), 10*SElength.size,
                               dtype=np.float32))
        q = np.hstack([[0], q])

        H0 = np.pi * (q[1:]**2 - q[:-1]**2) * (q[1:] - q[:-1])

        # repq = np.tile(q, (SElength.size, 1)).T
        H = np.outer(q, SElength)
        j1(H, out=H)
        H *= q.reshape((-1, 1))
        H = H[1:] - H[:-1]
        H *= 2 * np.pi / SElength

        lam = np.asarray(lam, dtype=np.float32)
        reptheta = np.outer(q[1:], lam)
        reptheta /= np.float32(2*np.pi)
        np.arcsin(reptheta, out=reptheta)
        # reptheta = np.arcsin(repq*replam/2*np.pi)
        mask = reptheta > zaccept
        # H[mask] = 0

        # H = np.zeros((q.size, SElength.size), dtype=np.float32)
        # H0 = q * 0
        assert(H.shape == (q.size-1, SElength.size))

        self.q_calc = q[1:]
        self._H, self._H0 = H, H0
