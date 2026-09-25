r"""
Special functions for implementing scattering models.

Documentation for the special function library is in the
:ref:`Special_Functions` section of the manual. These are all available
for import from `sasmodels.special` even though most of them are not listed
in the documentation for this module.

See section :ref:`Python_Functions` for differences from the C library.
"""
# pylint: disable=unused-import

from functools import cache

import numpy as np

# Functions to add to our standard set
# C99 standard math library functions
from numpy import cos, inf, nan, sin
from scipy.special import j1 as sas_J1

# Gaussians
try:
    from scipy.special import _fast_gl, roots_legendre
except ImportError:
    # CRUFT: use our copy of the new roots_legendre function for older scipy
    from ._fastgl import gauss_legendre as roots_legendre

# erf, erfc, tgamma, lgamma  **do not use**

# C99 standard math constants
M_PI, M_PI_2, M_PI_4, M_SQRT1_2, M_E = np.pi, np.pi/2, np.pi/4, np.sqrt(0.5), np.e
NAN = nan
INFINITY = inf

# non-standard constants
M_PI_180, M_4PI_3 = M_PI/180, 4*M_PI/3

# can't do SINCOS in python; use "s, c = SINCOS(x)" instead
def SINCOS(x):
    """return sin(x), cos(x)"""
    return sin(x), cos(x)
sincos = SINCOS

def square(x):
    """return x^2"""
    return x*x

def cube(x):
    """return x^3"""
    return x*x*x

def sas_sinx_x(x):
    """return sin(x)/x"""
    from numpy import sinc as _sinc
    return _sinc(x/M_PI)

def powr(x, y):
    """return x^y for x>0"""
    return x**y
def pown(x, n):
    """return x^n for n integer"""
    return x**n

FLOAT_SIZE = 8

def polevl(x, c, n):
    """return p(x) for polynomial p of degree n-1 with coefficients c"""
    return np.polyval(c[:n], x)

def p1evl(x, c, n):
    """return x^n + p(x) for polynomial p of degree n-1 with coefficients c"""
    return np.polyval(np.hstack(([1.], c))[:n], x)

def sas_Si(x):
    """return Si(x)"""
    from scipy.special import sici
    return sici(x)[0]

def sas_j1(x):
    """return j1(x)"""
    if np.isscalar(x):
        retvalue = (sin(x) - x*cos(x))/x**2 if x != 0. else 0.
    else:
        with np.errstate(all='ignore'):
            retvalue = (sin(x) - x*cos(x))/x**2
        retvalue[x == 0.] = 0.
    return retvalue

def sas_3j1x_x(x):
    """return 3*j1(x)/x"""
    if np.isscalar(x):
        retvalue = 3*(sin(x) - x*cos(x))/x**3 if x != 0. else 1.
    else:
        with np.errstate(all='ignore'):
            retvalue = 3*(sin(x) - x*cos(x))/x**3
        retvalue[x == 0.] = 1.
    return retvalue

def sas_2J1x_x(x):
    """return 2*J1(x)/x"""
    if np.isscalar(x):
        retvalue = 2*sas_J1(x)/x if x != 0 else 1.
    else:
        with np.errstate(all='ignore'):
            retvalue = 2*sas_J1(x)/x
        retvalue[x == 0] = 1.
    return retvalue



_cached_roots_legendre = cache(roots_legendre)
_ADAPTIVE_MAX_N = 100000
ADAPTIVE_MAX_76 = (_ADAPTIVE_MAX_N // 76)
ADAPTIVE_MAX_OUTER = (_ADAPTIVE_MAX_N / 500)

def gauss_weights(qr: float, outer_n: int) -> tuple[np.ndarray, np.ndarray]:
    """Select a Gauss‑Legendre rule based on *qr* and *outer_n*.

    Adaptive integration uses nested gaussian integration with order determined by the
    length scale qr, which drives the oscillation frequency of the scattering amplitude.
    Using polar coordinates, the oscillation with longitude appears to be independent of
    of latitude, so efficient spherical integration methods that sample roughly uniformly
    in solid angle don't perform as well as expected. The theta-phi uniform grid only
    uses 2x more points than the efficient grids, so the improvements from using
    Gauss-Legendre integration outweigh the increase in the number of integration points.

    The qr for the loop over phi may be different from that over theta. This is particularly
    the case for long rods, where the loop over phi depends on the cross section of the
    rod, so you can use far fewer phi points than the max dimension of the shape would
    suggest. Models such as the triaxial ellipsoid can relabel their axis so that c the
    long axis of the shape, allowing a smaller number of total evaluation points.

    If size of the grid is too large even after relabeling then the evaluation will be
    very slow. To prevent this, we pass the number of elements in the outer loop to
    limit the size of the inner loop, giving n_outer*n_inner < _ADAPTIVE_MAX_N. This works
    well enough for existing models, giving high accuracy results even for 20 micro shapes
    at q=1e-3 (the upper limit of qr for USANS). Even out to q=1e-1 the calculations are
    mostly within 10% of the full grid value (good enough to estimate slit resolution for
    USANS measurements) without being exceedingly slow.

    For a standard integral over a hemisphere the loops will look like::

        from sasmodels.special import SINCOS, M_PI, sqrt, gauss_weights, ADAPTIVE_MAX_OUTER

        # Calculations that are independent of q
        ... swap dimensions so that c is the long axis
        ... compute contrast, volume, radius_effective and other q-independent values ...
        ... max_outer = max dimension # for a box use max(sqrt(a**2 + b**2), c)/2 ...
        ... max_inner = max cross section # for a box use max(a, b)/2 ...

        z_outer, w_outer = gauss_weights(q * max_outer, ADAPTIVE_MAX_OUTER)
        outer_F1 = outer_F2 = 0.
        for j in range(len(z_outer)):
            # x in [-1, 1] => cos(theta) in [0, 1]
            cos_theta = 0.5*z_outer[j] + 0.5
            sin_theta = sqrt(1.0 - cos_theta*cos_theta)
            qc = q*cos_theta
            qab = q*sin_theta

            z_inner, w_inner = gauss_weights(qab * max_inner, outer.n)
            inner_F1 = inner_F2 = 0.
            for k in range(len(z_inner)):
                # x in [-1, 1] => phi in [0, 2 pi]
                phi = M_PI*(z_inner[k] + 1.0)
                sin_phi, cos_phi = SINCOS(phi)
                qa = qab * sin_phi
                qb = qab * cos_phi

                F = Fq_abc(qa, qb, qc)

                inner_F1 += w_inner[k] * F
                inner_F2 += w_inner[k] * F * F
            outer_F1 += w_outer[k] * inner_F1
            outer_F2 += w_outer[k] * inner_F2

        # correct [-1, 1] => [0, 1] for u and [-1, 1] => [0, 2 pi] for phi
        outer_F1 /= 2.0*M_PI
        outer_F2 /= 2.0*M_PI

        # s (cm) = contrast (1E-6 / Ang^2) * volume (Ang^3)
        s = 1e-2 * contrast * volume
        F1 = s * outer_F1
        F2 = s * s * outer_F2
        return F1, F2

    Over the whole sphere is similar, except it would use z_outer directly for cos_theta,
    and the final correction would be 1/pi

    Some models (barbell, capped cylinder, pringle) are not using spherical integration.
    Instead we find that limiting the outer integral to 76 points leads to okay results::

        n_outer = gauss_weights(qr_outer, ADAPTIVE_MAX_76)
    """
    # Note: these values were determined empirically.
    # Before changing, do extensive testing on a range of models with various
    # sizes (20 nm to 20 um) and aspect ratios (rods, disks and cubes). For the
    # C models this was done in sasmodels/explore/check_adaptive.py
    if qr < 10 or outer_n > _ADAPTIVE_MAX_N // 76:
        n = 20
    elif qr < 100 or outer_n > _ADAPTIVE_MAX_N // 500:
        n = 76
    elif qr < 750 or outer_n > _ADAPTIVE_MAX_N // 5000:
        n = 500
    else:
        n = 5000
    return _cached_roots_legendre(n)


# CRUFT: don't need gauss{20,76,150} anymore because roots_legendre is now fast
class Gauss:
    """Gauss-Legendre integration weights"""
    def __init__(self, n):
        z, w = _cached_roots_legendre(n)
        self.n = n
        self.w = w
        self.z = z

gauss20 = Gauss(20)
gauss76 = Gauss(76)
gauss150 = Gauss(150)
