import numpy as np
from numpy import exp, sin, degrees, radians, pi, sqrt

from sasmodels.weights import Dispersion as BaseDispersion

class Dispersion(BaseDispersion):
    r"""
    Maier-Saupe dispersion on orientation (equal weights).

    .. math:

        w(\theta) = e^{a{\cos^2 \theta}}

    This provides a close match to the gaussian distribution for
    low angles, but the tails are limited to $\pm 90^\circ$.  For $a \ll 1$
    the distribution is approximately uniform.  The usual polar coordinate
    projection applies, with $\theta$ weights scaled by $\cos \theta$
    and $\phi$ weights unscaled.

    This is equivalent to a cyclic gaussian distribution
    $w(\theta) = e^{-sin^2(\theta)/(2\sigma^2)}.

    The $\theta$ points are spaced such that each interval has an
    equal contribution to the distribution.

    This works surprisingly poorly.  Try::

        $ sascomp cylinder -2d theta=45 phi=20 phi_pd_type=maier_saupe_eq \
          phi_pd_n=100,1000 radius=50 length=2*radius -midq phi_pd=5

    Leaving it here for others to improve.
    """
    type = "maier_saupe_eq"
    default = dict(npts=35, width=1, nsigmas=None)

    # Note: center is always zero for orientation distributions
    def _weights(self, center, sigma, lb, ub):
        # use the width parameter as the value for Maier-Saupe "a"
        a = sigma
        sigma = 1./sqrt(2.*a)

        # Create a lookup table for finding n points equally spaced
        # in the cumulative density function.
        # Limit width to +/-90 degrees.
        width = min(self.nsigmas*sigma, pi/2)
        xp = np.linspace(-width, width, max(self.npts*10, 100))

        # Compute CDF. Since we normalized the sum of the weights to 1,
        # we can scale by an arbitrary scale factor c = exp(m) to get:
        #     w = exp(m*cos(x)**2)/c = exp(-m*sin(x)**2)
        yp = np.cumsum(exp(-a*sin(xp)**2))
        yp /= yp[-1]

        # Find the mid-points of the equal-weighted intervals in the CDF
        y = np.linspace(0, 1, self.npts+2)[1:-1]
        x = np.interp(y, yp, xp)
        wx = np.ones(self.npts)

        # Truncate the distribution in case the parameter value is limited
        index = (x >= lb) & (x <= ub)
        x, wx = x[index], wx[index]

        return degrees(x), wx
