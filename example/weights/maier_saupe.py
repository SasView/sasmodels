from __future__ import print_function
import numpy as np
from numpy import exp, sin, degrees, radians, pi, sqrt

from sasmodels.weights import Dispersion as BaseDispersion

class Dispersion(BaseDispersion):
    r"""
    Maier-Saupe dispersion on orientation.

    .. math:

        w(\theta) = e^{P_2{\cos^2 \theta}}

    This provides a close match to the gaussian distribution for
    low angles, but the tails are limited to $\pm 90^\circ$.  For $P_2 \ll 1$
    the distribution is approximately uniform.  The usual polar coordinate
    projection applies, with $\theta$ weights scaled by $\cos \theta$
    and $\phi$ weights unscaled.

    This is equivalent to a cyclic gaussian distribution
    $w(\theta) = e^{-sin^2(\theta)/(2\P_2^2)}.
    """
    type = "maier_saupe"
    default = dict(npts=35, width=1, nsigmas=None)

    # Note: center is always zero for orientation distributions
    def _weights(self, center, sigma, lb, ub):
        # use the width parameter as the value for Maier-Saupe P_2
        # and find the equivalent width sigma
        P2 = sigma
        sigma = 1./sqrt(2.*P2)

        # Limit width to +/- 90 degrees
        width = min(self.nsigmas*sigma, pi/2)
        x = np.linspace(-width, width, self.npts)

        # Truncate the distribution in case the parameter value is limited
        x[(x >= radians(lb)) & (x <= radians(ub))]

        # Return orientation in degrees with Maier-Saupe weights
        # Note: weights are normalized to sum to 1, so we can scale
        # by an arbitrary scale factor c = exp(m) to get:
        #     w = exp(m*cos(x)**2)/c = exp(-m*sin(x)**2)
        return degrees(x), exp(-P2*sin(x)**2)
