import numpy as np
from numpy import exp, sin, cos, pi, radians, degrees

from sasmodels.weights import Dispersion as BaseDispersion

class Dispersion(BaseDispersion):
    r"""
    Cyclic gaussian dispersion on orientation.

    .. math:

        w(\theta) = e^{-\frac{\sin^2 \theta}{2 \sigma^2}}

    This provides a close match to the gaussian distribution for
    low angles, but the tails are limited to $\pm 90^\circ$.  For $\sigma$
    large the distribution is approximately uniform.  The usual polar coordinate
    projection applies, with $\theta$ weights scaled by $\cos \theta$
    and $\phi$ weights unscaled.

    This is eqivalent to a Maier-Saupe distribution with order
    parameter $a = 1/(2 \sigma^2)$, with $\sigma$ in radians.
    """
    type = "cyclic_gaussian"
    default = dict(npts=35, width=1, nsigmas=3)

    # Note: center is always zero for orientation distributions
    def _weights(self, center, sigma, lb, ub):
        # Convert sigma in degrees to radians
        sigma = radians(sigma)

        # Limit width to +/- 90 degrees
        width = min(self.nsigmas*sigma, pi/2)
        x = np.linspace(-width, width, self.npts)

        # Truncate the distribution in case the parameter value is limited
        x[(x >= radians(lb)) & (x <= radians(ub))]

        # Return orientation in degrees with Maier-Saupe weights
        return degrees(x), exp(-0.5*sin(x)**2/sigma**2)
