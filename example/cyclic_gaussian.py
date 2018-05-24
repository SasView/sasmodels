import numpy as np
from numpy import exp, sin, cos, pi, radians, degrees

from sasmodels.weights import Dispersion as BaseDispersion

class Dispersion(BaseDispersion):
    r"""
    Cyclic gaussian dispersion on orientation.
    
    .. math:

        w(\theta) = e^{-\frac{\sin^2 \theta}{2 \sigma^2}}

    This provides a close match to the gaussian distribution for
    low angles (with $\sin \theta \approx \theta$), but the tails 
    are limited to $\pm 90^\circ$.  For $\sigma$ large the
    distribution is approximately uniform.  The usual polar coordinate
    projection applies, with $\theta$ weights scaled by $\cos \theta$
    and $\phi$ weights unscaled.

    This is closely related to a Maier-Saupe distribution with order
    parameter $P_2$ and appropriate scaling constants, and changes
    between $\sin$ and $\cos$ as appropriate for the coordinate system
    representation.
    """
    type = "cyclic_gaussian"
    default = dict(npts=35, width=1, nsigmas=3)

    # Note: center is always zero for orientation distributions
    def _weights(self, center, sigma, lb, ub):
        # Convert sigma in degrees to the approximately equivalent Maier-Saupe "a"
        sigma = radians(sigma)
        a = -0.5/sigma**2

        # Limit width to +/-90 degrees; use an open interval since the
        # pattern at +90 is the same as that at -90.
        width = min(self.nsigmas*sigma, pi/2)
        x = np.linspace(-width, width, self.npts+2)[1:-1]
        wx = P(x, a)

        # Return orientation in degrees with Maier-Saupe weights
        return degrees(x), exp(a*sin(x)**2)

