import numpy as np
import scipy.stats

from sasmodels.weights import Dispersion as BaseDispersion

class Dispersion(BaseDispersion):
    r"""
    Gaussian dispersion, with 1-$\sigma$ width.

    .. math::

        w = \exp\left(-\tfrac12 (x - c)^2/\sigma^2\right)

    $x$ points are chosen such that each interval has equal weight.

    This works surprisingly poorly.  Try::

        $ sascomp cylinder -2d theta=45 phi=20 phi_pd_type=gaussian_eq \
          phi_pd_n=100,1000 radius=50 length=2*radius -midq phi_pd=5

    Leaving it here for others to improve.
    """
    type = "gaussian_eq"
    default = dict(npts=35, width=0, nsigmas=3)

    def _weights(self, center, sigma, lb, ub):
        # Use the gaussian distribution from scipy.stats
        dist = scipy.stats.norm(center, sigma)

        # Find the mid-points of the cdf intervals
        cdf = np.linspace(0, 1, self.npts+2)[1:-1]
        x = dist.ppf(cdf)

        # Since we are equally spaced in cdf, all weights are the same
        wx = np.ones_like(x)

        # Truncate the distribution in case the parameter value is limited
        index = (x >= lb) & (x <= ub)
        x, wx = x[index], wx[index]

        return x, wx

