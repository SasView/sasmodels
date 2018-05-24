import numpy as np
from scipy.stats import laplace

from sasmodels import weights

class Dispersion(weights.Dispersion):
    r"""
    Laplace distribution

    .. math::

        w(x) = e^{-\sigma |x - \mu|}
    """
    type = "laplace"
    default = dict(npts=35, width=0, nsigmas=3)  # default values
    def _weights(self, center, sigma, lb, ub):
        x = self._linspace(center, sigma, lb, ub)
        wx = laplace.pdf(x, center, sigma)
        return x, wx

