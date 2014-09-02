"""
SAS distributions for polydispersity.
"""
# TODO: include dispersion docs with the disperser models
from math import sqrt
import numpy as np
from scipy.special import gammaln

class Dispersion(object):
    """
    Base dispersion object.

    Subclasses should define *_weights(center, sigma, lb, ub)*
    which returns the x points and their corresponding weights.
    """
    type = "base disperser"
    default = dict(npts=35, width=0, nsigmas=3)
    def __init__(self, npts=None, width=None, nsigmas=None):
        self.npts = self.default['npts'] if npts is None else npts
        self.width = self.default['width'] if width is None else width
        self.nsigmas = self.default['nsigmas'] if nsigmas is None else nsigmas

    def get_pars(self):
        pars = {'type': self.type}
        pars.update(self.__dict__)
        return pars

    def set_weights(self, values, weights):
        raise RuntimeError("set_weights is only available for ArrayDispersion")

    def get_weights(self, center, lb, ub, relative):
        """
        Return the weights for the distribution.

        *center* is the center of the distribution

        *lb*,*ub* are the min and max allowed values

        *relative* is True if the distribution width is proportional to the
        center value instead of absolute.  For polydispersity use relative.
        For orientation parameters use absolute.
        """
        sigma = self.width * center if relative else self.width
        if sigma == 0 or self.npts < 2:
            if lb <= center <= ub:
                return np.array([center], 'd'), np.array([1.], 'd')
            else:
                return np.array([], 'd'), np.array([], 'd')
        return self._weights(center, sigma, lb, ub)

    def _weights(self, center, sigma, lb, ub):
        """actual work of computing the weights"""
        raise NotImplementedError

    def _linspace(self, center, sigma, lb, ub):
        """helper function to provide linear spaced weight points within range"""
        npts, nsigmas = self.npts, self.nsigmas
        x = center + np.linspace(-nsigmas*sigma, +nsigmas*sigma, npts)
        x = x[(x >= lb) & (x <= ub)]
        return x


class GaussianDispersion(Dispersion):
    type = "gaussian"
    default = dict(npts=35, width=0, nsigmas=3)
    def _weights(self, center, sigma, lb, ub):
        x = self._linspace(center, sigma, lb, ub)
        px = np.exp((x-center)**2 / (-2.0 * sigma * sigma))
        return x, px


class RectangleDispersion(Dispersion):
    type = "rectangle"
    default = dict(npts=35, width=0, nsigmas=1.70325)
    def _weights(self, center, sigma, lb, ub):
        x = self._linspace(center, sigma*sqrt(3.0), lb, ub)
        px = np.ones_like(x)
        return x, px


class LogNormalDispersion(Dispersion):
    type = "lognormal"
    default = dict(npts=80, width=0, nsigmas=8)
    def _weights(self, center, sigma, lb, ub):
        x = self._linspace(center, sigma, max(lb,1e-8), max(ub,1e-8))
        px = np.exp(-0.5*(np.log(x)-center)**2)/sigma**2/(x*sigma)
        return x, px


class SchulzDispersion(Dispersion):
    type = "schulz"
    default = dict(npts=80, width=0, nsigmas=8)
    def _weights(self, center, sigma, lb, ub):
        x = self._linspace(center, sigma, max(lb,1e-8), max(ub,1e-8))
        R= x/center
        z = (center/sigma)**2
        arg = z*np.log(z) + (z-1)*np.log(R) - R*z - np.log(center) - gammaln(z)
        px = np.exp(arg)
        return x, px


class ArrayDispersion(Dispersion):
    type = "array"
    default = dict(npts=35, width=0, nsigmas=1)
    def __init__(self, npts=None, width=None, nsigmas=None):
        Dispersion.__init__(self, npts, width, nsigmas)
        self.values = np.array([0.], 'd')
        self.weights = np.array([1.], 'd')

    def set_weights(self, values, weights):
        self.values = np.ascontiguousarray(values, 'd')
        self.weights = np.ascontiguousarray(weights, 'd')
        self.npts = len(values)

    def _weights(self, center, sigma, lb, ub):
        # TODO: interpolate into the array dispersion using npts
        x = center + self.values*sigma
        idx = (x>=lb)&(x<=ub)
        x = x[idx]
        px = self.weights[idx]
        return x, px


# dispersion name -> disperser lookup table.
models = dict((d.type,d) for d in (
    GaussianDispersion, RectangleDispersion,
    ArrayDispersion, SchulzDispersion, LogNormalDispersion
))


def get_weights(disperser, n, width, nsigmas, value, limits, relative):
    """
    Return the set of values and weights for a polydisperse parameter.

    *disperser* is the name of the disperser.

    *n* is the number of points in the weight vector.

    *width* is the width of the disperser distribution.

    *nsigmas* is the number of sigmas to span for the dispersion convolution.

    *value* is the value of the parameter in the model.

    *limits* is [lb, ub], the lower and upper bound of the weight value.

    *relative* is true if *width* is defined in proportion to the value
    of the parameter, and false if it is an absolute width.

    Returns *(value,weight)*, where *value* and *weight* are vectors.
    """
    cls = models[disperser]
    obj = cls(n, width, nsigmas)
    v,w = obj.get_weights(value, limits[0], limits[1], relative)
    return v,w

# Hack to allow sasview dispersion objects to interoperate with sasmodels
dispersers = dict((v.__name__,k) for k,v in models.items())
dispersers['DispersionModel'] = RectangleDispersion.type

