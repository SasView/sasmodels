"""
SAS distributions for polydispersity.
"""
# TODO: include dispersion docs with the disperser models
from __future__ import division, print_function

from math import sqrt  # type: ignore
from collections import OrderedDict

import numpy as np  # type: ignore
from scipy.special import gammaln  # type: ignore

# TODO: include dispersion docs with the disperser models

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
        """
        Return the parameters to the disperser as a dictionary.
        """
        pars = {'type': self.type}
        pars.update(self.__dict__)
        return pars

    # pylint: disable=no-self-use
    def set_weights(self, values, weights):
        """
        Set the weights on the disperser if it is :class:`ArrayDispersion`.
        """
        raise RuntimeError("set_weights is only available for ArrayDispersion")

    def get_weights(self, center, lb, ub, relative):
        """
        Return the weights for the distribution.

        *center* is the center of the distribution

        *lb*, *ub* are the min and max allowed values

        *relative* is True if the distribution width is proportional to the
        center value instead of absolute.  For polydispersity use relative.
        For orientation parameters use absolute.
        """
        sigma = self.width * center if relative else self.width
        if not relative:
            # For orientation, the jitter is relative to 0 not the angle
            center = 0
            pass
        if sigma == 0 or self.npts < 2:
            if lb <= center <= ub:
                return np.array([center], 'd'), np.array([1.], 'd')
            else:
                return np.array([], 'd'), np.array([], 'd')
        x, px = self._weights(center, sigma, lb, ub)
        return x, px

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
    r"""
    Gaussian dispersion, with 1-$\sigma$ width.

    .. math::

        w = \exp\left(-\tfrac12 (x - c)^2/\sigma^2\right)
    """
    type = "gaussian"
    default = dict(npts=35, width=0, nsigmas=3)
    def _weights(self, center, sigma, lb, ub):
        # TODO: sample high probability regions more densely
        # i.e., step uniformly in cumulative density rather than x value
        # so weight = 1/Npts for all weights, but values are unevenly spaced
        x = self._linspace(center, sigma, lb, ub)
        px = np.exp((x-center)**2 / (-2.0 * sigma * sigma))
        return x, px

class UniformDispersion(Dispersion):
    r"""
    Uniform dispersion, with width $\sigma$.

    .. math::

        w = 1
    """
    type = "uniform"
    default = dict(npts=35, width=0, nsigmas=None)
    def _weights(self, center, sigma, lb, ub):
        x = np.linspace(center-sigma, center+sigma, self.npts)
        x = x[(x >= lb) & (x <= ub)]
        return x, np.ones_like(x)

class RectangleDispersion(Dispersion):
    r"""
    Uniform dispersion, with width $\sqrt{3}\sigma$.

    .. math::

        w = 1
    """
    type = "rectangle"
    default = dict(npts=35, width=0, nsigmas=1.73205)
    def _weights(self, center, sigma, lb, ub):
         x = self._linspace(center, sigma, lb, ub)
         x = x[np.fabs(x-center) <= np.fabs(sigma)*sqrt(3.0)]
         return x, np.ones_like(x)

class LogNormalDispersion(Dispersion):
    r"""
    log Gaussian dispersion, with 1-$\sigma$ width.

    .. math::

        w = \frac{\exp\left(-\tfrac12 (\ln x - c)^2/\sigma^2\right)}{x\sigma}
    """
    type = "lognormal"
    default = dict(npts=80, width=0, nsigmas=8)
    def _weights(self, center, sigma, lb, ub):
        x = self._linspace(center, sigma, max(lb, 1e-8), max(ub, 1e-8))
        # sigma in the lognormal function is in ln(R) space, thus needs converting
        sig = np.fabs(sigma/center)
        px = np.exp(-0.5*((np.log(x)-np.log(center))/sig)**2)/(x*sig)
        return x, px


class SchulzDispersion(Dispersion):
    r"""
    Schultz dispersion, with 1-$\sigma$ width.

    .. math::

        w = \frac{z^z\,R^{z-1}}{e^{Rz}\,c \Gamma(z)}

    where $c$ is the center of the distribution, $R = x/c$ and $z=(c/\sigma)^2$.

    This is evaluated using logarithms as

    .. math::

        w = \exp\left(z \ln z + (z-1)\ln R - Rz - \ln c - \ln \Gamma(z) \right)
    """
    type = "schulz"
    default = dict(npts=80, width=0, nsigmas=8)
    def _weights(self, center, sigma, lb, ub):
        x = self._linspace(center, sigma, max(lb, 1e-8), max(ub, 1e-8))
        R = x/center
        z = (center/sigma)**2
        arg = z*np.log(z) + (z-1)*np.log(R) - R*z - np.log(center) - gammaln(z)
        px = np.exp(arg)
        return x, px


class ArrayDispersion(Dispersion):
    r"""
    Empirical dispersion curve.

    Use :meth:`set_weights` to set $w = f(x)$.
    """
    type = "array"
    default = dict(npts=35, width=0, nsigmas=1)
    def __init__(self, npts=None, width=None, nsigmas=None):
        Dispersion.__init__(self, npts, width, nsigmas)
        self.values = np.array([0.], 'd')
        self.weights = np.array([1.], 'd')

    def set_weights(self, values, weights):
        """
        Set the weights for the given x values.
        """
        self.values = np.ascontiguousarray(values, 'd')
        self.weights = np.ascontiguousarray(weights, 'd')
        self.npts = len(values)

    def _weights(self, center, sigma, lb, ub):
        # TODO: rebin the array dispersion using npts
        # TODO: use a distribution that can be recentered and scaled
        x = self.values
        #x = center + self.values*sigma
        idx = (x >= lb) & (x <= ub)
        x = x[idx]
        px = self.weights[idx]
        return x, px

class BoltzmannDispersion(Dispersion):
    r"""
    Boltzmann dispersion, with $\sigma=k T/E$.

    .. math::

        w = \exp\left( -|x - c|/\sigma\right)
    """
    type = "boltzmann"
    default = dict(npts=35, width=0, nsigmas=3)
    def _weights(self, center, sigma, lb, ub):
        x = self._linspace(center, sigma, lb, ub)
        px = np.exp(-np.abs(x-center) / np.abs(sigma))
        return x, px

# dispersion name -> disperser lookup table.
# Maintain order since this is used by sasview GUI to order the options in
# the dispersion type combobox.
MODELS = OrderedDict((d.type, d) for d in (
    RectangleDispersion,
    UniformDispersion,
    ArrayDispersion,
    LogNormalDispersion,
    GaussianDispersion,
    SchulzDispersion,
    BoltzmannDispersion
))


def get_weights(disperser, n, width, nsigmas, value, limits, relative):
    """
    Return the set of values and weights for a polydisperse parameter.

    *disperser* is the name of the disperser.

    *n* is the number of points in the weight vector.

    *width* is the width of the disperser distribution.

    *nsigmas* is the number of sigmas to span for the dispersion convolution.

    *value* is the value of the parameter in the model.

    *limits* is [lb, ub], the lower and upper bound on the possible values.

    *relative* is true if *width* is defined in proportion to the value
    of the parameter, and false if it is an absolute width.

    Returns *(value, weight)*, where *value* and *weight* are vectors.
    """
    if disperser == "array":
        raise NotImplementedError("Don't handle arrays through get_weights;"
                                  " use values and weights directly")
    cls = MODELS[disperser]
    obj = cls(n, width, nsigmas)
    v, w = obj.get_weights(value, limits[0], limits[1], relative)
    return v, w/np.sum(w)


def plot_weights(model_info, mesh):
    # type: (ModelInfo, List[Tuple[float, np.ndarray, np.ndarray]]) -> None
    """
    Plot the weights returned by :func:`get_weights`.

    *model_info* defines model parameters, etc.

    *mesh* is a list of tuples containing (*value*, *dispersity*, *weights*)
    for each parameter, where (*dispersity*, *weights*) pairs are the
    distributions to be plotted.
    """
    import pylab

    if any(len(dispersity) > 1 for value, dispersity, weights in mesh):
        labels = [p.name for p in model_info.parameters.call_parameters]
        #pylab.interactive(True)
        pylab.figure()
        for (v, x, w), s in zip(mesh, labels):
            if len(x) > 1:
                pylab.plot(x, w, '-o', label=s)
        pylab.grid(True)
        pylab.legend()
        #pylab.show()
