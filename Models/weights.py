import numpy as np

class GaussianDispersion(object):
    def __init__(self, npts=35, width=0, nsigmas=3): #number want, percent deviation, #standard deviations from mean
        self.type = 'gaussian'
        self.npts = npts
        self.width = width
        self.nsigmas = nsigmas

    def get_pars(self):
        return self.__dict__

    def get_weights(self, center, min, max, relative):
        """ *center* is the center of the distribution
        *min*,*max* are the min, max allowed values
        *relative* is True if the width is relative to the center instead of absolute
        For polydispersity use relative.  For orientation parameters use absolute."""
        npts, width, nsigmas = self.npts, self.width, self.nsigmas
        sigma = width * center if relative else width
        if sigma == 0:
            return np.array([center],'d'), np.array([1.], 'd')
        x = center + np.linspace(-nsigmas * sigma, +nsigmas * sigma, npts)
        x = x[(x >= min) & (x <= max)]
        px = np.exp((x-center)**2 / (-2.0 * sigma * sigma))
        return x, px

