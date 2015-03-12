"""
Wrap sasmodels for direct use by bumps.
"""

import datetime

import numpy as np

from . import sesans

# CRUFT python 2.6
if not hasattr(datetime.timedelta, 'total_seconds'):
    def delay(dt):
        """Return number date-time delta as number seconds"""
        return dt.days * 86400 + dt.seconds + 1e-6 * dt.microseconds
else:
    def delay(dt):
        """Return number date-time delta as number seconds"""
        return dt.total_seconds()


def tic():
    """
    Timer function.

    Use "toc=tic()" to start the clock and "toc()" to measure
    a time interval.
    """
    then = datetime.datetime.now()
    return lambda: delay(datetime.datetime.now() - then)


def load_data(filename):
    """
    Load data using a sasview loader.
    """
    from sas.dataloader.loader import Loader
    loader = Loader()
    data = loader.load(filename)
    if data is None:
        raise IOError("Data %r could not be loaded" % filename)
    return data


def empty_data1D(q):
    """
    Create empty 1D data using the given *q* as the x value.

    Resolutions dq/q is 5%.
    """

    from sas.dataloader.data_info import Data1D

    Iq = 100 * np.ones_like(q)
    dIq = np.sqrt(Iq)
    data = Data1D(q, Iq, dx=0.05 * q, dy=dIq)
    data.filename = "fake data"
    data.qmin, data.qmax = q.min(), q.max()
    return data


def empty_data2D(qx, qy=None):
    """
    Create empty 2D data using the given mesh.

    If *qy* is missing, create a square mesh with *qy=qx*.

    Resolution dq/q is 5%.
    """
    from sas.dataloader.data_info import Data2D, Detector

    if qy is None:
        qy = qx
    Qx, Qy = np.meshgrid(qx, qy)
    Qx, Qy = Qx.flatten(), Qy.flatten()
    Iq = 100 * np.ones_like(Qx)
    dIq = np.sqrt(Iq)
    mask = np.ones(len(Iq), dtype='bool')

    data = Data2D()
    data.filename = "fake data"
    data.qx_data = Qx
    data.qy_data = Qy
    data.data = Iq
    data.err_data = dIq
    data.mask = mask

    # 5% dQ/Q resolution
    data.dqx_data = 0.05 * Qx
    data.dqy_data = 0.05 * Qy

    detector = Detector()
    detector.pixel_size.x = 5 # mm
    detector.pixel_size.y = 5 # mm
    detector.distance = 4 # m
    data.detector.append(detector)
    data.xbins = qx
    data.ybins = qy
    data.source.wavelength = 5 # angstroms
    data.source.wavelength_unit = "A"
    data.Q_unit = "1/A"
    data.I_unit = "1/cm"
    data.q_data = np.sqrt(Qx ** 2 + Qy ** 2)
    data.xaxis("Q_x", "A^{-1}")
    data.yaxis("Q_y", "A^{-1}")
    data.zaxis("Intensity", r"\text{cm}^{-1}")
    return data


def set_beam_stop(data, radius, outer=None):
    """
    Add a beam stop of the given *radius*.  If *outer*, make an annulus.
    """
    from sas.dataloader.manipulations import Ringcut
    if hasattr(data, 'qx_data'):
        data.mask = Ringcut(0, radius)(data)
        if outer is not None:
            data.mask += Ringcut(outer, np.inf)(data)
    else:
        data.mask = (data.x >= radius)
        if outer is not None:
            data.mask &= (data.x < outer)


def set_half(data, half):
    """
    Select half of the data, either "right" or "left".
    """
    from sas.dataloader.manipulations import Boxcut
    if half == 'right':
        data.mask += \
            Boxcut(x_min=-np.inf, x_max=0.0, y_min=-np.inf, y_max=np.inf)(data)
    if half == 'left':
        data.mask += \
            Boxcut(x_min=0.0, x_max=np.inf, y_min=-np.inf, y_max=np.inf)(data)


def set_top(data, cutoff):
    """
    Chop the top off the data, above *cutoff*.
    """
    from sas.dataloader.manipulations import Boxcut
    data.mask += \
        Boxcut(x_min=-np.inf, x_max=np.inf, y_min=-np.inf, y_max=cutoff)(data)


def plot_data(data, Iq, vmin=None, vmax=None, view='log'):
    """
    Plot the target value for the data.  This could be the data itself,
    the theory calculation, or the residuals.

    *scale* can be 'log' for log scale data, or 'linear'.
    """
    from numpy.ma import masked_array
    import matplotlib.pyplot as plt
    if hasattr(data, 'qx_data'):
        Iq = Iq + 0
        valid = np.isfinite(Iq)
        if view == 'log':
            valid[valid] = (Iq[valid] > 0)
            Iq[valid] = np.log10(Iq[valid])
        elif view == 'q4':
            Iq[valid] = Iq*(data.qx_data[valid]**2+data.qy_data[valid]**2)**2
        Iq[~valid | data.mask] = 0
        #plottable = Iq
        plottable = masked_array(Iq, ~valid | data.mask)
        xmin, xmax = min(data.qx_data), max(data.qx_data)
        ymin, ymax = min(data.qy_data), max(data.qy_data)
        try:
            if vmin is None: vmin = Iq[valid & ~data.mask].min()
            if vmax is None: vmax = Iq[valid & ~data.mask].max()
        except:
            vmin, vmax = 0, 1
        plt.imshow(plottable.reshape(128, 128),
                   interpolation='nearest', aspect=1, origin='upper',
                   extent=[xmin, xmax, ymin, ymax], vmin=vmin, vmax=vmax)
    else: # 1D data
        if view == 'linear' or view == 'q4':
            #idx = np.isfinite(Iq)
            scale = data.x**4 if view == 'q4' else 1.0
            plt.plot(data.x, scale*Iq) #, '.')
        else:
            # Find the values that are finite and positive
            idx = np.isfinite(Iq)
            idx[idx] = (Iq[idx] > 0)
            Iq[~idx] = np.nan
            plt.loglog(data.x, Iq)


def _plot_result1D(data, theory, view):
    """
    Plot the data and residuals for 1D data.
    """
    import matplotlib.pyplot as plt
    from numpy.ma import masked_array, masked
    #print "not a number",sum(np.isnan(data.y))
    #data.y[data.y<0.05] = 0.5
    mdata = masked_array(data.y, data.mask)
    mdata[np.isnan(mdata)] = masked
    if view is 'log':
        mdata[mdata <= 0] = masked
    mtheory = masked_array(theory, mdata.mask)
    mresid = masked_array((theory - data.y) / data.dy, mdata.mask)

    scale = data.x**4 if view == 'q4' else 1.0
    plt.subplot(121)
    plt.errorbar(data.x, scale*mdata, yerr=data.dy)
    plt.plot(data.x, scale*mtheory, '-', hold=True)
    plt.yscale('linear' if view == 'q4' else view)
    plt.subplot(122)
    plt.plot(data.x, mresid, 'x')

# pylint: disable=unused-argument
def _plot_sesans(data, theory, view):
    import matplotlib.pyplot as plt
    resid = (theory - data.y) / data.dy
    plt.subplot(121)
    plt.errorbar(data.x, data.y, yerr=data.dy)
    plt.plot(data.x, theory, '-', hold=True)
    plt.xlabel('spin echo length (A)')
    plt.ylabel('polarization')
    plt.subplot(122)
    plt.plot(data.x, resid, 'x')
    plt.xlabel('spin echo length (A)')
    plt.ylabel('residuals')

def _plot_result2D(data, theory, view):
    """
    Plot the data and residuals for 2D data.
    """
    import matplotlib.pyplot as plt
    resid = (theory - data.data) / data.err_data
    plt.subplot(131)
    plot_data(data, data.data, view=view)
    plt.colorbar()
    plt.subplot(132)
    plot_data(data, theory, view=view)
    plt.colorbar()
    plt.subplot(133)
    plot_data(data, resid, view='linear')
    plt.colorbar()

class BumpsModel(object):
    """
    Return a bumps wrapper for a SAS model.

    *data* is the data to be fitted.

    *model* is the SAS model from :func:`core.load_model`.

    *cutoff* is the integration cutoff, which avoids computing the
    the SAS model where the polydispersity weight is low.

    Model parameters can be initialized with additional keyword
    arguments, or by assigning to model.parameter_name.value.

    The resulting bumps model can be used directly in a FitProblem call.
    """
    def __init__(self, data, model, cutoff=1e-5, **kw):
        from bumps.names import Parameter

        # remember inputs so we can inspect from outside
        self.data = data
        self.model = model
        self.cutoff = cutoff
        if hasattr(data, 'lam'):
            self.data_type = 'sesans'
        elif hasattr(data, 'qx_data'):
            self.data_type = 'Iqxy'
        else:
            self.data_type = 'Iq'

        partype = model.info['partype']

        # interpret data
        if self.data_type == 'sesans':
            q = sesans.make_q(data.sample.zacceptance, data.Rmax)
            self.index = slice(None, None)
            self.Iq = data.y
            self.dIq = data.dy
            self._theory = np.zeros_like(q)
            q_vectors = [q]
        elif self.data_type == 'Iqxy':
            self.index = (data.mask == 0) & (~np.isnan(data.data))
            self.Iq = data.data[self.index]
            self.dIq = data.err_data[self.index]
            self._theory = np.zeros_like(data.data)
            if not partype['orientation'] and not partype['magnetic']:
                q_vectors = [np.sqrt(data.qx_data ** 2 + data.qy_data ** 2)]
            else:
                q_vectors = [data.qx_data, data.qy_data]
        elif self.data_type == 'Iq':
            self.index = (data.x >= data.qmin) & (data.x <= data.qmax) & ~np.isnan(data.y)
            self.Iq = data.y[self.index]
            self.dIq = data.dy[self.index]
            self._theory = np.zeros_like(data.y)
            q_vectors = [data.x]
        else:
            raise ValueError("Unknown data type") # never gets here

        # Remember function inputs so we can delay loading the function and
        # so we can save/restore state
        self._fn_inputs = [v[self.index] for v in q_vectors]
        self._fn = None

        # define bumps parameters
        pars = []
        for p in model.info['parameters']:
            name, default, limits = p[0], p[2], p[3]
            value = kw.pop(name, default)
            setattr(self, name, Parameter.default(value, name=name, limits=limits))
            pars.append(name)
        for name in partype['pd-2d']:
            for xpart, xdefault, xlimits in [
                    ('_pd', 0, limits),
                    ('_pd_n', 35, (0, 1000)),
                    ('_pd_nsigma', 3, (0, 10)),
                    ('_pd_type', 'gaussian', None),
                ]:
                xname = name + xpart
                xvalue = kw.pop(xname, xdefault)
                if xlimits is not None:
                    xvalue = Parameter.default(xvalue, name=xname, limits=xlimits)
                    pars.append(xname)
                setattr(self, xname, xvalue)
        self._parameter_names = pars
        if kw:
            raise TypeError("unexpected parameters: %s"
                            % (", ".join(sorted(kw.keys()))))
        self.update()

    def update(self):
        self._cache = {}

    def numpoints(self):
        """
            Return the number of points
        """
        return len(self.Iq)

    def parameters(self):
        """
            Return a dictionary of parameters
        """
        return dict((k, getattr(self, k)) for k in self._parameter_names)

    def theory(self):
        if 'theory' not in self._cache:
            if self._fn is None:
                q_input = self.model.make_input(self._fn_inputs)
                self._fn = self.model(q_input)

            fixed_pars = [getattr(self, p).value for p in self._fn.fixed_pars]
            pd_pars = [self._get_weights(p) for p in self._fn.pd_pars]
            #print fixed_pars,pd_pars
            self._theory[self.index] = self._fn(fixed_pars, pd_pars, self.cutoff)
            #self._theory[:] = self._fn.eval(pars, pd_pars)
            if self.data_type == 'sesans':
                result = sesans.hankel(self.data.x, self.data.lam * 1e-9,
                                       self.data.sample.thickness / 10,
                                       self._fn_inputs[0], self._theory)
                self._cache['theory'] = result
            else:
                self._cache['theory'] = self._theory
        return self._cache['theory']

    def residuals(self):
        #if np.any(self.err ==0): print "zeros in err"
        return (self.theory()[self.index] - self.Iq) / self.dIq

    def nllf(self):
        delta = self.residuals()
        #if np.any(np.isnan(R)): print "NaN in residuals"
        return 0.5 * np.sum(delta ** 2)

    #def __call__(self):
    #    return 2 * self.nllf() / self.dof

    def plot(self, view='log'):
        """
        Plot the data and residuals.
        """
        data, theory = self.data, self.theory()
        if self.data_type == 'Iq':
            _plot_result1D(data, theory, view)
        elif self.data_type == 'Iqxy':
            _plot_result2D(data, theory, view)
        elif self.data_type == 'sesans':
            _plot_sesans(data, theory, view)
        else:
            raise ValueError("Unknown data type")

    def simulate_data(self, noise=None):
        print "noise", noise
        if noise is None:
            noise = self.dIq[self.index]
        else:
            noise = 0.01 * noise
            self.dIq[self.index] = noise
        y = self.theory()
        y += y * np.random.randn(*y.shape) * noise
        if self.data_type == 'Iq':
            self.data.y[self.index] = y
        elif self.data_type == 'Iqxy':
            self.data.data[self.index] = y
        elif self.data_type == 'sesans':
            self.data.y[self.index] = y
        else:
            raise ValueError("Unknown model")

    def save(self, basename):
        pass

    def _get_weights(self, par):
        """
        Get parameter dispersion weights
        """
        from . import weights

        relative = self.model.info['partype']['pd-rel']
        limits = self.model.info['limits']
        disperser, value, npts, width, nsigma = [
            getattr(self, par + ext)
            for ext in ('_pd_type', '', '_pd_n', '_pd', '_pd_nsigma')]
        value, weight = weights.get_weights(
            disperser, int(npts.value), width.value, nsigma.value,
            value.value, limits[par], par in relative)
        return value, weight / np.sum(weight)

    def __getstate__(self):
        # Can't pickle gpu functions, so instead make them lazy
        state = self.__dict__.copy()
        state['_fn'] = None
        return state

    def __setstate__(self, state):
        # pylint: disable=attribute-defined-outside-init
        self.__dict__ = state


def demo():
    data = load_data('DEC07086.DAT')
    set_beam_stop(data, 0.004)
    plot_data(data, data.data)
    import matplotlib.pyplot as plt; plt.show()


if __name__ == "__main__":
    demo()
