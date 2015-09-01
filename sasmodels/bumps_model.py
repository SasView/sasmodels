"""
Wrap sasmodels for direct use by bumps.

:class:`Model` is a wrapper for the sasmodels kernel which defines a
bumps *Parameter* box for each kernel parameter.  *Model* accepts keyword
arguments to set the initial value for each parameter.

:class:`Experiment` combines the *Model* function with a data file loaded by the
sasview data loader.  *Experiment* takes a *cutoff* parameter controlling
how far the polydispersity integral extends.

A variety of helper functions are provided:

    :func:`load_data` loads a sasview data file.

    :func:`empty_data1D` creates an empty dataset, which is useful for plotting
    a theory function before the data is measured.

    :func:`empty_data2D` creates an empty 2D dataset.

    :func:`set_beam_stop` masks the beam stop from the data.

    :func:`set_half` selects the right or left half of the data, which can
    be useful for shear measurements which have not been properly corrected
    for path length and reflections.

    :func:`set_top` cuts the top part off the data.

    :func:`plot_data` plots the data file.

    :func:`plot_theory` plots a calculated result from the model.

"""

import datetime
import warnings

import numpy as np

from . import sesans
from .resolution import Perfect1D, Pinhole1D, Slit1D
from .resolution2d import Pinhole2D

# CRUFT python 2.6
if not hasattr(datetime.timedelta, 'total_seconds'):
    def delay(dt):
        """Return number date-time delta as number seconds"""
        return dt.days * 86400 + dt.seconds + 1e-6 * dt.microseconds
else:
    def delay(dt):
        """Return number date-time delta as number seconds"""
        return dt.total_seconds()


# CRUFT: old style bumps wrapper which doesn't separate data and model
def BumpsModel(data, model, cutoff=1e-5, **kw):
    warnings.warn("Use of BumpsModel is deprecated.  Use bumps_model.Experiment instead.")
    model = Model(model, **kw)
    experiment = Experiment(data=data, model=model, cutoff=cutoff)
    for k in model._parameter_names:
        setattr(experiment, k, getattr(model, k))
    return experiment



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

def plot_data(data, view='log'):
    """
    Plot data loaded by the sasview loader.
    """
    if hasattr(data, 'qx_data'):
        _plot_2d_signal(data, data.data, view=view)
    else:
        # Note: kind of weird using the _plot_result1D to plot just the
        # data, but it handles the masking and graph markup already, so
        # do not repeat.
        _plot_result1D(data, None, None, view)

def plot_theory(data, theory, view='log'):
    if hasattr(data, 'qx_data'):
        _plot_2d_signal(data, theory, view=view)
    else:
        _plot_result1D(data, theory, None, view, include_data=False)


def empty_data1D(q, resolution=0.05):
    """
    Create empty 1D data using the given *q* as the x value.

    *resolution* dq/q defaults to 5%.
    """

    from sas.dataloader.data_info import Data1D

    Iq = 100 * np.ones_like(q)
    dIq = np.sqrt(Iq)
    data = Data1D(q, Iq, dx=resolution * q, dy=dIq)
    data.filename = "fake data"
    data.qmin, data.qmax = q.min(), q.max()
    data.mask = np.zeros(len(Iq), dtype='bool')
    return data


def empty_data2D(qx, qy=None, resolution=0.05):
    """
    Create empty 2D data using the given mesh.

    If *qy* is missing, create a square mesh with *qy=qx*.

    *resolution* dq/q defaults to 5%.
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
    data.qmin = 1e-16
    data.qmax = np.inf

    # 5% dQ/Q resolution
    if resolution != 0:
        # https://www.ncnr.nist.gov/staff/hammouda/distance_learning/chapter_15.pdf
        # Should have an additional constant which depends on distances and
        # radii of the aperture, pixel dimensions and wavelength spread
        # Instead, assume radial dQ/Q is constant, and perpendicular matches
        # radial (which instead it should be inverse).
        Q = np.sqrt(Qx**2 + Qy**2)
        data.dqx_data = resolution * Q
        data.dqy_data = resolution * Q

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


def _plot_result1D(data, theory, resid, view, include_data=True):
    """
    Plot the data and residuals for 1D data.
    """
    import matplotlib.pyplot as plt
    from numpy.ma import masked_array, masked
    #print "not a number",sum(np.isnan(data.y))
    #data.y[data.y<0.05] = 0.5
    mdata = masked_array(data.y, data.mask)
    mdata[~np.isfinite(mdata)] = masked
    if view is 'log':
        mdata[mdata <= 0] = masked

    scale = data.x**4 if view == 'q4' else 1.0
    if resid is not None:
        plt.subplot(121)

    if include_data:
        plt.errorbar(data.x, scale*mdata, yerr=data.dy, fmt='.')
    if theory is not None:
        mtheory = masked_array(theory, mdata.mask)
        plt.plot(data.x, scale*mtheory, '-', hold=True)
    plt.xscale(view)
    plt.yscale('linear' if view == 'q4' else view)
    plt.xlabel('Q')
    plt.ylabel('I(Q)')
    if resid is not None:
        mresid = masked_array(resid, mdata.mask)
        plt.subplot(122)
        plt.plot(data.x, mresid, 'x')
        plt.ylabel('residuals')
        plt.xlabel('Q')
        plt.xscale(view)

# pylint: disable=unused-argument
def _plot_sesans(data, theory, resid, view):
    import matplotlib.pyplot as plt
    plt.subplot(121)
    plt.errorbar(data.x, data.y, yerr=data.dy)
    plt.plot(data.x, theory, '-', hold=True)
    plt.xlabel('spin echo length (nm)')
    plt.ylabel('polarization (P/P0)')
    plt.subplot(122)
    plt.plot(data.x, resid, 'x')
    plt.xlabel('spin echo length (nm)')
    plt.ylabel('residuals (P/P0)')

def _plot_result2D(data, theory, resid, view):
    """
    Plot the data and residuals for 2D data.
    """
    import matplotlib.pyplot as plt
    target = data.data[~data.mask]
    if view == 'log':
        vmin = min(target[target>0].min(), theory[theory>0].min())
        vmax = max(target.max(), theory.max())
    else:
        vmin = min(target.min(), theory.min())
        vmax = max(target.max(), theory.max())
    #print vmin, vmax
    plt.subplot(131)
    _plot_2d_signal(data, target, view=view, vmin=vmin, vmax=vmax)
    plt.title('data')
    plt.colorbar()
    plt.subplot(132)
    _plot_2d_signal(data, theory, view=view, vmin=vmin, vmax=vmax)
    plt.title('theory')
    plt.colorbar()
    plt.subplot(133)
    _plot_2d_signal(data, resid, view='linear')
    plt.title('residuals')
    plt.colorbar()

def _plot_2d_signal(data, signal, vmin=None, vmax=None, view='log'):
    """
    Plot the target value for the data.  This could be the data itself,
    the theory calculation, or the residuals.

    *scale* can be 'log' for log scale data, or 'linear'.
    """
    import matplotlib.pyplot as plt
    from numpy.ma import masked_array

    image = np.zeros_like(data.qx_data)
    image[~data.mask] = signal
    valid = np.isfinite(image)
    if view == 'log':
        valid[valid] = (image[valid] > 0)
        image[valid] = np.log10(image[valid])
    elif view == 'q4':
        image[valid] *= (data.qx_data[valid]**2+data.qy_data[valid]**2)**2
    image[~valid | data.mask] = 0
    #plottable = Iq
    plottable = masked_array(image, ~valid | data.mask)
    xmin, xmax = min(data.qx_data), max(data.qx_data)
    ymin, ymax = min(data.qy_data), max(data.qy_data)
    # TODO: fix vmin, vmax so it is shared for theory/resid
    vmin = vmax = None
    try:
        if vmin is None: vmin = image[valid & ~data.mask].min()
        if vmax is None: vmax = image[valid & ~data.mask].max()
    except:
        vmin, vmax = 0, 1
    #print vmin,vmax
    plt.imshow(plottable.reshape(128, 128),
               interpolation='nearest', aspect=1, origin='upper',
               extent=[xmin, xmax, ymin, ymax], vmin=vmin, vmax=vmax)


class Model(object):
    def __init__(self, kernel, **kw):
        from bumps.names import Parameter

        self.kernel = kernel
        partype = kernel.info['partype']

        pars = []
        for p in kernel.info['parameters']:
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

    def parameters(self):
        """
        Return a dictionary of parameters
        """
        return dict((k, getattr(self, k)) for k in self._parameter_names)

class Experiment(object):
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
    def __init__(self, data, model, cutoff=1e-5):

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

        # interpret data
        partype = model.kernel.info['partype']
        if self.data_type == 'sesans':
            q = sesans.make_q(data.sample.zacceptance, data.Rmax)
            self.index = slice(None, None)
            self.Iq = data.y
            self.dIq = data.dy
            #self._theory = np.zeros_like(q)
            q_vectors = [q]
        elif self.data_type == 'Iqxy':
            q = np.sqrt(data.qx_data**2 + data.qy_data**2)
            qmin = getattr(data, 'qmin', 1e-16)
            qmax = getattr(data, 'qmax', np.inf)
            accuracy = getattr(data, 'accuracy', 'Low')
            self.index = (~data.mask) & (~np.isnan(data.data)) \
                         & (q >= qmin) & (q <= qmax)
            self.Iq = data.data[self.index]
            self.dIq = data.err_data[self.index]
            self.resolution = Pinhole2D(data=data, index=self.index,
                                        nsigma=3.0, accuracy=accuracy)
            #self._theory = np.zeros_like(self.Iq)
            if not partype['orientation'] and not partype['magnetic']:
                raise ValueError("not 2D without orientation or magnetic parameters")
                #qx,qy = self.resolution.q_calc
                #q_vectors = [np.sqrt(qx**2 + qy**2)]
            else:
                q_vectors = self.resolution.q_calc
        elif self.data_type == 'Iq':
            self.index = (data.x >= data.qmin) & (data.x <= data.qmax) & ~np.isnan(data.y)
            self.Iq = data.y[self.index]
            self.dIq = data.dy[self.index]
            if getattr(data, 'dx', None) is not None:
                q, dq = data.x[self.index], data.dx[self.index]
                if (dq>0).any():
                    self.resolution = Pinhole1D(q, dq)
                else:
                    self.resolution = Perfect1D(q)
            elif (getattr(data, 'dxl', None) is not None and
                  getattr(data, 'dxw', None) is not None):
                q = data.x[self.index]
                width = data.dxh[self.index]  # Note: dx
                self.resolution = Slit1D(data.x[self.index],
                                         width=data.dxh[self.index],
                                         height=data.dxw[self.index])
            else:
                self.resolution = Perfect1D(data.x[self.index])

            #self._theory = np.zeros_like(self.Iq)
            q_vectors = [self.resolution.q_calc]
        else:
            raise ValueError("Unknown data type") # never gets here

        # Remember function inputs so we can delay loading the function and
        # so we can save/restore state
        self._fn_inputs = [v for v in q_vectors]
        self._fn = None

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
        return self.model.parameters()

    def theory(self):
        if 'theory' not in self._cache:
            if self._fn is None:
                q_input = self.model.kernel.make_input(self._fn_inputs)
                self._fn = self.model.kernel(q_input)

            fixed_pars = [getattr(self.model, p).value for p in self._fn.fixed_pars]
            pd_pars = [self._get_weights(p) for p in self._fn.pd_pars]
            #print fixed_pars,pd_pars
            Iq_calc = self._fn(fixed_pars, pd_pars, self.cutoff)
            #self._theory[:] = self._fn.eval(pars, pd_pars)
            if self.data_type == 'sesans':
                result = sesans.hankel(self.data.x, self.data.lam * 1e-9,
                                       self.data.sample.thickness / 10,
                                       self._fn_inputs[0], Iq_calc)
                self._cache['theory'] = result
            else:
                Iq = self.resolution.apply(Iq_calc)
                self._cache['theory'] = Iq
        return self._cache['theory']

    def residuals(self):
        #if np.any(self.err ==0): print "zeros in err"
        return (self.theory() - self.Iq) / self.dIq

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
        data, theory, resid = self.data, self.theory(), self.residuals()
        if self.data_type == 'Iq':
            _plot_result1D(data, theory, resid, view)
        elif self.data_type == 'Iqxy':
            _plot_result2D(data, theory, resid, view)
        elif self.data_type == 'sesans':
            _plot_sesans(data, theory, resid, view)
        else:
            raise ValueError("Unknown data type")

    def simulate_data(self, noise=None):
        theory = self.theory()
        if noise is not None:
            self.dIq = theory*noise*0.01
        dy = self.dIq
        y = theory + np.random.randn(*dy.shape) * dy
        self.Iq = y
        if self.data_type == 'Iq':
            self.data.dy[self.index] = dy
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

        relative = self.model.kernel.info['partype']['pd-rel']
        limits = self.model.kernel.info['limits']
        disperser, value, npts, width, nsigma = [
            getattr(self.model, par + ext)
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
    plot_data(data)
    import matplotlib.pyplot as plt; plt.show()


if __name__ == "__main__":
    demo()
