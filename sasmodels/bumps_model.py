"""
Sasmodels core.
"""
import sys, os
import datetime

import numpy as np

try:
    from .gpu import load_model as _loader
except ImportError,exc:
    import warnings
    warnings.warn(str(exc))
    warnings.warn("OpenCL not available --- using ctypes instead")
    from .dll import load_model as _loader

def load_model(modelname, dtype='single'):
    """
    Load model by name.
    """
    sasmodels = __import__('sasmodels.models.'+modelname)
    module = getattr(sasmodels.models, modelname, None)
    model = _loader(module, dtype=dtype)
    return model


def tic():
    """
    Timer function.

    Use "toc=tic()" to start the clock and "toc()" to measure
    a time interval.
    """
    then = datetime.datetime.now()
    return lambda: (datetime.datetime.now()-then).total_seconds()


def load_data(filename):
    """
    Load data using a sasview loader.
    """
    from sans.dataloader.loader import Loader
    loader = Loader()
    data = loader.load(filename)
    if data is None:
        raise IOError("Data %r could not be loaded"%filename)
    return data


def empty_data1D(q):
    """
    Create empty 1D data using the given *q* as the x value.

    Resolutions dq/q is 5%.
    """

    from sans.dataloader.data_info import Data1D

    Iq = 100*np.ones_like(q)
    dIq = np.sqrt(Iq)
    data = Data1D(q, Iq, dx=0.05*q, dy=dIq)
    data.filename = "fake data"
    data.qmin, data.qmax = q.min(), q.max()
    return data


def empty_data2D(qx, qy=None):
    """
    Create empty 2D data using the given mesh.

    If *qy* is missing, create a square mesh with *qy=qx*.

    Resolution dq/q is 5%.
    """
    from sans.dataloader.data_info import Data2D, Detector

    if qy is None:
        qy = qx
    Qx,Qy = np.meshgrid(qx,qy)
    Qx,Qy = Qx.flatten(), Qy.flatten()
    Iq = 100*np.ones_like(Qx)
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
    data.dqx_data = 0.05*Qx
    data.dqy_data = 0.05*Qy

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
    data.q_data = np.sqrt(Qx**2 + Qy**2)
    data.xaxis("Q_x", "A^{-1}")
    data.yaxis("Q_y", "A^{-1}")
    data.zaxis("Intensity", r"\text{cm}^{-1}")
    return data


def set_beam_stop(data, radius, outer=None):
    """
    Add a beam stop of the given *radius*.  If *outer*, make an annulus.
    """
    from sans.dataloader.manipulations import Ringcut
    if hasattr(data, 'qx_data'):
        data.mask = Ringcut(0, radius)(data)
        if outer is not None:
            data.mask += Ringcut(outer,np.inf)(data)
    else:
        data.mask = (data.x>=radius)
        if outer is not None:
            data.mask &= (data.x<outer)


def set_half(data, half):
    """
    Select half of the data, either "right" or "left".
    """
    from sans.dataloader.manipulations import Boxcut
    if half == 'right':
        data.mask += Boxcut(x_min=-np.inf, x_max=0.0, y_min=-np.inf, y_max=np.inf)(data)
    if half == 'left':
        data.mask += Boxcut(x_min=0.0, x_max=np.inf, y_min=-np.inf, y_max=np.inf)(data)


def set_top(data, max):
    """
    Chop the top off the data, above *max*.
    """
    from sans.dataloader.manipulations import Boxcut
    data.mask += Boxcut(x_min=-np.inf, x_max=np.inf, y_min=-np.inf, y_max=max)(data)


def plot_data(data, iq, vmin=None, vmax=None, scale='log'):
    """
    Plot the target value for the data.  This could be the data itself,
    the theory calculation, or the residuals.

    *scale* can be 'log' for log scale data, or 'linear'.
    """
    from numpy.ma import masked_array, masked
    import matplotlib.pyplot as plt
    if hasattr(data, 'qx_data'):
        iq = iq[:]
        valid = np.isfinite(iq)
        if scale == 'log':
            valid[valid] = (iq[valid] > 0)
            iq[valid] = np.log10(iq[valid])
        iq[~valid|data.mask] = 0
        #plottable = iq
        plottable = masked_array(iq, ~valid|data.mask)
        xmin, xmax = min(data.qx_data), max(data.qx_data)
        ymin, ymax = min(data.qy_data), max(data.qy_data)
        plt.imshow(plottable.reshape(128,128),
                   interpolation='nearest', aspect=1, origin='upper',
                   extent=[xmin, xmax, ymin, ymax], vmin=vmin, vmax=vmax)
    else: # 1D data
        if scale == 'linear':
            idx = np.isfinite(iq)
            plt.plot(data.x[idx], iq[idx])
        else:
            idx = np.isfinite(iq)
            idx[idx] = (iq[idx]>0)
            plt.loglog(data.x[idx], iq[idx])


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
    mresid = masked_array((theory-data.y)/data.dy, mdata.mask)

    plt.subplot(121)
    plt.errorbar(data.x, mdata, yerr=data.dy)
    plt.plot(data.x, mtheory, '-', hold=True)
    plt.yscale(view)
    plt.subplot(122)
    plt.plot(data.x, mresid, 'x')
    #plt.axhline(1, color='black', ls='--',lw=1, hold=True)
    #plt.axhline(0, color='black', lw=1, hold=True)
    #plt.axhline(-1, color='black', ls='--',lw=1, hold=True)


def _plot_result2D(data, theory, view):
    """
    Plot the data and residuals for 2D data.
    """
    import matplotlib.pyplot as plt
    resid = (theory-data.data)/data.err_data
    plt.subplot(131)
    plot_data(data, data.data, scale=view)
    plt.colorbar()
    plt.subplot(132)
    plot_data(data, theory, scale=view)
    plt.colorbar()
    plt.subplot(133)
    plot_data(data, resid, scale='linear')
    plt.colorbar()

def plot_result(data, theory, view='log'):
    """
    Plot the data and residuals.
    """
    if hasattr(data, 'qx_data'):
        _plot_result2D(data, theory, view)
    else:
        _plot_result1D(data, theory, view)


class BumpsModel(object):
    """
    Return a bumps wrapper for a SAS model.

    *data* is the data to be fitted.

    *model* is the SAS model, e.g., from :func:`gen.opencl_model`.

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

        partype = model.info['partype']

        # interpret data
        if hasattr(data, 'qx_data'):
            self.index = (data.mask==0) & (~np.isnan(data.data))
            self.iq = data.data[self.index]
            self.diq = data.err_data[self.index]
            self._theory = np.zeros_like(data.data)
            if not partype['orientation'] and not partype['magnetic']:
                q_vectors = [np.sqrt(data.qx_data**2+data.qy_data**2)]
            else:
                q_vectors = [data.qx_data, data.qy_data]
        else:
            self.index = (data.x>=data.qmin) & (data.x<=data.qmax) & ~np.isnan(data.y)
            self.iq = data.y[self.index]
            self.diq = data.dy[self.index]
            self._theory = np.zeros_like(data.y)
            q_vectors = [data.x]

        # Remember function inputs so we can delay loading the function and
        # so we can save/restore state
        self._fn_inputs = [v[self.index] for v in q_vectors]
        self._fn = None

        # define bumps parameters
        pars = []
        for p in model.info['parameters']:
            name, default, limits, ptype = p[0], p[2], p[3], p[4]
            value = kw.pop(name, default)
            setattr(self, name, Parameter.default(value, name=name, limits=limits))
            pars.append(name)
        for name in partype['pd-2d']:
            for xpart,xdefault,xlimits in [
                    ('_pd', 0, limits),
                    ('_pd_n', 35, (0,1000)),
                    ('_pd_nsigma', 3, (0, 10)),
                    ('_pd_type', 'gaussian', None),
                ]:
                xname = name+xpart
                xvalue = kw.pop(xname, xdefault)
                if xlimits is not None:
                    xvalue = Parameter.default(xvalue, name=xname, limits=xlimits)
                    pars.append(xname)
                setattr(self, xname, xvalue)
        self._parameter_names = pars
        if kw:
            raise TypeError("unexpected parameters: %s"%(", ".join(sorted(kw.keys()))))
        self.update()

    def update(self):
        self._cache = {}

    def numpoints(self):
        return len(self.iq)

    def parameters(self):
        return dict((k,getattr(self,k)) for k in self._parameter_names)

    def theory(self):
        if 'theory' not in self._cache:
            if self._fn is None:
                input = self.model.make_input(self._fn_inputs)
                self._fn = self.model(input)

            pars = [getattr(self,p).value for p in self._fn.fixed_pars]
            pd_pars = [self._get_weights(p) for p in self._fn.pd_pars]
            #print pars
            self._theory[self.index] = self._fn(pars, pd_pars, self.cutoff)
            #self._theory[:] = self._fn.eval(pars, pd_pars)
            self._cache['theory'] = self._theory
        return self._cache['theory']

    def residuals(self):
        #if np.any(self.err ==0): print "zeros in err"
        return (self.theory()[self.index]-self.iq)/self.diq

    def nllf(self):
        R = self.residuals()
        #if np.any(np.isnan(R)): print "NaN in residuals"
        return 0.5*np.sum(R**2)

    def __call__(self):
        return 2*self.nllf()/self.dof

    def plot(self, view='log'):
        plot_result(self.data, self.theory(), view=view)

    def save(self, basename):
        pass

    def _get_weights(self, par):
        from . import weights

        relative = self.model.info['partype']['pd-rel']
        limits = self.model.info['limits']
        disperser,value,npts,width,nsigma = [getattr(self, par+ext)
                for ext in ('_pd_type','','_pd_n','_pd','_pd_nsigma')]
        v,w = weights.get_weights(
            disperser, int(npts.value), width.value, nsigma.value,
            value.value, limits[par], par in relative)
        return v,w/w.max()

    def __getstate__(self):
        # Can't pickle gpu functions, so instead make them lazy
        state = self.__dict__.copy()
        state['_fn'] = None
        return state

    def __setstate__(self, state):
        self.__dict__ = state


def demo():
    data = load_data('DEC07086.DAT')
    set_beam_stop(data, 0.004)
    plot_data(data)
    import matplotlib.pyplot as plt; plt.show()


if __name__ == "__main__":
    demo()
