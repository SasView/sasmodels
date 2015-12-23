"""
SAS data representations.

Plotting functions for data sets:

    :func:`plot_data` plots the data file.

    :func:`plot_theory` plots a calculated result from the model.

Wrappers for the sasview data loader and data manipulations:

    :func:`load_data` loads a sasview data file.

    :func:`set_beam_stop` masks the beam stop from the data.

    :func:`set_half` selects the right or left half of the data, which can
    be useful for shear measurements which have not been properly corrected
    for path length and reflections.

    :func:`set_top` cuts the top part off the data.


Empty data sets for evaluating models without data:

    :func:`empty_data1D` creates an empty dataset, which is useful for plotting
    a theory function before the data is measured.

    :func:`empty_data2D` creates an empty 2D dataset.

Note that the empty datasets use a minimal representation of the SasView
objects so that models can be run without SasView on the path.  You could
also use these for your own data loader.

"""
import traceback

import numpy as np

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
        data.mask = (data.x < radius)
        if outer is not None:
            data.mask |= (data.x >= outer)


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


class Data1D(object):
    def __init__(self, x=None, y=None, dx=None, dy=None):
        self.x, self.y, self.dx, self.dy = x, y, dx, dy
        self.dxl = None

    def xaxis(self, label, unit):
        """
        set the x axis label and unit
        """
        self._xaxis = label
        self._xunit = unit

    def yaxis(self, label, unit):
        """
        set the y axis label and unit
        """
        self._yaxis = label
        self._yunit = unit



class Data2D(object):
    def __init__(self):
        self.detector = []
        self.source = Source()

    def xaxis(self, label, unit):
        """
        set the x axis label and unit
        """
        self._xaxis = label
        self._xunit = unit

    def yaxis(self, label, unit):
        """
        set the y axis label and unit
        """
        self._yaxis = label
        self._yunit = unit

    def zaxis(self, label, unit):
        """
        set the y axis label and unit
        """
        self._zaxis = label
        self._zunit = unit


class Vector(object):
    def __init__(self, x=None, y=None, z=None):
        self.x, self.y, self.z = x, y, z

class Detector(object):
    def __init__(self):
        self.pixel_size = Vector()

class Source(object):
    pass


def empty_data1D(q, resolution=0.05):
    """
    Create empty 1D data using the given *q* as the x value.

    *resolution* dq/q defaults to 5%.
    """

    #Iq = 100 * np.ones_like(q)
    #dIq = np.sqrt(Iq)
    Iq, dIq = None, None
    data = Data1D(q, Iq, dx=resolution * q, dy=dIq)
    data.filename = "fake data"
    data.qmin, data.qmax = q.min(), q.max()
    data.mask = np.zeros(len(q), dtype='bool')
    return data


def empty_data2D(qx, qy=None, resolution=0.05):
    """
    Create empty 2D data using the given mesh.

    If *qy* is missing, create a square mesh with *qy=qx*.

    *resolution* dq/q defaults to 5%.
    """
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
    else:
        data.dqx_data = data.dqy_data = None

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


def plot_data(data, view='log'):
    """
    Plot data loaded by the sasview loader.
    """
    # Note: kind of weird using the plot result functions to plot just the
    # data, but they already handle the masking and graph markup already, so
    # do not repeat.
    if hasattr(data, 'lam'):
        _plot_result_sesans(data, None, None, plot_data=True)
    elif hasattr(data, 'qx_data'):
        _plot_result2D(data, None, None, view, plot_data=True)
    else:
        _plot_result1D(data, None, None, view, plot_data=True)


def plot_theory(data, theory, resid=None, view='log', plot_data=True):
    if hasattr(data, 'lam'):
        _plot_result_sesans(data, theory, resid, plot_data=True)
    elif hasattr(data, 'qx_data'):
        _plot_result2D(data, theory, resid, view, plot_data)
    else:
        _plot_result1D(data, theory, resid, view, plot_data)


def protect(fn):
    def wrapper(*args, **kw):
        try:
            return fn(*args, **kw)
        except:
            traceback.print_exc()
            pass

    return wrapper


@protect
def _plot_result1D(data, theory, resid, view, plot_data):
    """
    Plot the data and residuals for 1D data.
    """
    import matplotlib.pyplot as plt
    from numpy.ma import masked_array, masked

    plot_theory = theory is not None
    plot_resid = resid is not None

    if data.y is None:
        plot_data = False
    scale = data.x**4 if view == 'q4' else 1.0

    if plot_data or plot_theory:
        if plot_resid:
            plt.subplot(121)

        #print(vmin, vmax)
        all_positive = True
        some_present = False
        if plot_data:
            mdata = masked_array(data.y, data.mask.copy())
            mdata[~np.isfinite(mdata)] = masked
            if view is 'log':
                mdata[mdata <= 0] = masked
            plt.errorbar(data.x/10, scale*mdata, yerr=data.dy, fmt='.')
            all_positive = all_positive and (mdata>0).all()
            some_present = some_present or (mdata.count() > 0)


        if plot_theory:
            mtheory = masked_array(theory, data.mask.copy())
            mtheory[~np.isfinite(mtheory)] = masked
            if view is 'log':
                mtheory[mtheory<=0] = masked
            plt.plot(data.x/10, scale*mtheory, '-', hold=True)
            all_positive = all_positive and (mtheory>0).all()
            some_present = some_present or (mtheory.count() > 0)

        plt.xscale('linear' if not some_present else view)
        plt.yscale('linear'
                   if view == 'q4' or not some_present or not all_positive
                   else view)
        plt.xlabel("$q$/nm$^{-1}$")
        plt.ylabel('$I(q)$')

    if plot_resid:
        if plot_data or plot_theory:
            plt.subplot(122)

        mresid = masked_array(resid, data.mask.copy())
        mresid[~np.isfinite(mresid)] = masked
        some_present = (mresid.count() > 0)
        plt.plot(data.x/10, mresid, '-')
        plt.xlabel("$q$/nm$^{-1}$")
        plt.ylabel('residuals')
        plt.xscale('linear' if not some_present else view)


@protect
def _plot_result_sesans(data, theory, resid, plot_data):
    import matplotlib.pyplot as plt
    if data.y is None:
        plot_data = False
    plot_theory = theory is not None
    plot_resid = resid is not None

    if plot_data or plot_theory:
        if plot_resid:
            plt.subplot(121)
        if plot_data:
            plt.errorbar(data.x, data.y, yerr=data.dy)
        if theory is not None:
            plt.plot(data.x, theory, '-', hold=True)
        plt.xlabel('spin echo length (nm)')
        plt.ylabel('polarization (P/P0)')

    if resid is not None:
        if plot_data or plot_theory:
            plt.subplot(122)

        plt.plot(data.x, resid, 'x')
        plt.xlabel('spin echo length (nm)')
        plt.ylabel('residuals (P/P0)')


@protect
def _plot_result2D(data, theory, resid, view, plot_data):
    """
    Plot the data and residuals for 2D data.
    """
    import matplotlib.pyplot as plt
    if data.data is None:
        plot_data = False
    plot_theory = theory is not None
    plot_resid = resid is not None

    # Put theory and data on a common colormap scale
    vmin, vmax = np.inf, -np.inf
    if plot_data:
        target = data.data[~data.mask]
        datamin = target[target>0].min() if view == 'log' else target.min()
        datamax = target.max()
        vmin = min(vmin, datamin)
        vmax = max(vmax, datamax)
    if plot_theory:
        theorymin = theory[theory>0].min() if view == 'log' else theory.min()
        theorymax = theory.max()
        vmin = min(vmin, theorymin)
        vmax = max(vmax, theorymax)

    if plot_data:
        if plot_theory and plot_resid:
            plt.subplot(131)
        elif plot_theory or plot_resid:
            plt.subplot(121)
        _plot_2d_signal(data, target, view=view, vmin=vmin, vmax=vmax)
        plt.title('data')
        h = plt.colorbar()
        h.set_label('$I(q)$')

    if plot_theory:
        if plot_data and plot_resid:
            plt.subplot(132)
        elif plot_data:
            plt.subplot(122)
        elif plot_resid:
            plt.subplot(121)
        _plot_2d_signal(data, theory, view=view, vmin=vmin, vmax=vmax)
        plt.title('theory')
        h = plt.colorbar()
        h.set_label('$I(q)$')

    #if plot_data or plot_theory:
    #    plt.colorbar()

    if plot_resid:
        if plot_data and plot_theory:
            plt.subplot(133)
        elif plot_data or plot_theory:
            plt.subplot(122)
        _plot_2d_signal(data, resid, view='linear')
        plt.title('residuals')
        h = plt.colorbar()
        h.set_label('$\Delta I(q)$')


@protect
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
    xmin, xmax = min(data.qx_data)/10, max(data.qx_data)/10
    ymin, ymax = min(data.qy_data)/10, max(data.qy_data)/10
    # TODO: fix vmin, vmax so it is shared for theory/resid
    vmin = vmax = None
    try:
        if vmin is None: vmin = image[valid & ~data.mask].min()
        if vmax is None: vmax = image[valid & ~data.mask].max()
    except:
        vmin, vmax = 0, 1
    plt.imshow(plottable.reshape(len(data.xbins), len(data.ybins)),
               interpolation='nearest', aspect=1, origin='upper',
               extent=[xmin, xmax, ymin, ymax], vmin=vmin, vmax=vmax)
    plt.xlabel("$q_x$/nm$^{-1}$")
    plt.ylabel("$q_y$/nm$^{-1}$")


def demo():
    data = load_data('DEC07086.DAT')
    set_beam_stop(data, 0.004)
    plot_data(data)
    import matplotlib.pyplot as plt; plt.show()


if __name__ == "__main__":
    demo()
