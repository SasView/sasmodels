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

import numpy as np  # type: ignore

try:
    from typing import Union, Dict, List, Optional
except ImportError:
    pass
else:
    Data = Union["Data1D", "Data2D", "SesansData"]

def load_data(filename):
    # type: (str) -> Data
    """
    Load data using a sasview loader.
    """
    from sas.sascalc.dataloader.loader import Loader  # type: ignore
    loader = Loader()
    data = loader.load(filename)
    if data is None:
        raise IOError("Data %r could not be loaded" % filename)
    return data


def set_beam_stop(data, radius, outer=None):
    # type: (Data, float, Optional[float]) -> None
    """
    Add a beam stop of the given *radius*.  If *outer*, make an annulus.
    """
    from sas.sascalc.dataloader.manipulations import Ringcut
    if hasattr(data, 'qx_data'):
        data.mask = Ringcut(0, radius)(data)
        if outer is not None:
            data.mask += Ringcut(outer, np.inf)(data)
    else:
        data.mask = (data.x < radius)
        if outer is not None:
            data.mask |= (data.x >= outer)


def set_half(data, half):
    # type: (Data, str) -> None
    """
    Select half of the data, either "right" or "left".
    """
    from sas.sascalc.dataloader.manipulations import Boxcut
    if half == 'right':
        data.mask += \
            Boxcut(x_min=-np.inf, x_max=0.0, y_min=-np.inf, y_max=np.inf)(data)
    if half == 'left':
        data.mask += \
            Boxcut(x_min=0.0, x_max=np.inf, y_min=-np.inf, y_max=np.inf)(data)


def set_top(data, cutoff):
    # type: (Data, float) -> None
    """
    Chop the top off the data, above *cutoff*.
    """
    from sas.sascalc.dataloader.manipulations import Boxcut
    data.mask += \
        Boxcut(x_min=-np.inf, x_max=np.inf, y_min=-np.inf, y_max=cutoff)(data)


class Data1D(object):
    """
    1D data object.

    Note that this definition matches the attributes from sasview, with
    some generic 1D data vectors and some SAS specific definitions.  Some
    refactoring to allow consistent naming conventions between 1D, 2D and
    SESANS data would be helpful.

    **Attributes**

    *x*, *dx*: $q$ vector and gaussian resolution

    *y*, *dy*: $I(q)$ vector and measurement uncertainty

    *mask*: values to include in plotting/analysis

    *dxl*: slit widths for slit smeared data, with *dx* ignored

    *qmin*, *qmax*: range of $q$ values in *x*

    *filename*: label for the data line

    *_xaxis*, *_xunit*: label and units for the *x* axis

    *_yaxis*, *_yunit*: label and units for the *y* axis
    """
    def __init__(self, x=None, y=None, dx=None, dy=None):
        # type: (Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]) -> None
        self.x, self.y, self.dx, self.dy = x, y, dx, dy
        self.dxl = None
        self.filename = None
        self.qmin = x.min() if x is not None else np.NaN
        self.qmax = x.max() if x is not None else np.NaN
        # TODO: why is 1D mask False and 2D mask True?
        self.mask = (np.isnan(y) if y is not None
                     else np.zeros_like(x, 'b') if x is not None
                     else None)
        self._xaxis, self._xunit = "x", ""
        self._yaxis, self._yunit = "y", ""

    def xaxis(self, label, unit):
        # type: (str, str) -> None
        """
        set the x axis label and unit
        """
        self._xaxis = label
        self._xunit = unit

    def yaxis(self, label, unit):
        # type: (str, str) -> None
        """
        set the y axis label and unit
        """
        self._yaxis = label
        self._yunit = unit

class SesansData(Data1D):
    """
    SESANS data object.

    This is just :class:`Data1D` with a wavelength parameter.

    *x* is spin echo length and *y* is polarization (P/P0).
    """
    def __init__(self, **kw):
        Data1D.__init__(self, **kw)
        self.lam = None # type: Optional[np.ndarray]

class Data2D(object):
    """
    2D data object.

    Note that this definition matches the attributes from sasview. Some
    refactoring to allow consistent naming conventions between 1D, 2D and
    SESANS data would be helpful.

    **Attributes**

    *qx_data*, *dqx_data*: $q_x$ matrix and gaussian resolution

    *qy_data*, *dqy_data*: $q_y$ matrix and gaussian resolution

    *data*, *err_data*: $I(q)$ matrix and measurement uncertainty

    *mask*: values to exclude from plotting/analysis

    *qmin*, *qmax*: range of $q$ values in *x*

    *filename*: label for the data line

    *_xaxis*, *_xunit*: label and units for the *x* axis

    *_yaxis*, *_yunit*: label and units for the *y* axis

    *_zaxis*, *_zunit*: label and units for the *y* axis

    *Q_unit*, *I_unit*: units for Q and intensity

    *x_bins*, *y_bins*: grid steps in *x* and *y* directions
    """
    def __init__(self, x=None, y=None, z=None, dx=None, dy=None, dz=None):
        # type: (Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]) -> None
        self.qx_data, self.dqx_data = x, dx
        self.qy_data, self.dqy_data = y, dy
        self.data, self.err_data = z, dz
        self.mask = (np.isnan(z) if z is not None
                     else np.zeros_like(x, dtype='bool') if x is not None
                     else None)
        self.q_data = np.sqrt(x**2 + y**2)
        self.qmin = 1e-16
        self.qmax = np.inf
        self.detector = []
        self.source = Source()
        self.Q_unit = "1/A"
        self.I_unit = "1/cm"
        self.xaxis("Q_x", "1/A")
        self.yaxis("Q_y", "1/A")
        self.zaxis("Intensity", "1/cm")
        self._xaxis, self._xunit = "x", ""
        self._yaxis, self._yunit = "y", ""
        self._zaxis, self._zunit = "z", ""
        self.x_bins, self.y_bins = None, None
        self.filename = None

    def xaxis(self, label, unit):
        # type: (str, str) -> None
        """
        set the x axis label and unit
        """
        self._xaxis = label
        self._xunit = unit

    def yaxis(self, label, unit):
        # type: (str, str) -> None
        """
        set the y axis label and unit
        """
        self._yaxis = label
        self._yunit = unit

    def zaxis(self, label, unit):
        # type: (str, str) -> None
        """
        set the y axis label and unit
        """
        self._zaxis = label
        self._zunit = unit


class Vector(object):
    """
    3-space vector of *x*, *y*, *z*
    """
    def __init__(self, x=None, y=None, z=None):
        # type: (float, float, Optional[float]) -> None
        self.x, self.y, self.z = x, y, z

class Detector(object):
    """
    Detector attributes.
    """
    def __init__(self, pixel_size=(None, None), distance=None):
        # type: (Tuple[float, float], float) -> None
        self.pixel_size = Vector(*pixel_size)
        self.distance = distance

class Source(object):
    """
    Beam attributes.
    """
    def __init__(self):
        # type: () -> None
        self.wavelength = np.NaN
        self.wavelength_unit = "A"


def empty_data1D(q, resolution=0.0):
    # type: (np.ndarray, float) -> Data1D
    """
    Create empty 1D data using the given *q* as the x value.

    *resolution* dq/q defaults to 5%.
    """

    #Iq = 100 * np.ones_like(q)
    #dIq = np.sqrt(Iq)
    Iq, dIq = None, None
    q = np.asarray(q)
    data = Data1D(q, Iq, dx=resolution * q, dy=dIq)
    data.filename = "fake data"
    return data


def empty_data2D(qx, qy=None, resolution=0.0):
    # type: (np.ndarray, Optional[np.ndarray], float) -> Data2D
    """
    Create empty 2D data using the given mesh.

    If *qy* is missing, create a square mesh with *qy=qx*.

    *resolution* dq/q defaults to 5%.
    """
    if qy is None:
        qy = qx
    qx, qy = np.asarray(qx), np.asarray(qy)
    # 5% dQ/Q resolution
    Qx, Qy = np.meshgrid(qx, qy)
    Qx, Qy = Qx.flatten(), Qy.flatten()
    Iq = 100 * np.ones_like(Qx)  # type: np.ndarray
    dIq = np.sqrt(Iq)
    if resolution != 0:
        # https://www.ncnr.nist.gov/staff/hammouda/distance_learning/chapter_15.pdf
        # Should have an additional constant which depends on distances and
        # radii of the aperture, pixel dimensions and wavelength spread
        # Instead, assume radial dQ/Q is constant, and perpendicular matches
        # radial (which instead it should be inverse).
        Q = np.sqrt(Qx**2 + Qy**2)
        dqx = resolution * Q
        dqy = resolution * Q
    else:
        dqx = dqy = None

    data = Data2D(x=Qx, y=Qy, z=Iq, dx=dqx, dy=dqy, dz=dIq)
    data.x_bins = qx
    data.y_bins = qy
    data.filename = "fake data"

    # pixel_size in mm, distance in m
    detector = Detector(pixel_size=(5, 5), distance=4)
    data.detector.append(detector)
    data.source.wavelength = 5 # angstroms
    data.source.wavelength_unit = "A"
    return data


def plot_data(data, view='log', limits=None):
    # type: (Data, str, Optional[Tuple[float, float]]) -> None
    """
    Plot data loaded by the sasview loader.

    *data* is a sasview data object, either 1D, 2D or SESANS.

    *view* is log or linear.

    *limits* sets the intensity limits on the plot; if None then the limits
    are inferred from the data.
    """
    # Note: kind of weird using the plot result functions to plot just the
    # data, but they already handle the masking and graph markup already, so
    # do not repeat.
    if hasattr(data, 'lam'):
        _plot_result_sesans(data, None, None, use_data=True, limits=limits)
    elif hasattr(data, 'qx_data'):
        _plot_result2D(data, None, None, view, use_data=True, limits=limits)
    else:
        _plot_result1D(data, None, None, view, use_data=True, limits=limits)


def plot_theory(data, theory, resid=None, view='log',
                use_data=True, limits=None, Iq_calc=None):
    # type: (Data, Optional[np.ndarray], Optional[np.ndarray], str, bool, Optional[Tuple[float,float]], Optional[np.ndarray]) -> None
    """
    Plot theory calculation.

    *data* is needed to define the graph properties such as labels and
    units, and to define the data mask.

    *theory* is a matrix of the same shape as the data.

    *view* is log or linear

    *use_data* is True if the data should be plotted as well as the theory.

    *limits* sets the intensity limits on the plot; if None then the limits
    are inferred from the data.

    *Iq_calc* is the raw theory values without resolution smearing
    """
    if hasattr(data, 'lam'):
        _plot_result_sesans(data, theory, resid, use_data=True, limits=limits)
    elif hasattr(data, 'qx_data'):
        _plot_result2D(data, theory, resid, view, use_data, limits=limits)
    else:
        _plot_result1D(data, theory, resid, view, use_data,
                       limits=limits, Iq_calc=Iq_calc)


def protect(func):
    # type: (Callable) -> Callable
    """
    Decorator to wrap calls in an exception trapper which prints the
    exception and continues.  Keyboard interrupts are ignored.
    """
    def wrapper(*args, **kw):
        """
        Trap and print errors from function.
        """
        try:
            return func(*args, **kw)
        except Exception:
            traceback.print_exc()

    return wrapper


@protect
def _plot_result1D(data, theory, resid, view, use_data,
                   limits=None, Iq_calc=None):
    # type: (Data1D, Optional[np.ndarray], Optional[np.ndarray], str, bool, Optional[Tuple[float, float]], Optional[np.ndarray]) -> None
    """
    Plot the data and residuals for 1D data.
    """
    import matplotlib.pyplot as plt  # type: ignore
    from numpy.ma import masked_array, masked  # type: ignore

    use_data = use_data and data.y is not None
    use_theory = theory is not None
    use_resid = resid is not None
    use_calc = use_theory and Iq_calc is not None
    num_plots = (use_data or use_theory) + use_calc + use_resid
    non_positive_x = (data.x <= 0.0).any()

    scale = data.x**4 if view == 'q4' else 1.0

    if use_data or use_theory:
        if num_plots > 1:
            plt.subplot(1, num_plots, 1)

        #print(vmin, vmax)
        all_positive = True
        some_present = False
        if use_data:
            mdata = masked_array(data.y, data.mask.copy())
            mdata[~np.isfinite(mdata)] = masked
            if view is 'log':
                mdata[mdata <= 0] = masked
            plt.errorbar(data.x, scale*mdata, yerr=data.dy, fmt='.')
            all_positive = all_positive and (mdata > 0).all()
            some_present = some_present or (mdata.count() > 0)


        if use_theory:
            # Note: masks merge, so any masked theory points will stay masked,
            # and the data mask will be added to it.
            mtheory = masked_array(theory, data.mask.copy())
            mtheory[~np.isfinite(mtheory)] = masked
            if view is 'log':
                mtheory[mtheory <= 0] = masked
            plt.plot(data.x, scale*mtheory, '-', hold=True)
            all_positive = all_positive and (mtheory > 0).all()
            some_present = some_present or (mtheory.count() > 0)

        if limits is not None:
            plt.ylim(*limits)

        plt.xscale('linear' if not some_present or non_positive_x  else view)
        plt.yscale('linear'
                   if view == 'q4' or not some_present or not all_positive
                   else view)
        plt.xlabel("$q$/A$^{-1}$")
        plt.ylabel('$I(q)$')

    if use_calc:
        # Only have use_calc if have use_theory
        plt.subplot(1, num_plots, 2)
        qx, qy, Iqxy = Iq_calc
        plt.pcolormesh(qx, qy[qy > 0], np.log10(Iqxy[qy > 0, :]))
        plt.xlabel("$q_x$/A$^{-1}$")
        plt.xlabel("$q_y$/A$^{-1}$")
        plt.xscale('log')
        plt.yscale('log')
        #plt.axis('equal')

    if use_resid:
        mresid = masked_array(resid, data.mask.copy())
        mresid[~np.isfinite(mresid)] = masked
        some_present = (mresid.count() > 0)

        if num_plots > 1:
            plt.subplot(1, num_plots, use_calc + 2)
        plt.plot(data.x, mresid, '-')
        plt.xlabel("$q$/A$^{-1}$")
        plt.ylabel('residuals')
        plt.xscale('linear' if not some_present or non_positive_x else view)


@protect
def _plot_result_sesans(data, theory, resid, use_data, limits=None):
    # type: (SesansData, Optional[np.ndarray], Optional[np.ndarray], bool, Optional[Tuple[float, float]]) -> None
    """
    Plot SESANS results.
    """
    import matplotlib.pyplot as plt  # type: ignore
    use_data = use_data and data.y is not None
    use_theory = theory is not None
    use_resid = resid is not None
    num_plots = (use_data or use_theory) + use_resid

    if use_data or use_theory:
        is_tof = (data.lam != data.lam[0]).any()
        if num_plots > 1:
            plt.subplot(1, num_plots, 1)
        if use_data:
            if is_tof:
                plt.errorbar(data.x, np.log(data.y)/(data.lam*data.lam),
                             yerr=data.dy/data.y/(data.lam*data.lam))
            else:
                plt.errorbar(data.x, data.y, yerr=data.dy)
        if theory is not None:
            if is_tof:
                plt.plot(data.x, np.log(theory)/(data.lam*data.lam), '-', hold=True)
            else:
                plt.plot(data.x, theory, '-', hold=True)
        if limits is not None:
            plt.ylim(*limits)

        plt.xlabel('spin echo length ({})'.format(data._xunit))
        if is_tof:
            plt.ylabel(r'(Log (P/P$_0$))/$\lambda^2$')
        else:
            plt.ylabel('polarization (P/P0)')


    if resid is not None:
        if num_plots > 1:
            plt.subplot(1, num_plots, (use_data or use_theory) + 1)
        plt.plot(data.x, resid, 'x')
        plt.xlabel('spin echo length ({})'.format(data._xunit))
        plt.ylabel('residuals (P/P0)')


@protect
def _plot_result2D(data, theory, resid, view, use_data, limits=None):
    # type: (Data2D, Optional[np.ndarray], Optional[np.ndarray], str, bool, Optional[Tuple[float,float]]) -> None
    """
    Plot the data and residuals for 2D data.
    """
    import matplotlib.pyplot as plt  # type: ignore
    use_data = use_data and data.data is not None
    use_theory = theory is not None
    use_resid = resid is not None
    num_plots = use_data + use_theory + use_resid

    # Put theory and data on a common colormap scale
    vmin, vmax = np.inf, -np.inf
    target = None # type: Optional[np.ndarray]
    if use_data:
        target = data.data[~data.mask]
        datamin = target[target > 0].min() if view == 'log' else target.min()
        datamax = target.max()
        vmin = min(vmin, datamin)
        vmax = max(vmax, datamax)
    if use_theory:
        theorymin = theory[theory > 0].min() if view == 'log' else theory.min()
        theorymax = theory.max()
        vmin = min(vmin, theorymin)
        vmax = max(vmax, theorymax)

    # Override data limits from the caller
    if limits is not None:
        vmin, vmax = limits

    # Plot data
    if use_data:
        if num_plots > 1:
            plt.subplot(1, num_plots, 1)
        _plot_2d_signal(data, target, view=view, vmin=vmin, vmax=vmax)
        plt.title('data')
        h = plt.colorbar()
        h.set_label('$I(q)$')

    # plot theory
    if use_theory:
        if num_plots > 1:
            plt.subplot(1, num_plots, use_data+1)
        _plot_2d_signal(data, theory, view=view, vmin=vmin, vmax=vmax)
        plt.title('theory')
        h = plt.colorbar()
        h.set_label(r'$\log_{10}I(q)$' if view == 'log'
                    else r'$q^4 I(q)$' if view == 'q4'
                    else '$I(q)$')

    # plot resid
    if use_resid:
        if num_plots > 1:
            plt.subplot(1, num_plots, use_data+use_theory+1)
        _plot_2d_signal(data, resid, view='linear')
        plt.title('residuals')
        h = plt.colorbar()
        h.set_label(r'$\Delta I(q)$')


@protect
def _plot_2d_signal(data, signal, vmin=None, vmax=None, view='log'):
    # type: (Data2D, np.ndarray, Optional[float], Optional[float], str) -> Tuple[float, float]
    """
    Plot the target value for the data.  This could be the data itself,
    the theory calculation, or the residuals.

    *scale* can be 'log' for log scale data, or 'linear'.
    """
    import matplotlib.pyplot as plt  # type: ignore
    from numpy.ma import masked_array  # type: ignore

    image = np.zeros_like(data.qx_data)
    image[~data.mask] = signal
    valid = np.isfinite(image)
    if view == 'log':
        valid[valid] = (image[valid] > 0)
        if vmin is None: vmin = image[valid & ~data.mask].min()
        if vmax is None: vmax = image[valid & ~data.mask].max()
        image[valid] = np.log10(image[valid])
    elif view == 'q4':
        image[valid] *= (data.qx_data[valid]**2+data.qy_data[valid]**2)**2
        if vmin is None: vmin = image[valid & ~data.mask].min()
        if vmax is None: vmax = image[valid & ~data.mask].max()
    else:
        if vmin is None: vmin = image[valid & ~data.mask].min()
        if vmax is None: vmax = image[valid & ~data.mask].max()

    image[~valid | data.mask] = 0
    #plottable = Iq
    plottable = masked_array(image, ~valid | data.mask)
    # Divide range by 10 to convert from angstroms to nanometers
    xmin, xmax = min(data.qx_data), max(data.qx_data)
    ymin, ymax = min(data.qy_data), max(data.qy_data)
    if view == 'log':
        vmin, vmax = np.log10(vmin), np.log10(vmax)
    plt.imshow(plottable.reshape(len(data.x_bins), len(data.y_bins)),
               interpolation='nearest', aspect=1, origin='lower',
               extent=[xmin, xmax, ymin, ymax], vmin=vmin, vmax=vmax)
    plt.xlabel("$q_x$/A$^{-1}$")
    plt.ylabel("$q_y$/A$^{-1}$")
    return vmin, vmax

def demo():
    # type: () -> None
    """
    Load and plot a SAS dataset.
    """
    data = load_data('DEC07086.DAT')
    set_beam_stop(data, 0.004)
    plot_data(data)
    import matplotlib.pyplot as plt  # type: ignore
    plt.show()


if __name__ == "__main__":
    demo()
