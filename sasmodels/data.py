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
from numpy import sqrt, sin, cos, pi

# pylint: disable=unused-import
try:
    from typing import Union, Dict, List, Optional
except ImportError:
    pass
else:
    Data = Union["Data1D", "Data2D", "SesansData"]
# pylint: enable=unused-import

def load_data(filename, index=0):
    # type: (str) -> Data
    """
    Load data using a sasview loader.
    """
    from sas.sascalc.dataloader.loader import Loader  # type: ignore
    loader = Loader()
    # Allow for one part in multipart file
    if '[' in filename:
        filename, indexstr = filename[:-1].split('[')
        index = int(indexstr)
    datasets = loader.load(filename)
    if not datasets:  # None or []
        raise IOError("Data %r could not be loaded" % filename)
    if not isinstance(datasets, list):
        datasets = [datasets]
    for data in datasets:
        if hasattr(data, 'x'):
            data.qmin, data.qmax = data.x.min(), data.x.max()
            data.mask = (np.isnan(data.y) if data.y is not None
                         else np.zeros_like(data.x, dtype='bool'))
        elif hasattr(data, 'qx_data'):
            data.mask = ~data.mask
    return datasets[index] if index != 'all' else datasets


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
    def __init__(self,
                 x=None,  # type: Optional[np.ndarray]
                 y=None,  # type: Optional[np.ndarray]
                 dx=None, # type: Optional[np.ndarray]
                 dy=None  # type: Optional[np.ndarray]
                ):
        # type: (...) -> None
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
    isSesans = True
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
    def __init__(self,
                 x=None,   # type: Optional[np.ndarray]
                 y=None,   # type: Optional[np.ndarray]
                 z=None,   # type: Optional[np.ndarray]
                 dx=None,  # type: Optional[np.ndarray]
                 dy=None,  # type: Optional[np.ndarray]
                 dz=None   # type: Optional[np.ndarray]
                ):
        # type: (...) -> None
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

class Sample(object):
    """
    Sample attributes.
    """
    def __init__(self):
        # type: () -> None
        pass

def empty_data1D(q, resolution=0.0, L=0., dL=0.):
    # type: (np.ndarray, float) -> Data1D
    r"""
    Create empty 1D data using the given *q* as the x value.

    rms *resolution* $\Delta q/q$ defaults to 0%.  If wavelength *L* and rms
    wavelength divergence *dL* are defined, then *resolution* defines
    rms $\Delta \theta/\theta$ for the lowest *q*, with $\theta$ derived from
    $q = 4\pi/\lambda \sin(\theta)$.
    """

    #Iq = 100 * np.ones_like(q)
    #dIq = np.sqrt(Iq)
    Iq, dIq = None, None
    q = np.asarray(q)
    if L != 0 and resolution != 0:
        theta = np.arcsin(q*L/(4*pi))
        dtheta = theta[0]*resolution
        ## Solving Gaussian error propagation from
        ##   Dq^2 = (dq/dL)^2 DL^2 + (dq/dtheta)^2 Dtheta^2
        ## gives
        ##   (Dq/q)^2 = (DL/L)**2 + (Dtheta/tan(theta))**2
        ## Take the square root and multiply by q, giving
        ##   Dq = (4*pi/L) * sqrt((sin(theta)*dL/L)**2 + (cos(theta)*dtheta)**2)
        dq = (4*pi/L) * sqrt((sin(theta)*dL/L)**2 + (cos(theta)*dtheta)**2)
    else:
        dq = resolution * q

    data = Data1D(q, Iq, dx=dq, dy=dIq)
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
    if hasattr(data, 'isSesans') and data.isSesans:
        _plot_result_sesans(data, None, None, use_data=True, limits=limits)
    elif hasattr(data, 'qx_data') and not getattr(data, 'radial', False):
        _plot_result2D(data, None, None, view, use_data=True, limits=limits)
    else:
        _plot_result1D(data, None, None, view, use_data=True, limits=limits)


def plot_theory(data,          # type: Data
                theory,        # type: Optional[np.ndarray]
                resid=None,    # type: Optional[np.ndarray]
                view='log',    # type: str
                use_data=True, # type: bool
                limits=None,   # type: Optional[np.ndarray]
                Iq_calc=None   # type: Optional[np.ndarray]
               ):
    # type: (...) -> None
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
    if hasattr(data, 'isSesans') and data.isSesans:
        _plot_result_sesans(data, theory, resid, use_data=True, limits=limits)
    elif hasattr(data, 'qx_data') and not getattr(data, 'radial', False):
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
def _plot_result1D(data,         # type: Data1D
                   theory,       # type: Optional[np.ndarray]
                   resid,        # type: Optional[np.ndarray]
                   view,         # type: str
                   use_data,     # type: bool
                   limits=None,  # type: Optional[Tuple[float, float]]
                   Iq_calc=None  # type: Optional[np.ndarray]
                  ):
    # type: (...) -> None
    """
    Plot the data and residuals for 1D data.
    """
    import matplotlib.pyplot as plt  # type: ignore
    from numpy.ma import masked_array, masked  # type: ignore

    if getattr(data, 'radial', False):
        data.x = data.q_data
        data.y = data.data

    use_data = use_data and data.y is not None
    use_theory = theory is not None
    use_resid = resid is not None
    use_calc = use_theory and Iq_calc is not None
    num_plots = (use_data or use_theory) + use_calc + use_resid
    non_positive_x = (data.x <= 0.0).any()

    scale = data.x**4 if view == 'q4' else 1.0
    xscale = yscale = 'linear' if view == 'linear' else 'log'

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
            #mtheory = masked_array(theory, data.mask.copy())
            theory_x = data.x[data.mask == 0]
            mtheory = masked_array(theory)
            mtheory[~np.isfinite(mtheory)] = masked
            if view is 'log':
                mtheory[mtheory <= 0] = masked
            plt.plot(theory_x, scale*mtheory, '-')
            all_positive = all_positive and (mtheory > 0).all()
            some_present = some_present or (mtheory.count() > 0)

        if limits is not None:
            plt.ylim(*limits)


        xscale = ('linear' if not some_present or non_positive_x
                  else view if view is not None
                  else 'log')
        yscale = ('linear'
                  if view == 'q4' or not some_present or not all_positive
                  else view if view is not None
                  else 'log')
        plt.xscale(xscale)
        plt.xlabel("$q$/A$^{-1}$")
        plt.yscale(yscale)
        plt.ylabel('$I(q)$')
        title = ("data and model" if use_theory and use_data
                 else "data" if use_data
                 else "model")
        plt.title(title)

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
        theory_x = data.x[data.mask == 0]
        mresid = masked_array(resid)
        mresid[~np.isfinite(mresid)] = masked
        some_present = (mresid.count() > 0)

        if num_plots > 1:
            plt.subplot(1, num_plots, use_calc + 2)
        plt.plot(theory_x, mresid, '.')
        plt.xlabel("$q$/A$^{-1}$")
        plt.ylabel('residuals')
        plt.title('(model - Iq)/dIq')
        plt.xscale(xscale)
        plt.yscale('linear')


@protect
def _plot_result_sesans(data,        # type: SesansData
                        theory,      # type: Optional[np.ndarray]
                        resid,       # type: Optional[np.ndarray]
                        use_data,    # type: bool
                        limits=None  # type: Optional[Tuple[float, float]]
                       ):
    # type: (...) -> None
    """
    Plot SESANS results.
    """
    import matplotlib.pyplot as plt  # type: ignore
    use_data = use_data and data.y is not None
    use_theory = theory is not None
    use_resid = resid is not None
    num_plots = (use_data or use_theory) + use_resid

    if use_data or use_theory:
        is_tof = data.lam is not None and (data.lam != data.lam[0]).any()
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
                plt.plot(data.x, np.log(theory)/(data.lam*data.lam), '-')
            else:
                plt.plot(data.x, theory, '-')
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
def _plot_result2D(data,         # type: Data2D
                   theory,       # type: Optional[np.ndarray]
                   resid,        # type: Optional[np.ndarray]
                   view,         # type: str
                   use_data,     # type: bool
                   limits=None   # type: Optional[Tuple[float, float]]
                  ):
    # type: (...) -> None
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
def _plot_2d_signal(data,       # type: Data2D
                    signal,     # type: np.ndarray
                    vmin=None,  # type: Optional[float]
                    vmax=None,  # type: Optional[float]
                    view='log'  # type: str
                   ):
    # type: (...) -> Tuple[float, float]
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
        if vmin is None:
            vmin = image[valid & ~data.mask].min()
        if vmax is None:
            vmax = image[valid & ~data.mask].max()
        image[valid] = np.log10(image[valid])
    elif view == 'q4':
        image[valid] *= (data.qx_data[valid]**2+data.qy_data[valid]**2)**2
        if vmin is None:
            vmin = image[valid & ~data.mask].min()
        if vmax is None:
            vmax = image[valid & ~data.mask].max()
    else:
        if vmin is None:
            vmin = image[valid & ~data.mask].min()
        if vmax is None:
            vmax = image[valid & ~data.mask].max()

    image[~valid | data.mask] = 0
    #plottable = Iq
    plottable = masked_array(image, ~valid | data.mask)
    # Divide range by 10 to convert from angstroms to nanometers
    xmin, xmax = min(data.qx_data), max(data.qx_data)
    ymin, ymax = min(data.qy_data), max(data.qy_data)
    if view == 'log':
        vmin_scaled, vmax_scaled = np.log10(vmin), np.log10(vmax)
    else:
        vmin_scaled, vmax_scaled = vmin, vmax
    #nx, ny = len(data.x_bins), len(data.y_bins)
    x_bins, y_bins, image = _build_matrix(data, plottable)
    plt.imshow(image,
               interpolation='nearest', aspect=1, origin='lower',
               extent=[xmin, xmax, ymin, ymax],
               vmin=vmin_scaled, vmax=vmax_scaled)
    plt.xlabel("$q_x$/A$^{-1}$")
    plt.ylabel("$q_y$/A$^{-1}$")
    return vmin, vmax


# === The following is modified from sas.sasgui.plottools.PlotPanel
def _build_matrix(self, plottable):
    """
    Build a matrix for 2d plot from a vector
    Returns a matrix (image) with ~ square binning
    Requirement: need 1d array formats of
    self.data, self.qx_data, and self.qy_data
    where each one corresponds to z, x, or y axis values

    """
    # No qx or qy given in a vector format
    if self.qx_data is None or self.qy_data is None \
            or self.qx_data.ndim != 1 or self.qy_data.ndim != 1:
        return self.x_bins, self.y_bins, plottable

    # maximum # of loops to fillup_pixels
    # otherwise, loop could never stop depending on data
    max_loop = 1
    # get the x and y_bin arrays.
    x_bins, y_bins = _get_bins(self)
    # set zero to None

    #Note: Can not use scipy.interpolate.Rbf:
    # 'cause too many data points (>10000)<=JHC.
    # 1d array to use for weighting the data point averaging
    #when they fall into a same bin.
    weights_data = np.ones([self.data.size])
    # get histogram of ones w/len(data); this will provide
    #the weights of data on each bins
    weights, xedges, yedges = np.histogram2d(x=self.qy_data,
                                             y=self.qx_data,
                                             bins=[y_bins, x_bins],
                                             weights=weights_data)
    # get histogram of data, all points into a bin in a way of summing
    image, xedges, yedges = np.histogram2d(x=self.qy_data,
                                           y=self.qx_data,
                                           bins=[y_bins, x_bins],
                                           weights=plottable)
    # Now, normalize the image by weights only for weights>1:
    # If weight == 1, there is only one data point in the bin so
    # that no normalization is required.
    image[weights > 1] = image[weights > 1] / weights[weights > 1]
    # Set image bins w/o a data point (weight==0) as None (was set to zero
    # by histogram2d.)
    image[weights == 0] = None

    # Fill empty bins with 8 nearest neighbors only when at least
    #one None point exists
    loop = 0

    # do while loop until all vacant bins are filled up up
    #to loop = max_loop
    while (weights == 0).any():
        if loop >= max_loop:  # this protects never-ending loop
            break
        image = _fillup_pixels(image=image, weights=weights)
        loop += 1

    return x_bins, y_bins, image

def _get_bins(self):
    """
    get bins
    set x_bins and y_bins into self, 1d arrays of the index with
    ~ square binning
    Requirement: need 1d array formats of
    self.qx_data, and self.qy_data
    where each one corresponds to  x, or y axis values
    """
    # find max and min values of qx and qy
    xmax = self.qx_data.max()
    xmin = self.qx_data.min()
    ymax = self.qy_data.max()
    ymin = self.qy_data.min()

    # calculate the range of qx and qy: this way, it is a little
    # more independent
    x_size = xmax - xmin
    y_size = ymax - ymin

    # estimate the # of pixels on each axes
    npix_y = int(np.floor(np.sqrt(len(self.qy_data))))
    npix_x = int(np.floor(len(self.qy_data) / npix_y))

    # bin size: x- & y-directions
    xstep = x_size / (npix_x - 1)
    ystep = y_size / (npix_y - 1)

    # max and min taking account of the bin sizes
    xmax = xmax + xstep / 2.0
    xmin = xmin - xstep / 2.0
    ymax = ymax + ystep / 2.0
    ymin = ymin - ystep / 2.0

    # store x and y bin centers in q space
    x_bins = np.linspace(xmin, xmax, npix_x)
    y_bins = np.linspace(ymin, ymax, npix_y)

    return x_bins, y_bins

def _fillup_pixels(image=None, weights=None):
    """
    Fill z values of the empty cells of 2d image matrix
    with the average over up-to next nearest neighbor points

    :param image: (2d matrix with some zi = None)

    :return: image (2d array )

    :TODO: Find better way to do for-loop below

    """
    # No image matrix given
    if image is None or np.ndim(image) != 2 \
            or np.isfinite(image).all() \
            or weights is None:
        return image
    # Get bin size in y and x directions
    len_y = len(image)
    len_x = len(image[1])
    temp_image = np.zeros([len_y, len_x])
    weit = np.zeros([len_y, len_x])
    # do for-loop for all pixels
    for n_y in range(len(image)):
        for n_x in range(len(image[1])):
            # find only null pixels
            if weights[n_y][n_x] > 0 or np.isfinite(image[n_y][n_x]):
                continue
            else:
                # find 4 nearest neighbors
                # check where or not it is at the corner
                if n_y != 0 and np.isfinite(image[n_y - 1][n_x]):
                    temp_image[n_y][n_x] += image[n_y - 1][n_x]
                    weit[n_y][n_x] += 1
                if n_x != 0 and np.isfinite(image[n_y][n_x - 1]):
                    temp_image[n_y][n_x] += image[n_y][n_x - 1]
                    weit[n_y][n_x] += 1
                if n_y != len_y - 1 and np.isfinite(image[n_y + 1][n_x]):
                    temp_image[n_y][n_x] += image[n_y + 1][n_x]
                    weit[n_y][n_x] += 1
                if n_x != len_x - 1 and np.isfinite(image[n_y][n_x + 1]):
                    temp_image[n_y][n_x] += image[n_y][n_x + 1]
                    weit[n_y][n_x] += 1
                # go 4 next nearest neighbors when no non-zero
                # neighbor exists
                if n_y != 0 and n_x != 0 and \
                        np.isfinite(image[n_y - 1][n_x - 1]):
                    temp_image[n_y][n_x] += image[n_y - 1][n_x - 1]
                    weit[n_y][n_x] += 1
                if n_y != len_y - 1 and n_x != 0 and \
                        np.isfinite(image[n_y + 1][n_x - 1]):
                    temp_image[n_y][n_x] += image[n_y + 1][n_x - 1]
                    weit[n_y][n_x] += 1
                if n_y != len_y and n_x != len_x - 1 and \
                        np.isfinite(image[n_y - 1][n_x + 1]):
                    temp_image[n_y][n_x] += image[n_y - 1][n_x + 1]
                    weit[n_y][n_x] += 1
                if n_y != len_y - 1 and n_x != len_x - 1 and \
                        np.isfinite(image[n_y + 1][n_x + 1]):
                    temp_image[n_y][n_x] += image[n_y + 1][n_x + 1]
                    weit[n_y][n_x] += 1

    # get it normalized
    ind = (weit > 0)
    image[ind] = temp_image[ind] / weit[ind]

    return image


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
