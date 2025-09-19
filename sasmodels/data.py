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
from functools import wraps

import numpy as np  # type: ignore
from numpy import cos, ma, pi, sin, sqrt

# pylint: disable=unused-import
try:
    from typing import Callable, Optional, Union
    Data = Union["Data1D", "Data2D", "SesansData"]
    OptArray = Optional[np.ndarray]
    OptLimits = Optional[tuple[float, float]]
    OptString = Optional[str]
except ImportError:
    pass
# pylint: enable=unused-import

def load_data(filename, index=0):
    # type: (str, int) -> Data
    """
    Load data using a sasview loader.
    """
    try:
        from sasdata.dataloader.loader import Loader  # type: ignore
    except ImportError as ie:
        raise ImportError(f"{ie.name} is not available. Add sasdata to the python path.")
    loader = Loader()
    filename = str(filename)  # In case a Path was given.
    # Allow for one part in multipart file
    if '[' in filename:
        filename, indexstr = filename[:-1].split('[')
        index = int(indexstr)
    datasets = loader.load(filename)
    if not datasets:  # None or []
        raise OSError("Data %r could not be loaded" % filename)
    if not isinstance(datasets, list):
        datasets = [datasets]
    for data in datasets:
        if getattr(data, 'isSesans', False):
            pass
        elif hasattr(data, 'x'):
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
    try:
        from sasdata.data_util.manipulations import Ringcut
    except ImportError as ie:
        raise ImportError(f"{ie.name} is not available. Add sasdata to the python path.")
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
    try:
        from sasdata.data_util.manipulations import Boxcut
    except ImportError as ie:
        raise ImportError(f"{ie.name} is not available. Add sasdata to the python path.")
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
    try:
        from sasdata.data_util.manipulations import Boxcut
    except ImportError as ie:
        raise ImportError(f"{ie.name} is not available. Add sasdata to the python path.")
    data.mask += \
        Boxcut(x_min=-np.inf, x_max=np.inf, y_min=-np.inf, y_max=cutoff)(data)


def _as_numpy(data):
    return None if data is None else np.asarray(data)

class Data1D:
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
        # type: (OptArray, OptArray, OptArray, OptArray) -> None
        self.x, self.dx = _as_numpy(x), _as_numpy(dx)
        self.y, self.dy = _as_numpy(y), _as_numpy(dy)
        self.dxl = None
        self.filename = None
        self.qmin = self.x.min() if self.x is not None else np.nan
        self.qmax = self.x.max() if self.x is not None else np.nan
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
        self.lam = None  # type: OptArray
        self.xaxis("SE length", "A")
        self.yaxis("log(P)/(t L^2)", "1/A^2 1/cm")

class Data2D:
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
        # type: (OptArray, OptArray, OptArray, OptArray, OptArray, OptArray) -> None
        self.qx_data, self.dqx_data = _as_numpy(x), _as_numpy(dx)
        self.qy_data, self.dqy_data = _as_numpy(y), _as_numpy(dy)
        self.data, self.err_data = _as_numpy(z), _as_numpy(dz)
        self.mask = (np.isnan(z) if z is not None
                     else np.zeros_like(x, dtype='bool') if x is not None
                     else None)
        self.q_data = np.sqrt(self.qx_data**2 + self.qy_data**2)
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


class Vector:
    """
    3-space vector of *x*, *y*, *z*
    """
    def __init__(self, x=None, y=None, z=None):
        # type: (float, float, Optional[float]) -> None
        self.x, self.y, self.z = x, y, z

class Detector:
    """
    Detector attributes.
    """
    def __init__(self, pixel_size=(None, None), distance=None):
        # type: (tuple[float, float], float) -> None
        self.pixel_size = Vector(*pixel_size)
        self.distance = distance

class Source:
    """
    Beam attributes.
    """
    def __init__(self):
        # type: () -> None
        self.wavelength = np.nan
        self.wavelength_unit = "A"

class Sample:
    """
    Sample attributes.
    """
    def __init__(self):
        # type: () -> None
        pass

def empty_sesans(z, wavelength=None, zacceptance=None):
    data = SesansData(x=z, y=None, dx=None, dy=None)
    data.filename = "fake data"
    DEFAULT_WAVELENGTH = 5
    if wavelength is None:
        wavelength = DEFAULT_WAVELENGTH
    if np.isscalar(wavelength):
        wavelength = np.full_like(z, wavelength)
    if zacceptance is None:
        zacceptance = (np.pi/2, 'radians')
    source = Source()
    source.wavelength = wavelength
    source.wavelength_unit = "A"
    sample = Sample()
    sample.zacceptance = zacceptance
    data.source = source
    data.sample = sample
    return data

def empty_data1D(q, resolution=0.0, L=0., dL=0.):
    # type: (np.ndarray, float, float, float) -> Data1D
    r"""
    Create empty 1D data using the given *q* as the x value.

    rms *resolution* $\Delta q/q$ defaults to 0.0. Note that this is
    expressed as a fraction rather than a percentage.

    If wavelength *L* and rms wavelength divergence *dL* are given, then
    angle *theta* is infered from $q = 4\pi/\lambda \sin(\theta)$, and
    the *resolution* term applies to rms $\Delta \theta/\theta$ instead
    of $\Delta q/q$. Resolution $\Delta q$ is then calculated from wavelength
    and angular resolution.
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
    # type: (np.ndarray, OptArray, float) -> Data2D
    """
    Create empty 2D data using the given mesh.

    If *qy* is missing, create a square mesh with *qy=qx*.

    *resolution* dq/q defaults to 0.0.
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


def plot_data(data, view=None, limits=None, backend='matplotlib'):
    # type: (Data, str, OptLimits, str) -> None
    """
    Plot data loaded by the sasview loader.

    *data* is a sasview data object, either 1D, 2D or SESANS.

    *view* is log, linear or normed.

    *limits* sets the intensity limits on the plot; if None then the limits
    are inferred from the data.
    """
    # Note: kind of weird using the plot result functions to plot just the
    # data, but they already handle the masking and graph markup already, so
    # do not repeat.
    if hasattr(data, 'isSesans') and data.isSesans:
        _plot_result_sesans(data, None, None, view, use_data=True, limits=limits)
    elif hasattr(data, 'qx_data') and not getattr(data, 'radial', False):
        _plot_result2D(data, None, None, view, use_data=True, limits=limits)
    else:
        return _plot_result1D(
            data, None, None, view, use_data=True,
            limits=limits, backend=backend)


def plot_theory(data, theory, resid=None, view=None, use_data=True,
                limits=None, Iq_calc=None, label='theory', backend='matplotlib'):
    # type: (Data, OptArray, OptArray, OptString, bool, OptLimits, OptArray, str, str) -> None
    """
    Plot theory calculation.

    *data* is needed to define the graph properties such as labels and
    units, and to define the data mask.

    *theory* is a matrix of the same shape as the data.

    *view* is log, linear or normed

    *use_data* is True if the data should be plotted as well as the theory.

    *limits* sets the intensity limits on the plot; if None then the limits
    are inferred from the data. If (-inf, inf) then use auto limits.

    *Iq_calc* is the raw theory values without resolution smearing

    *label* is the label for the theory curve
    """
    if limits is not None and np.isinf(limits[0]):
        limits = None
    if hasattr(data, 'isSesans') and data.isSesans:
        fig = _plot_result_sesans(data, theory, resid, view, use_data=True, limits=limits, label=label, backend=backend)
    elif hasattr(data, 'qx_data') and not getattr(data, 'radial', False):
        fig = _plot_result2D(data, theory, resid, view, use_data, limits=limits, label=label, backend=backend)
    else:
        fig = _plot_result1D(
            data, theory, resid, view, use_data,
            limits=limits, Iq_calc=Iq_calc, label=label, backend=backend)
    return fig

def protect(func):
    # type: (Callable) -> Callable
    """
    Decorator to wrap calls in an exception trapper which prints the
    exception and continues.  Keyboard interrupts are ignored.
    """
    @wraps(func)
    def wrapper(*args, **kw):
        """
        Trap and print errors from function.
        """
        try:
            return func(*args, **kw)
        except Exception:
            # TODO: fix traceback printing
            # The first bit prints the traceback up to the protect call
            # The remaining bit prints the traceback within the protect call
            # Ideally we would join these into one traceback
            #print("Traceback (most recent call last):")
            stack = traceback.extract_stack(limit=4)
            lines = traceback.format_list(stack[:-1])
            print("".join(lines)[:-1])
            #print(f"{exc.__class__.__name__}: {exc}")
            traceback.print_exc()

    return wrapper

def get_data_label(data):
    # TODO: what do we want as data label?
    # example_data/1d_data/latex_smeared.xml[0]
    #     data.filename: "latex_smeared.xml"
    #     data.run: ["latex_sans"]
    #     data.title: "latex particles 0.5micron diameter in D2O Qdev"
    # example_data/1d_data/latex_smeared.xml[1]
    #     data.filename: "latex_smeared.xml"
    #     data.run: ["latex_usans"]
    #     data.title: "latex particles 0.5micron diameter in D2O slit"
    # For this dataset I'm using:
    #     data.run[0] > data.filename > data.title > "data"
    title = getattr(data, 'title', None)
    filename = getattr(data, 'filename', None)
    run = getattr(data, 'run', [])
    if isinstance(run, list) and run: # at least one run name
        run = run[0]
    return run if run else filename if filename else title if title else "data"

RESID_RATIO = 0.15
def _plot_result1D(
        data, theory, resid, view, use_data,
        limits=None, Iq_calc=None, label='theory',
        backend='matplotlib'):
    # type: (Data1D, OptArray, OptArray, str, bool, OptLimits, OptArray, str, str) -> None
    """
    Plot the data and residuals for 1D data.
    """
    # Default to 'log' view
    if view is None:
        view = 'log'

    if getattr(data, 'radial', False):
        data.x = data.q_data
        data.y = data.data

    use_data = use_data and data.y is not None
    use_theory = theory is not None
    use_resid = resid is not None
    use_calc = Iq_calc is not None and len(Iq_calc) == 2  # (q, Iq)
    use_calc2D = Iq_calc is not None and len(Iq_calc) == 3  # (qx, qy, Iq)

    def ytransform(x):
        return 1e8*x**4 if view == 'q4' else 1.0

    graph_ratio = 1.0 if resid is None else 1 - RESID_RATIO

    xlabel = "q/Å"
    ylabel = '10⁸q⁴I(q)' if view == 'q4' else 'I(q)'
    xscale = 'log'
    yscale = 'log' if view == 'log' else 'linear'

    if use_data:
        data_label = get_data_label(data)
        xdata = data.x # not a copy
        if (xdata <= 0).any():
            xscale = 'linear'
        ydata = ma.masked_array(data.y*ytransform(xdata), data.mask.copy()) # copy
        ydata[~np.isfinite(ydata)] = ma.masked
        if yscale == 'log':
            ydata[ydata <= 0] = ma.masked
        yerr = data.dy * ytransform(xdata)

    if use_theory:
        # TODO: provide model name to the plotter
        theory_label = label
        xtheory = data.x[data.mask == 0] #copy
        ytheory = ma.masked_array(theory*ytransform(xtheory)) # copy
        if xscale == 'log':
            xtheory[xtheory <= 0] = ma.masked
        ytheory[~np.isfinite(ytheory)] = ma.masked
        if yscale == 'log':
            ytheory[ytheory <= 0] = ma.masked

    if use_calc:
        calc_label = 'no resolution'
        xcalc, ycalc = Iq_calc # not a copy
        ycalc = ma.masked_array(ycalc*ytransform(xcalc)) # copy
        if xscale == 'log':
            xcalc = xcalc.copy() # copy
            xcalc[xcalc <= 0] = ma.masked
        ycalc[~np.isfinite(ycalc)] = ma.masked
        if yscale == 'log':
            ycalc[ycalc <= 0] = ma.masked

    if use_resid:
        resid_label = '(I(q) - y)/Δy'
        xresid = data.x[data.mask == 0] # copy
        yresid = ma.masked_array(resid.copy()) # copy
        yresid[~np.isfinite(yresid)] = ma.masked


    if backend == 'plotly':
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        rows = 2 if graph_ratio > 0 else 1
        heights = [graph_ratio, 1 - graph_ratio] if graph_ratio > 0 else None
        fig = make_subplots(
            cols=1, rows=rows, row_heights=heights,
            # TODO: x-axis sharing isn't working
            shared_xaxes=True,
            vertical_spacing=0.01)

        if use_data:
            fig.add_trace(
                go.Scatter(
                    x=xdata, y=ydata, error_y=dict(type='data', array=yerr),
                    name=data_label, mode='markers'),
                row=1, col=1)

        if use_theory:
            fig.add_trace(
                go.Scatter(x=xtheory, y=ytheory, name=theory_label, mode='lines'),
                row=1, col=1)

        if use_calc:  # (q, Iq), not (qx, qy, Iqxy)
            fig.add_trace(
                go.Scatter(x=xcalc, y=ycalc, name=calc_label, mode='lines', line=dict(dash='dot')),
                row=1, col=1)

        if limits is not None:
            fig.update_yaxes(range=limits, row=1, col=1)

        fig.update_xaxes(type=xscale, row=1, col=1)
        fig.update_yaxes(type=yscale, row=1, col=1)
        fig.update_layout(yaxis_title=ylabel)

        if use_resid:
            fig.add_trace(
                go.Scatter(x=xtheory, y=yresid, mode='markers', name=resid_label),
                row=2, col=1)
            fig.add_hline(y=0, line=dict(color='black'), row=2, col=1)
            fig.add_hline(y=1, line=dict(color='black', dash='dot'), row=2, col=1)
            fig.add_hline(y=-1, line=dict(color='black', dash='dot'), row=2, col=1)
            fig.update_yaxes(title_text='resid', type='linear', row=2, col=1)
            fig.update_xaxes(title_text=xlabel, row=2, col=1)
        else:
            fig.update_xaxes(title_text=xlabel, row=1, col=1)

        # Configure legend with automatic positioning
        fig.update_layout(
            legend=dict(
                x=0.99, y=0.99,  # top right
                yanchor="top",   # Anchor to middle of legend
                xanchor="right",   # Anchor to center of legend
                orientation="v",    # Horizontal layout
                bordercolor="black",
                borderwidth=1,
                bgcolor="rgba(255,255,255,0.8)",  # Semi-transparent background
                #draggable=True      # Allow manual adjustment if needed
            ),
            margin=dict(l=50, r=50, t=50, b=50)  # Add space for legend
        )

        margin = fig.layout.margin
        margin["t"] *= 1.2
        margin["r"] /= 4
        fig.update_layout(margin=margin)

        return fig

    else: # backend defaults to matplotlib
        import matplotlib.gridspec as gridspec
        import matplotlib.pyplot as plt

        # Note: matplotlib inset plots do not work properly in mpld3, so
        # use the more cumbersome gridspec interface directly.
        #   ax0 = plt.gca()  # for inset
        fig = plt.gcf()
        if use_resid:
            if fig.axes:
                ax = fig.gca()
                gs = gridspec.GridSpecFromSubplotSpec(
                    2, 1, subplot_spec=ax.get_subplotspec(),
                    height_ratios=[graph_ratio, 1 - graph_ratio])
            else:
                gs = gridspec.GridSpec(
                    2, 1, height_ratios=[graph_ratio, 1 - graph_ratio])
            ax0 = fig.add_subplot(gs[0])
        else:
            ax0 = fig.gca()

        if use_data:
            ax0.errorbar(xdata, ydata, yerr=yerr, fmt='.', label=data_label)

        if use_theory:
            ax0.plot(xtheory, ytheory, '-', label=theory_label)

        if use_calc:  # (q, Iq), not (qx, qy, Iqxy)
            ax0.plot(xcalc, ycalc, ':', alpha=0.8, label=calc_label)

        if limits is not None:
            ax0.set_ylim(*limits)

        ax0.set_xscale(xscale)
        ax0.set_yscale(yscale)
        ax0.set_ylabel(ylabel)
        ax0.legend()
        if use_resid:
            # Inset plot would use the following
            #   location = [0, -resid_ratio/(1-resid_ratio), 1, resid_ratio/(1-resid_ratio)]
            #   ax1 = ax0.inset_axes(location, transform=ax0.transAxes, sharex=ax0)
            # TODO: sharex isn't working right for mpld3; it is also sharing y
            ax1 = fig.add_subplot(gs[1]) #, sharex=ax0)
            ax1.plot(xresid, yresid, '.', label=resid_label)
            # Broken in mpld3
            #ax1.axhline(0, color='black')
            #ax1.axhline(1, color='black', linestyle=':')
            #ax1.axhline(-1, color='black', linestyle=':')
            xlim = min(xresid), max(xresid)
            #ax1.hlines([0, 1, -1], xlim[0], xlim[1], linestyles=['-',':',':'], color='black', alpha=0.5)
            ax1.plot(xlim, [0, 0], 'k-', label="_none_", alpha=0.5)
            ax1.plot(xlim, [-1, -1], 'k:', label="_none_", alpha=0.5)
            ax1.plot(xlim, [1, 1], 'k:', label="_none_", alpha=0.5)
            # ax1.set_ylabel('(I(q)-y)/Δy')
            ax1.set_ylabel('resid')
            ax1.set_yscale('linear')
            ax1.set_xlabel(xlabel)
            # TODO: reduce the tick density in the center when z scores are too high
            # The following are useful ticks for residuals. They are a little ugly in the
            # middle when the fit is bad, but at least this allows us to work around
            # flakiness in mpld3: The number of tick marks in the locator is not being
            # honoured by javascript.
            forward = [1, 3, 10, 30, 100, 300, 1000, 3000]
            reverse = [-v for v in forward[::-1]]
            ticks = [*reverse, *forward]
            ax1.yaxis.set_major_locator(plt.FixedLocator(ticks))
            #ax1.yaxis.set_major_locator(plt.LinearLocator(numticks=5))
            #ax1.yaxis.set_major_locator(plt.MaxNLocator(5))  # Show maximum 5 ticks
            plt.setp(ax0.get_xticklabels(), visible=False)
        else:
            ax0.set_xlabel(xlabel)

        if use_calc2D:  # (qx, qy, Iqxy)
            inset_calc_2D(ax0, Iq_calc)

        plt.tight_layout()


def inset_calc_2D(ax0, Iqxy_calc, alpha=0.8, portion=0.3, logz=True):
    qx, qy, Iqxy = Iqxy_calc
    # place in the top right
    location = [1-portion, 1-portion, portion, portion]
    ax_inset = ax0.inset_axes(location, transform=ax0.transAxes)
    ax_inset.pcolormesh(qx, qy, np.log10(Iqxy), alpha=alpha)
    ax_inset.set_aspect('equal')
    ax_inset.set_xlim(qx.min(), qx.max())
    ax_inset.set_ylim(qy.min(), qy.max())
    ax_inset.axis('off')
    return ax_inset



@protect
def _plot_result_sesans(data, theory, resid, view, use_data, limits=None, label='theory', backend="matplotlib"):
    # type: (SesansData, OptArray, OptArray, OptString, bool, OptLimits, str, str) -> None
    """
    Plot SESANS results.

    view is "normed" or not "normed"
    """
    import matplotlib.pyplot as plt  # type: ignore
    use_data = use_data and data.y is not None
    use_theory = theory is not None
    use_resid = resid is not None

    graph_ratio = 1.0 if resid is None else 1 - RESID_RATIO

    normed = (view == "normed")
    #normed = True
    offset, scale = 0, 1
    if normed and use_theory:
        offset, scale = theory[-1], theory[0] - theory[-1]

    is_tof = data.lam is not None and (data.lam != data.lam[0]).any()
    xdata = data.x
    ydata = data.y
    yerr = data.dy
    lamsq = data.source.wavelength**2
    if use_data:
        data_label = get_data_label(data)
        ydata = np.log10(ydata) / lamsq if is_tof else (ydata - offset) / scale
        yerr = yerr / ydata / lamsq if is_tof else yerr / scale
    if use_theory:
        ytheory = np.log10(theory) / lamsq if is_tof else (theory - offset) / scale

    theory_label = label
    xlabel = f"spin echo length ({data._xunit})"
    ylabel = "log(P)/(tλ²) / (Å² cm)"
    resid_label = 'resid'

    if backend == 'plotly':
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        rows = 2 if graph_ratio > 0 else 1
        heights = [graph_ratio, 1 - graph_ratio] if graph_ratio > 0 else None
        fig = make_subplots(
            cols=1, rows=rows, row_heights=heights,
            # TODO: x-axis sharing isn't working
            shared_xaxes=True,
            vertical_spacing=0.01)

        if use_data:
            fig.add_trace(
                go.Scatter(
                    x=xdata, y=ydata, error_y=dict(type='data', array=yerr),
                    name=data_label, mode='markers'),
                row=1, col=1)

        if use_theory:
            fig.add_trace(
                go.Scatter(x=xdata, y=ytheory, name=theory_label, mode='lines'),
                row=1, col=1)

        if limits is not None:
            fig.update_yaxes(range=limits, row=1, col=1)

        fig.update_xaxes(type='log', row=1, col=1)
        fig.update_yaxes(type='linear', row=1, col=1)
        fig.update_layout(yaxis_title=ylabel)

        if use_resid:
            fig.add_trace(
                go.Scatter(x=xdata, y=resid, mode='markers', name=resid_label),
                row=2, col=1)
            fig.add_hline(y=0, line=dict(color='black'), row=2, col=1)
            fig.add_hline(y=1, line=dict(color='black', dash='dot'), row=2, col=1)
            fig.add_hline(y=-1, line=dict(color='black', dash='dot'), row=2, col=1)
            fig.update_yaxes(title_text='resid', type='linear', row=2, col=1)
            fig.update_xaxes(title_text=xlabel, row=2, col=1)
        else:
            fig.update_xaxes(title_text=xlabel, row=1, col=1)

        # Configure legend with automatic positioning
        fig.update_layout(
            legend=dict(
                x=0.99, y=0.99,  # top right
                yanchor="top",   # Anchor to middle of legend
                xanchor="right",   # Anchor to center of legend
                orientation="v",    # Horizontal layout
                bordercolor="black",
                borderwidth=1,
                bgcolor="rgba(255,255,255,0.8)",  # Semi-transparent background
                #draggable=True      # Allow manual adjustment if needed
            ),
            margin=dict(l=50, r=50, t=50, b=50)  # Add space for legend
        )

        margin = fig.layout.margin
        margin["t"] *= 1.2
        margin["r"] /= 4
        fig.update_layout(margin=margin)

        return fig

    else: # backend == "matplotlib":
        import matplotlib.gridspec as gridspec
        import matplotlib.pyplot as plt

        # Note: matplotlib inset plots do not work properly in mpld3, so
        # use the more cumbersome gridspec interface directly.
        #   ax0 = fig.gca()  # for inset
        fig = plt.gcf()
        if use_resid:
            if fig.axes:
                ax = fig.gca()
                gs = gridspec.GridSpecFromSubplotSpec(
                    2, 1, subplot_spec=ax.get_subplotspec(),
                    height_ratios=[graph_ratio, 1 - graph_ratio])
            else:
                gs = gridspec.GridSpec(
                    2, 1, height_ratios=[graph_ratio, 1 - graph_ratio])
            ax0 = fig.add_subplot(gs[0])
        else:
            ax0 = fig.gca()

        if use_data:
            ax0.errorbar(xdata, ydata, yerr=yerr)
        if use_theory:
            ax0.plot(xdata, ytheory, '-')
        if limits is not None:
            ax0.set_ylim(*limits)
        ax0.set_ylabel(ylabel)
        ax0.set_xscale('log')

        if use_resid:
            ax1 = fig.add_subplot(gs[1]) #, sharex=ax0)
            ax1.plot(data.x, resid, 'x')
            ax1.set_xlabel(xlabel)
            ax1.set_ylabel(resid_label)
            ax1.set_xscale('log')

            # Broken in mpld3
            #ax1.axhline(0, color='black')
            #ax1.axhline(1, color='black', linestyle=':')
            #ax1.axhline(-1, color='black', linestyle=':')
            xlim = min(data.x), max(data.x)
            #ax1.hlines([0, 1, -1], xlim[0], xlim[1], linestyles=['-',':',':'], color='black', alpha=0.5)
            ax1.plot(xlim, [0, 0], 'k-', label="_none_", alpha=0.5)
            ax1.plot(xlim, [-1, -1], 'k:', label="_none_", alpha=0.5)
            ax1.plot(xlim, [1, 1], 'k:', label="_none_", alpha=0.5)

            # TODO: reduce the tick density in the center when z scores are too high
            # The following are useful ticks for residuals. They are a little ugly in the
            # middle when the fit is bad, but at least this allows us to work around
            # flakiness in mpld3: The number of tick marks in the locator is not being
            # honoured by javascript.
            forward = [1, 3, 10, 30, 100, 300, 1000, 3000]
            reverse = [-v for v in forward[::-1]]
            ticks = [*reverse, *forward]
            ax1.yaxis.set_major_locator(plt.FixedLocator(ticks))
            #ax1.yaxis.set_major_locator(plt.LinearLocator(numticks=5))
            #ax1.yaxis.set_major_locator(plt.MaxNLocator(5))  # Show maximum 5 ticks
            plt.setp(ax0.get_xticklabels(), visible=False)
        else:
            ax0.set_xlabel(xlabel)


@protect
def _plot_result2D(data, theory, resid, view, use_data, limits=None, label='theory', backend='matplotlib'):
    # type: (Data2D, OptArray, OptArray, str, bool, OptLimits, str, str) -> None
    """
    Plot the data and residuals for 2D data.
    """
    if view is None:
        view = 'log'

    use_data = use_data and data.data is not None
    use_theory = theory is not None
    use_resid = resid is not None
    num_plots = use_data + use_theory + use_resid

    # Note: leaving data masked since it may already be gridded.

    # Put theory and data on a common scale. Use the range on the data to set
    # the map, with 1/4 of the top of the error bars used as the basis of the
    # log scale colourmap to avoid spurious near-zero values after background
    # subtraction.
    if view == "q4":
        scale = 1e8*(data.qx_data**2+data.qy_data**2)**2
    else:
        scale = 1
    if view == "log" and limits is not None:
        limits = np.log10(limits[0]), np.log10(limits[1])
    if use_data and limits is None:
        # TODO: do we want log(q^4 Iq) = 4 log q + log Iq
        if view == 'log':
            upper = (abs(data.data) + data.err_data)[~data.mask]
            vmax = np.log10(upper.max())
            vmin = np.log10(upper[upper > 0].min()/4)
        elif view == 'q4':
            vdata = (data.data*scale)[~data.mask]
            vmin, vmax = vdata.min(), vdata.max()
        else:
            vdata = data.data[~data.mask]
            vmin, vmax = vdata.min(), vdata.max()

        # Allow 20% above the data to compare theory
        # For log data and q^4 data this is applying after the transform
        window = 0.2 if use_theory else 1
        vmin, vmax = vmin - (vmax-vmin)*(1+window), vmax + (vmax-vmin)*(1+window)
        if limits is None:
            limits = vmin, vmax

    if use_theory and limits is None:
        if view == 'log':
            vtheory = np.log10(theory)
        elif view == 'q4':
            vtheory = scale[~data.mask]*theory
        else:
            vtheory = theory

        if limits is None:
            limits = vtheory.min(), vtheory.max()

    # Use the full data limits for the residuals data
    rlimits = (resid.min(), resid.max()) if use_resid else None

    # _plot_2d_signal knows that the theory is only computed for the unmasked
    # data points. When using it to plot data, need to turn the data into a signal
    # by applying the mask. We need to keep using the mask so that data that is
    # stored as a 2D grid stays as a 2D grid, otherwise we need to apply a complicated
    # histogramming and interpolation algorithm to rebuild a grid.
    active_data = data.data[~data.mask]

    zlabel = 'log I(q)' if view == 'log' else '10⁸q⁴I(q)' if view == 'q4' else 'I(q)'
    theory_label = label
    data_label = f"{zlabel} data"
    resid_label = '(I(q) - y)/Δy'

    if backend == "plotly":
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        def add_trace(trace, row, col, label):
            colorbar = dict(
                title=dict(text=label, side='right'),
                yanchor="bottom",
                xanchor="left",
            )
            if num_plots == 1:
                trace.update(colorbar=colorbar)
                fig.add_trace(trace)
            else:
                trace.update(colorbar=dict(x=col/2-0.05, y=1-(row/2-0.05), len=0.45, **colorbar))
                fig.add_trace(trace, row=row, col=col)
            fig.update_xaxes(row=row, col=col, scaleanchor=f"y{row}", scaleratio=1)
            fig.update_yaxes(row=row, col=col, scaleanchor=f"x{col}", scaleratio=1)

        fig = make_subplots(rows=2, cols=2) if num_plots > 1 else go.Figure()
        if use_data:
            trace = _plot_2d_signal(data, active_data, view=view, limits=limits, backend=backend)
            add_trace(trace, 1, 1, data_label)

        if use_theory:
            trace = _plot_2d_signal(data, theory, view=view, limits=limits, backend=backend)
            add_trace(trace, 1, 2, theory_label)

        if use_resid:
            trace = _plot_2d_signal(data, resid, view='linear', limits=rlimits, backend=backend)
            add_trace(trace, 2, 2, resid_label)

        #from pprint import pprint; pprint(fig)
        return fig

    else: # backend == "matplotlib"
        import matplotlib.pyplot as plt  # type: ignore
        def maybe_hide_labels(row, col):
            if num_plots == 1:
                return
            if row == 1 and col == 2:
                plt.gca().set_xticks([])
                plt.xlabel('')
                plt.gca().set_yticks([])
                plt.ylabel('')

        # plot theory
        if use_theory:
            if num_plots > 1:
                plt.subplot(2, 2, 2)
            _plot_2d_signal(data, theory, label=theory_label, view=view, limits=limits, backend=backend)
            maybe_hide_labels(1, 2)

        # Plot data
        if use_data:
            if num_plots > 1:
                plt.subplot(2, 2, 1)
            _plot_2d_signal(data, active_data, label=data_label, view=view, limits=limits, backend=backend)
            maybe_hide_labels(1, 1)

        # plot resid
        if use_resid:
            if num_plots > 1:
                plt.subplot(2, 2, 4)
            _plot_2d_signal(data, resid, label=resid_label, view='linear', limits=rlimits, backend=backend)
            maybe_hide_labels(2, 2)


@protect
def _plot_2d_signal(data, signal, limits=None, view=None, label=None, backend='matplotlib'):
    # type: (Data2D, np.ndarray, Optional[float], Optional[float], str) -> tuple[float, float]
    """
    Plot the target value for the data.  This could be the data itself,
    the theory calculation, or the residuals.

    *scale* can be 'log' for log scale data, or 'linear'.
    """
    assert view is not None

    # Using signal evaluated at data.data[~data.mask], create an image on the regular
    # grid defined uniform vectors x, y
    # TODO: Divide range by 10 to convert from angstroms to nanometers
    xbins, ybins, masked_z = _build_image(data, signal)

    if view == 'log':
        masked_z[masked_z <= 0] = ma.masked
        masked_z = np.ma.log10(masked_z)
    elif view == 'q4':
        q4 = (ybins[:, None]**2 + xbins[None, :]**2)**2
        masked_z *= 1e8*q4

    if backend == 'plotly':
        import plotly.graph_objects as go

        xc = (xbins[1:] + xbins[:-1])/2
        yc = (ybins[1:] + ybins[:-1])/2
        trace = go.Heatmap(
            z=masked_z.filled(np.nan), y=yc, x=xc,
            zmin=limits[0], zmax=limits[1],
            colorscale='Viridis', showscale=True,
            )
        return trace

    else: # backend == 'matplotlib':
        import matplotlib.pyplot as plt  # type: ignore

        xmin, xmax = xbins[0], xbins[-1]
        ymin, ymax = ybins[0], ybins[-1]

        plt.imshow(masked_z,
                interpolation='nearest', aspect='equal', origin='lower',
                extent=[xmin, xmax, ymin, ymax],
                vmin=limits[0], vmax=limits[1])
        plt.xlabel("q_x/Å")
        plt.ylabel("q_y/Å")
        h = plt.colorbar(location='right')
        h.set_label(label)


# === The following is modified from sas.sasgui.plottools.PlotPanel
def _build_image(data, signal, radius=1):
    """
    Build a matrix for 2d plot from a vector
    Returns a matrix (image) with ~ square binning
    Requirement: need 1d array formats of
    self.data, self.qx_data, and self.qy_data
    where each one corresponds to z, x, or y axis values

    Holes in the image where there were no data points are filled
    using interpolation. Set *radius* to the maximum number
    of fill iterations.

    Returns qx, qy, masked image
    """

    # Check if data is already gridded. If so, return the gridding.
    # No qx or qy given in a vector format
    qx, qy, mask = data.qx_data, data.qy_data, data.mask
    if (qx is None or qy is None or qx.ndim != 1 or qy.ndim != 1):
        # build a masked array to return
        mask = mask if mask.ndim == 2 else mask.reshape((len(qy), len(qx)))
        image = np.zeros(mask.shape)
        image[~mask] = signal
        return data.x_bins, data.y_bins, ma.masked_array(image, mask)

    # get the x and y_bin arrays.
    x_bins, y_bins = _get_bins(qx, qy)
    # set zero to None

    # TODO: should masking apply before or after x, y step size estimate?
    # TODO: don't use number of points to estimate x, y step size
    # Apply masking to non-gridded data (signal is already masked)
    qx, qy = qx[~mask], qy[~mask]

    #Note: Can not use scipy.interpolate.Rbf:
    # 'cause too many data points (>10000)<=JHC.
    # get histogram of ones w/len(data); this will provide
    #the weights of data on each bins
    # TODO: check that x=qy, y=qx is correct.
    weights, _, _ = np.histogram2d(x=qy, y=qx, bins=[y_bins, x_bins])
    # get histogram of data, all points into a bin in a way of summing
    image, _, _ = np.histogram2d(x=qy, y=qx, bins=[y_bins, x_bins], weights=signal)
    # Now, normalize the image by weights only for weights>1:
    # If weight == 1, there is only one data point in the bin so
    # that no normalization is required.
    image[weights > 1] = image[weights > 1] / weights[weights > 1]
    # Set image bins w/o a data point (weight==0) as None (was set to zero
    # by histogram2d.)
    image[weights == 0] = np.nan

    # Fill empty bins with 8 nearest neighbors only when at least
    # one None point exists. Loop until all vacant bins are filled.
    for _ in range(radius):
        if (weights != 0).all():
            break
        image = _fillup_pixels(image=image, weights=weights)

    # Data is interpolated to fill the empty pixels
    return x_bins, y_bins, ma.masked_array(image, weights == 0)

def _get_bins(qx, qy):
    """
    get bins
    set x_bins and y_bins into self, 1d arrays of the index with
    ~ square binning
    Requirement: need 1d array formats for qx, qy
    where each one corresponds to x, or y axis values
    """
    # find max and min values of qx and qy
    xmin, xmax = qx.min(), qx.max()
    ymin, ymax = qy.min(), qy.max()

    # calculate the range of qx and qy: this way, it is a little
    # more independent
    x_size = xmax - xmin
    y_size = ymax - ymin

    # estimate the # of pixels on each axes
    npix_y = int(np.floor(np.sqrt(len(qx))))
    npix_x = int(np.floor(len(qy) / npix_y))

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
    if (image is None or np.ndim(image) != 2
            or np.isfinite(image).all()
            or weights is None):
        return image
    # Get bin size in y and x directions
    len_x, len_y = image.shape[1], image.shape[0]
    temp_image = np.zeros([len_y, len_x])
    weit = np.zeros([len_y, len_x])
    # do for-loop for all pixels
    for n_y in range(len_y):
        for n_x in range(len_x):
            # find only null pixels
            if weights[n_y][n_x] > 0 or np.isfinite(image[n_y][n_x]):
                continue
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
            if (n_y != 0 and n_x != 0
                    and np.isfinite(image[n_y - 1][n_x - 1])):
                temp_image[n_y][n_x] += image[n_y - 1][n_x - 1]
                weit[n_y][n_x] += 1
            if (n_y != len_y - 1 and n_x != 0
                    and np.isfinite(image[n_y + 1][n_x - 1])):
                temp_image[n_y][n_x] += image[n_y + 1][n_x - 1]
                weit[n_y][n_x] += 1
            if (n_y != len_y and n_x != len_x - 1
                    and np.isfinite(image[n_y - 1][n_x + 1])):
                temp_image[n_y][n_x] += image[n_y - 1][n_x + 1]
                weit[n_y][n_x] += 1
            if (n_y != len_y - 1 and n_x != len_x - 1
                    and np.isfinite(image[n_y + 1][n_x + 1])):
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
