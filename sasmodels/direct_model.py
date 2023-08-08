"""
Class interface to the model calculator.

Calling a model is somewhat non-trivial since the functions called depend
on the data type.  For 1D data the *Iq* kernel needs to be called, for
2D data the *Iqxy* kernel needs to be called, and for SESANS data the
*Iq* kernel needs to be called followed by a Hankel transform.  Before
the kernel is called an appropriate *q* calculation vector needs to be
constructed.  This is not the simple *q* vector where you have measured
the data since the resolution calculation will require values beyond the
range of the measured data.  After the calculation the resolution calculator
must be called to return the predicted value for each measured data point.

:class:`DirectModel` is a callable object that takes *parameter=value*
keyword arguments and returns the appropriate theory values for the data.

:class:`DataMixin` does the real work of interpreting the data and calling
the model calculator.  This is used by :class:`DirectModel`, which uses
direct parameter values and by :class:`.bumps_model.Experiment` which wraps
the parameter values in boxes so that the user can set fitting ranges, etc.
on the individual parameters and send the model to the Bumps optimizers.
"""
from __future__ import print_function, division

import os

import numpy as np  # type: ignore

# TODO: fix sesans module
from . import sesans  # type: ignore
from . import weights
from . import resolution
from . import resolution2d
from .details import make_kernel_args, dispersion_mesh
from .product import RADIUS_MODE_ID

# pylint: disable=unused-import
from typing import Optional, Dict, Tuple, List, Callable
from collections import OrderedDict
from .data import Data
from .kernel import Kernel, KernelModel
from .modelinfo import Parameter, ParameterSet, ModelInfo
# pylint: enable=unused-import

def call_kernel(calculator, pars, cutoff=0., mono=False):
    # type: (Kernel, ParameterSet, float, bool) -> np.ndarray
    """
    Call *kernel* returned from *model.make_kernel* with parameters *pars*.

    *cutoff* is the limiting value for the product of dispersion weights used
    to perform the multidimensional dispersion calculation more quickly at a
    slight cost to accuracy. The default value of *cutoff=0* integrates over
    the entire dispersion cube.  Using *cutoff=1e-5* can be 50% faster, but
    with an error of about 1%, which is usually less than the measurement
    uncertainty.

    *mono* is True if polydispersity should be set to none on all parameters.
    """
    mesh = get_mesh(calculator.info, pars, dim=calculator.dim, mono=mono)
    #print("in call_kernel: pars:", list(zip(*mesh))[0])
    call_details, values, is_magnetic = make_kernel_args(calculator, mesh)
    #print("in call_kernel: values:", values)
    return calculator(call_details, values, cutoff, is_magnetic)

def call_Fq(calculator, pars, cutoff=0., mono=False):
    # type: (Kernel, ParameterSet, float, bool) -> np.ndarray
    """
    Like :func:`call_kernel`, but returning F, F^2, R_eff, V_shell, V_form/V_shell.

    For solid objects V_shell is equal to V_form and the volume ratio is 1.

    Use parameter *radius_effective_mode* to select the effective radius
    calculation to use amongst the *radius_effective_modes* list given in the
    model.
    """
    R_eff_type = int(pars.pop(RADIUS_MODE_ID, 1.0))
    mesh = get_mesh(calculator.info, pars, dim=calculator.dim, mono=mono)
    #print("in call_Fq: pars", list(zip(*mesh))[0])
    call_details, values, is_magnetic = make_kernel_args(calculator, mesh)
    #print("in call_Fq: values:", values)
    return calculator.Fq(call_details, values, cutoff, is_magnetic, R_eff_type)

def call_profile(model_info, pars=None):
    # type: (ModelInfo, ParameterSet) -> Tuple[np.ndarray, np.ndarray, Tuple[str, str]]
    """
    Returns the profile *x, y, (xlabel, ylabel)* representing the model.
    """
    if pars is None:
        pars = {}
    args = {}
    for p in model_info.parameters.kernel_parameters:
        if p.length > 1:
            value = np.array([pars.get(p.id+str(j), p.default)
                              for j in range(1, p.length+1)])
        else:
            value = pars.get(p.id, p.default)
        args[p.id] = value
    x, y = model_info.profile(**args)
    return x, y, model_info.profile_axes

def get_mesh(model_info, values, dim='1d', mono=False):
    # type: (ModelInfo, Dict[str, float], str, bool) -> List[Tuple[float, np.ndarray, np.ndarray]]
    """
    Retrieve the dispersity mesh described by the parameter set.

    Returns a list of *(value, dispersity, weights)* with one tuple for each
    parameter in the model call parameters.  Inactive parameters return the
    default value with a weight of 1.0.
    """
    parameters = model_info.parameters
    if mono:
        active = lambda name: False
    elif dim == '1d':
        active = lambda name: name in parameters.pd_1d
    elif dim == '2d':
        active = lambda name: name in parameters.pd_2d
    else:
        active = lambda name: True

    #print("in get_mesh: pars:",[p.id for p in parameters.call_parameters])
    values = values.copy()
    mesh = [_pop_par_weights(p, values, active(p.name))
            for p in parameters.call_parameters]

    if values:
        raise TypeError(f"Unused parameters in call: {', '.join(values.keys())}")

    return mesh


def _pop_par_weights(parameter, values, active=True):
    # type: (Parameter, Dict[str, float], bool) -> Tuple[float, np.ndarray, np.ndarray]
    """
    Generate the distribution for parameter *name* given the parameter values
    in *pars*.

    Uses "name", "name_pd", "name_pd_type", "name_pd_n", "name_pd_sigma"
    from the *pars* dictionary for parameter value and parameter dispersion.
    """
    value = float(values.pop(parameter.name, parameter.default))
    if parameter.polydisperse:
        npts = values.pop(parameter.name+'_pd_n', 0)
        width = values.pop(parameter.name+'_pd', 0.0)
        nsigma = values.pop(parameter.name+'_pd_nsigma', 3.0)
        distribution = values.pop(parameter.name+'_pd_type', 'gaussian')
        relative = parameter.relative_pd
        if npts == 0 or width == 0.0 or not active:
            # Note: orientation parameters have the viewing angle as the parameter
            # value and the jitter in the distribution, so be sure to set the
            # empty pd for orientation parameters to 0.
            pd = [value if relative else 0.0], [1.0]
        else:
            limits = parameter.limits
            pd = weights.get_weights(distribution, npts, width, nsigma,
                                    value, limits, relative)
    else:
        pd = [value], [1.0]
    return value, pd[0], pd[1]


def _make_sesans_transform(data):
    # Pre-compute the Hankel matrix (H)
    SElength, SEunits = data.x, data._xunit
    wavelength, wunits = data.source.wavelength, data.source.wavelength_unit
    theta_max, theta_units = data.sample.zacceptance
    if SEunits != "A" or wunits != "A" or theta_units != "radians":
        try:
            from sasdata.data_util.nxsunit import Converter
        except ImportError as ie:
            raise ImportError(f"{ie.name} is not available. Add sasdata to the python path.")

        SElength = Converter("A")(SElength, units=SEunits)
        wavelength = Converter("A")(wavelength, units=wunits)
        theta_max = Converter("radian")(theta_max, units=theta_units)

    Rmax = 10000000
    zaccept = 2 * np.pi / np.max(wavelength) * np.sin(theta_max)
    hankel = sesans.SesansTransform(data.x, SElength, wavelength, zaccept, Rmax)
    return hankel


class DataMixin(object):
    """
    DataMixin captures the common aspects of evaluating a SAS model for a
    particular data set, including calculating Iq and evaluating the
    resolution function.  It is used in particular by :class:`DirectModel`,
    which evaluates a SAS model parameters as key word arguments to the
    calculator method, and by :class:`.bumps_model.Experiment`, which wraps the
    model and data for use with the Bumps fitting engine.  It is not
    currently used by :class:`.sasview_model.SasviewModel` since this will
    require a number of changes to SasView before we can do it.

    *_interpret_data* initializes the data structures necessary
    to manage the calculations.  This sets attributes in the child class
    such as *data_type* and *resolution*.

    *_calc_theory* evaluates the model at the given control values.

    *_set_data* sets the intensity data in the data object,
    possibly with random noise added.  This is useful for simulating a
    dataset with the results from *_calc_theory*.
    """
    def _interpret_data(self, data: Data, model: KernelModel) -> None:
        # not type: (Data, KernelModel) -> None
        # pylint: disable=attribute-defined-outside-init

        self._data = data
        self._model = model

        # interpret data
        if getattr(data, 'isSesans', False):
            self.data_type = 'sesans'
        elif hasattr(data, 'qx_data'):
            self.data_type = 'Iqxy'
        elif getattr(data, 'oriented', False):
            self.data_type = 'Iq-oriented'
        else:
            self.data_type = 'Iq'

        if self.data_type == 'sesans':
            res = _make_sesans_transform(data)
            index = slice(None, None)
            if data.y is not None:
                Iq, dIq = data.y, data.dy
            else:
                Iq, dIq = None, None
        elif self.data_type == 'Iqxy':
            #if not model.info.parameters.has_2d:
            #    raise ValueError("not 2D without orientation or magnetic parameters")
            q = np.sqrt(data.qx_data**2 + data.qy_data**2)
            qmin = getattr(data, 'qmin', 1e-16)
            qmax = getattr(data, 'qmax', np.inf)
            accuracy = getattr(data, 'accuracy', 'Low')
            index = (data.mask == 0) & (q >= qmin) & (q <= qmax)
            if data.data is not None:
                index &= ~np.isnan(data.data)
                Iq = data.data[index]
                dIq = data.err_data[index]
            else:
                Iq, dIq = None, None
            res = resolution2d.Pinhole2D(data=data, index=index,
                                         nsigma=3.0, accuracy=accuracy)
        elif self.data_type == 'Iq':
            index = (data.x >= data.qmin) & (data.x <= data.qmax)
            mask = getattr(data, 'mask', None)
            if mask is not None:
                index &= (mask == 0)
            if data.y is not None:
                index &= ~np.isnan(data.y)
                Iq = data.y[index]
                dIq = data.dy[index]
            else:
                Iq, dIq = None, None
            if getattr(data, 'dx', None) is not None:
                q, dq = data.x[index], data.dx[index]
                if (dq > 0).any():
                    res = resolution.Pinhole1D(q, dq)
                else:
                    res = resolution.Perfect1D(q)
            elif (getattr(data, 'dxl', None) is not None
                  or getattr(data, 'dxw', None) is not None):
                res = resolution.Slit1D(
                    data.x[index],
                    q_length=None if data.dxl is None else data.dxl[index],
                    q_width=None if data.dxw is None else data.dxw[index])
            else:
                res = resolution.Perfect1D(data.x[index])
        elif self.data_type == 'Iq-oriented':
            index = (data.x >= data.qmin) & (data.x <= data.qmax)
            if data.y is not None:
                index &= ~np.isnan(data.y)
                Iq = data.y[index]
                dIq = data.dy[index]
            else:
                Iq, dIq = None, None
            if (getattr(data, 'dxl', None) is None
                    or getattr(data, 'dxw', None) is None):
                raise ValueError("oriented sample with 1D data needs slit resolution")

            res = resolution2d.Slit2D(
                data.x[index],
                qx_width=data.dxw[index],
                qy_width=data.dxl[index])
        else:
            raise ValueError("Unknown data type") # never gets here

        # Remember function inputs so we can delay loading the function and
        # so we can save/restore state
        self._kernel = None
        self.Iq, self.dIq, self.index = Iq, dIq, index
        self.resolution = res
        self.results = None  # type: Optional[Callable[[], OrderedDict]]

    def _set_data(self, Iq, noise=None):
        # type: (np.ndarray, Optional[float]) -> None
        # pylint: disable=attribute-defined-outside-init
        if noise is not None:
            self.dIq = Iq*noise*0.01
        dy = self.dIq
        y = Iq + np.random.randn(*dy.shape) * dy
        self.Iq = y
        if self.data_type in ('Iq', 'Iq-oriented'):
            if self._data.y is None:
                self._data.y = np.empty(len(self._data.x), 'd')
            if self._data.dy is None:
                self._data.dy = np.empty(len(self._data.x), 'd')
            self._data.dy[self.index] = dy
            self._data.y[self.index] = y
        elif self.data_type == 'Iqxy':
            if self._data.data is None:
                self._data.data = np.empty_like(self._data.qx_data, 'd')
            if self._data.err_data is None:
                self._data.err_data = np.empty_like(self._data.qx_data, 'd')
            self._data.data[self.index] = y
            self._data.err_data[self.index] = dy
        elif self.data_type == 'sesans':
            if self._data.y is None:
                self._data.y = np.empty(len(self._data.x), 'd')
            self._data.y[self.index] = y
        else:
            raise ValueError("Unknown model")

    def _calc_theory(self, pars, cutoff=0.0):
        # type: (ParameterSet, float) -> np.ndarray
        if self._kernel is None:
            # TODO: change interfaces so that resolution returns kernel inputs
            # Maybe have resolution always return a tuple, or maybe have
            # make_kernel accept either an ndarray or a pair of ndarrays.
            kernel_inputs = self.resolution.q_calc
            if isinstance(kernel_inputs, np.ndarray):
                kernel_inputs = (kernel_inputs,)
            self._kernel = self._model.make_kernel(kernel_inputs)

        # Need to pull background out of resolution for multiple scattering
        default_background = self._model.info.parameters.common_parameters[1].default
        background = (
            pars.get('background', default_background)
            if self.data_type != 'sesans' else 0.)
        pars = pars.copy()
        pars['background'] = 0.

        Iq_calc = call_kernel(self._kernel, pars, cutoff=cutoff)
        self.results = getattr(self._kernel, 'results', None)
        # Storing the calculated Iq values so that they can be plotted.
        # Only applies to oriented USANS data for now.
        # TODO: extend plotting of calculate Iq to other measurement types
        # TODO: refactor so we don't store the result in the model
        self.Iq_calc = Iq_calc
        result = self.resolution.apply(Iq_calc)
        if hasattr(self.resolution, 'nx'):
            self.Iq_calc = (
                self.resolution.qx_calc, self.resolution.qy_calc,
                np.reshape(Iq_calc, (self.resolution.ny, self.resolution.nx))
            )
        return result + background


class DirectModel(DataMixin):
    """
    Create a calculator object for a model.

    *data* is 1D SAS, 2D SAS or SESANS data

    *model* is a model calculator return from :func:`.core.load_model`

    *cutoff* is the polydispersity weight cutoff.
    """
    def __init__(self, data: Data, model: KernelModel, cutoff: float=1e-5) -> None:
        # not type: (Data, KernelModel, float) -> None
        self.model = model
        self.cutoff = cutoff
        # Note: _interpret_data defines the model attributes
        self._interpret_data(data, model)

    def __call__(self, **pars):
        # type: (**float) -> np.ndarray
        return self._calc_theory(pars, cutoff=self.cutoff)

    def simulate_data(self, noise=None, **pars):
        # type: (Optional[float], **float) -> None
        """
        Generate simulated data for the model.
        """
        Iq = self.__call__(**pars)
        self._set_data(Iq, noise=noise)

    def profile(self, **pars):
        # type: (**float) -> None
        """
        Generate a plottable profile.
        """
        return call_profile(self.model.info, pars)

def test_reparameterize():
    # type: () -> None
    """Check simple reparameterized models will load and build"""
    from numpy import inf, pi
    from numpy.linalg import norm
    from .data import empty_data1D
    from .core import load_model_info, build_model, reparameterize
    parameters = [
        ["volume", "Ang^3", 1e5, [0, inf], "volume", "ellipsoid volume"],
        ["eccentricity", "", 1, [0, inf], "volume", "polar:equatorial radius"],
    ]
    translation = """
        Re = cbrt(volume/eccentricity/M_4PI_3)
        radius_polar = eccentricity*Re
        radius_equatorial = Re  # python style comments allowed
        """
    info = reparameterize('ellipsoid', parameters, translation)

    # Check translated parameters
    actual = [p.name for p in info.parameters.kernel_parameters]
    target = ["sld", "sld_solvent", "volume", "eccentricity", "theta", "phi"]
    assert target == actual, "%s != %s"%(target, actual)

    data = empty_data1D([0.01, 0.1])
    Rp, Re = 30, 50

    base_model = build_model(load_model_info('ellipsoid'), dtype='d')
    base_calculator = DirectModel(data, base_model)
    base_Iq = base_calculator(radius_polar=Rp, radius_equatorial=Re)

    derived_model = build_model(info, dtype='d')
    derived_calculator = DirectModel(data, derived_model)
    derived_Iq = derived_calculator(volume=4*pi/3*Rp*Re**2, eccentricity=Rp/Re)
    assert norm((base_Iq - derived_Iq)/base_Iq) < 1e-13

    # clean up the dll that was created in the cache as part of the test
    try:
        os.remove(derived_model.dllpath)
    except Exception:
        pass

    # Again using insert_after.
    insert_after = {
        '': 'eccentricity',
        'radius_equatorial': 'volume',
        }
    info = reparameterize('ellipsoid', parameters, translation,
                          insert_after=insert_after)
    # Check translated parameters
    actual = [p.name for p in info.parameters.kernel_parameters]
    target = ["eccentricity", "sld", "sld_solvent", "volume", "theta", "phi"]
    assert target == actual, "%s != %s"%(target, actual)

    # Check calculation.
    derived_model = build_model(info, dtype='d')
    derived_calculator = DirectModel(data, derived_model)
    derived_Iq = derived_calculator(volume=4*pi/3*Rp*Re**2, eccentricity=Rp/Re)
    assert norm((base_Iq - derived_Iq)/base_Iq) < 1e-13

    # clean up the dll that was created in the cache as part of the test
    try:
        os.remove(derived_model.dllpath)
    except Exception:
        pass

def _direct_calculate(model, data, pars):
    from .core import load_model_info, build_model
    model_info = load_model_info(model)
    kernel = build_model(model_info)
    calculator = DirectModel(data, kernel)
    return calculator(**pars)

def Iq(model, q, dq=None, ql=None, qw=None, **pars):
    """
    Compute I(q) for *model*. Resolution is *dq* for pinhole or *ql* and *qw*
    for slit geometry. Use 0 or None for infinite slits.

    Model is the name of a builtin or custom model, or a model expression, such
    as sphere+sphere for a mixture of spheres of different radii, or
    sphere@hardsphere for concentrated solutions where the dilute approximation
    no longer applies.

    Use additional keywords for model parameters, tagged with *_pd*, *_pd_n*,
    *_pd_nsigma*, *_pd_type* to set polydispersity parameters, or *_M0*,
    *_mphi*, *_mtheta* for magnetic parameters.

    This is not intended for use when the same I(q) is evaluated many times
    with different parameter values. For that you should set up the model
    with `model = build_model(load_model_info(model_name))`, set up a data
    object to define q values and resolution, then use
    `calculator = DirectModel(data, model)` to set up a calculator, or
    `problem = bumps.FitProblem(sasmodels.bumps_model.Experiment(data, model))`
    to define a fit problem for uses with the bumps optimizer. Data can be
    loaded using the `sasdata` package, or use one of the empty data generators
    from `sasmodels.data`.

    Models are cached. Custom models will not be reloaded even if the
    underlying files have changed. If you are using this in a long running
    application then you will need to call
    `sasmodels.direct_model._model_cache.clear()` to reset the cache and force
    custom model reload.
    """
    from .data import Data1D, _as_numpy
    data = Data1D(x=q, dx=dq)
    def broadcast(v):
        return (
            None if v is None
            else np.full(len(q), v) if np.isscalar(v)
            else _as_numpy(v))
    data.dxl, data.dxw = broadcast(ql), broadcast(qw)
    return _direct_calculate(model, data, pars)

def Iqxy(model, qx, qy, dqx=None, dqy=None, **pars):
    """
    Compute I(qx, qy) for *model*. Resolution is *dqx* and *dqy*.
    See :func:`Iq` for details on model and parameters.
    """
    from .data import Data2D
    data = Data2D(x=qx, y=qy, dx=dqx, dy=dqy)
    return _direct_calculate(model, data, pars)

def Gxi(model, xi, **pars):
    """
    Compute SESANS correlation G' = G(xi) - G(0) for *model*.
    See :func:`Iq` for details on model and parameters.
    """
    from .data import empty_sesans
    data = empty_sesans(z=xi)
    return _direct_calculate(model, data, pars)

def main():
    # type: () -> None
    """
    Program to evaluate a particular model at a set of q values.
    """
    import sys

    if len(sys.argv) < 3:
        print("usage: python -m sasmodels.direct_model modelname (q|qx,qy) par=val ...")
        sys.exit(1)
    model = sys.argv[1]
    call = sys.argv[2].upper()
    pars = dict((k, (float(v) if not k.endswith("_pd_type") else v))
                for pair in sys.argv[3:]
                for k, v in [pair.split('=')])
    try:
        values = [float(v) for v in call.split(',')]
    except ValueError:
        values = []
    if len(values) == 1:
        q, = values
        dq = dqw = dql = None
        #dq = [q*0.05] # 5% pinhole resolution
        #dqw, dql = [q*0.05], [1.0] # 5% horizontal slit resolution
        print(Iq(model, [q], dq=dq, qw=dqw, ql=dql, **pars)[0])
        #print(Gxi(model, [q], **pars)[0])
    elif len(values) == 2:
        qx, qy = values
        dq = None
        #dq = [0.005] # 5% pinhole resolution at q = 0.1
        print(Iqxy(model, [qx], [qy], dqx=dq, dqy=dq, **pars)[0])
    else:
        print("use q or qx,qy")
        sys.exit(1)

def test_simple_interface():
    def near(value, target):
        """Close enough in single precision"""
        #print(f"value: {value}, target: {target}")
        return np.allclose(value, target, rtol=1e-6, atol=0, equal_nan=True)
    # Note: target values taken from running main() on parameters.
    # Resolution was 5% dq/q.
    pars = dict(radius=200, background=0)  # default background=1e-3, scale=1
    # simple sphere in 1D (perfect, pinhole, slit)
    perfect_target = 0.6190146273894904
    assert near(Iq('sphere', [0.1], **pars), [perfect_target])
    assert near(Iq('sphere', [0.1], dq=[0.005], **pars), [2.3009224683980215])
    assert near(Iq('sphere', [0.1], qw=[0.005], ql=[1.0], **pars), [0.3663431784535172])
    # simple sphere in 2D (perfect, pinhole)
    assert near(Iqxy('sphere', [0.1], [0.1], **pars), [1.1771532874802199])
    assert near(Iqxy('sphere', [0.1], [0.1], dqx=[0.005], dqy=[0.005], **pars), 
        [0.8167780778578667])
    # sesans (no background or scale)
    assert near(Gxi('sphere', [100], **pars), [-0.19146959126623486])
    # Check that single point sesans matches value in an array
    xi = np.logspace(1, 3, 100)
    y = Gxi('sphere', xi, **pars)
    for k in (0, len(xi)//5, len(xi)//2, len(xi)-1):
        ysingle = Gxi('sphere', [xi[k]], **pars)[0]
        print(f"SESANS point check {k}: xi={xi[k]:.1f} single={ysingle:.4f} vector={y[k]:.4f}")
        assert abs((ysingle-y[k])/y[k]) < 0.1, "SESANS point value not matching vector value within 10%"
    # magnetic 2D
    pars = dict(radius=200, sld_M0=3, sld_mtheta=30)
    assert near(Iqxy('sphere', [0.1], [0.1], **pars), [1.5577852226925908])
    # polydisperse 1D
    pars = dict(
        radius=200, radius_pd=0.1, radius_pd_n=15, radius_pd_nsigma=2.5,
        radius_pd_type="uniform")
    assert near(Iq('sphere', [0.1], **pars), [2.703169824954617])
    # background and scale
    background, scale = 1e-4, 0.1
    pars = dict(radius=200, background=background, scale=scale)
    assert near(Iq('sphere', [0.1], **pars), [perfect_target*scale + background])


if __name__ == "__main__":
    import logging
    logging.disable(logging.ERROR)
    main()
    #test_simple_interface()
