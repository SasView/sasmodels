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
direct parameter values and by :class:`bumps_model.Experiment` which wraps
the parameter values in boxes so that the user can set fitting ranges, etc.
on the individual parameters and send the model to the Bumps optimizers.
"""
from __future__ import print_function

import numpy as np  # type: ignore

# TODO: fix sesans module
from . import sesans  # type: ignore
from . import weights
from . import resolution
from . import resolution2d
from .details import make_kernel_args, dispersion_mesh
from .modelinfo import DEFAULT_BACKGROUND

# pylint: disable=unused-import
try:
    from typing import Optional, Dict, Tuple
except ImportError:
    pass
else:
    from .data import Data
    from .kernel import Kernel, KernelModel
    from .modelinfo import Parameter, ParameterSet
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
    #print("pars", list(zip(*mesh))[0])
    call_details, values, is_magnetic = make_kernel_args(calculator, mesh)
    #print("values:", values)
    return calculator(call_details, values, cutoff, is_magnetic)

def call_ER(model_info, pars):
    # type: (ModelInfo, ParameterSet) -> float
    """
    Call the model ER function using *values*.

    *model_info* is either *model.info* if you have a loaded model,
    or *kernel.info* if you have a model kernel prepared for evaluation.
    """
    if model_info.ER is None:
        return 1.0
    elif not model_info.parameters.form_volume_parameters:
        # handle the case where ER is provided but model is not polydisperse
        return model_info.ER()
    else:
        value, weight = _vol_pars(model_info, pars)
        individual_radii = model_info.ER(*value)
        return np.sum(weight*individual_radii) / np.sum(weight)


def call_VR(model_info, pars):
    # type: (ModelInfo, ParameterSet) -> float
    """
    Call the model VR function using *pars*.

    *model_info* is either *model.info* if you have a loaded model,
    or *kernel.info* if you have a model kernel prepared for evaluation.
    """
    if model_info.VR is None:
        return 1.0
    elif not model_info.parameters.form_volume_parameters:
        # handle the case where ER is provided but model is not polydisperse
        return model_info.VR()
    else:
        value, weight = _vol_pars(model_info, pars)
        whole, part = model_info.VR(*value)
        return np.sum(weight*part)/np.sum(weight*whole)


def call_profile(model_info, **pars):
    # type: (ModelInfo, ...) -> Tuple[np.ndarray, np.ndarray, Tuple[str, str]]
    """
    Returns the profile *x, y, (xlabel, ylabel)* representing the model.
    """
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
    # type: (ModelInfo, Dict[str, float], str, bool) -> List[Tuple[float, np.ndarray, np.ndarry]]
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

    #print("pars",[p.id for p in parameters.call_parameters])
    mesh = [_get_par_weights(p, values, active(p.name))
            for p in parameters.call_parameters]
    return mesh


def _get_par_weights(parameter, values, active=True):
    # type: (Parameter, Dict[str, float]) -> Tuple[float, np.ndarray, np.ndarray]
    """
    Generate the distribution for parameter *name* given the parameter values
    in *pars*.

    Uses "name", "name_pd", "name_pd_type", "name_pd_n", "name_pd_sigma"
    from the *pars* dictionary for parameter value and parameter dispersion.
    """
    value = float(values.get(parameter.name, parameter.default))
    npts = values.get(parameter.name+'_pd_n', 0)
    width = values.get(parameter.name+'_pd', 0.0)
    relative = parameter.relative_pd
    if npts == 0 or width == 0.0 or not active:
        # Note: orientation parameters have the viewing angle as the parameter
        # value and the jitter in the distribution, so be sure to set the
        # empty pd for orientation parameters to 0.
        pd = [value if relative or not parameter.polydisperse else 0.0], [1.0]
    else:
        limits = parameter.limits
        disperser = values.get(parameter.name+'_pd_type', 'gaussian')
        nsigma = values.get(parameter.name+'_pd_nsigma', 3.0)
        pd = weights.get_weights(disperser, npts, width, nsigma,
                                 value, limits, relative)
    return value, pd[0], pd[1]


def _vol_pars(model_info, values):
    # type: (ModelInfo, ParameterSet) -> Tuple[np.ndarray, np.ndarray]
    vol_pars = [_get_par_weights(p, values)
                for p in model_info.parameters.call_parameters
                if p.type == 'volume']
    #import pylab; pylab.plot(vol_pars[0][0],vol_pars[0][1]); pylab.show()
    dispersity, weight = dispersion_mesh(model_info, vol_pars)
    return dispersity, weight


def _make_sesans_transform(data):
    from sas.sascalc.data_util.nxsunit import Converter

    # Pre-compute the Hankel matrix (H)
    SElength = Converter(data._xunit)(data.x, "A")

    theta_max = Converter("radians")(data.sample.zacceptance)[0]
    q_max = 2 * np.pi / np.max(data.source.wavelength) * np.sin(theta_max)
    zaccept = Converter("1/A")(q_max, "1/" + data.source.wavelength_unit),

    Rmax = 10000000
    hankel = sesans.SesansTransform(data.x, SElength,
                                    data.source.wavelength,
                                    zaccept, Rmax)
    return hankel


class DataMixin(object):
    """
    DataMixin captures the common aspects of evaluating a SAS model for a
    particular data set, including calculating Iq and evaluating the
    resolution function.  It is used in particular by :class:`DirectModel`,
    which evaluates a SAS model parameters as key word arguments to the
    calculator method, and by :class:`bumps_model.Experiment`, which wraps the
    model and data for use with the Bumps fitting engine.  It is not
    currently used by :class:`sasview_model.SasviewModel` since this will
    require a number of changes to SasView before we can do it.

    :meth:`_interpret_data` initializes the data structures necessary
    to manage the calculations.  This sets attributes in the child class
    such as *data_type* and *resolution*.

    :meth:`_calc_theory` evaluates the model at the given control values.

    :meth:`_set_data` sets the intensity data in the data object,
    possibly with random noise added.  This is useful for simulating a
    dataset with the results from :meth:`_calc_theory`.
    """
    def _interpret_data(self, data, model):
        # type: (Data, KernelModel) -> None
        # pylint: disable=attribute-defined-outside-init

        self._data = data
        self._model = model

        # interpret data
        if hasattr(data, 'isSesans') and data.isSesans:
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
                  and getattr(data, 'dxw', None) is not None):
                res = resolution.Slit1D(data.x[index],
                                        qx_width=data.dxl[index],
                                        qy_width=data.dxw[index])
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

            res = resolution2d.Slit2D(data.x[index],
                                      qx_width=data.dxw[index],
                                      qy_width=data.dxl[index])
        else:
            raise ValueError("Unknown data type") # never gets here

        # Remember function inputs so we can delay loading the function and
        # so we can save/restore state
        self._kernel = None
        self.Iq, self.dIq, self.index = Iq, dIq, index
        self.resolution = res

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
        background = pars.get('background', DEFAULT_BACKGROUND)
        pars = pars.copy()
        pars['background'] = 0.

        Iq_calc = call_kernel(self._kernel, pars, cutoff=cutoff)
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

    *model* is a model calculator return from :func:`generate.load_model`

    *cutoff* is the polydispersity weight cutoff.
    """
    def __init__(self, data, model, cutoff=1e-5):
        # type: (Data, KernelModel, float) -> None
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
        return call_profile(self.model.info, **pars)

def main():
    # type: () -> None
    """
    Program to evaluate a particular model at a set of q values.
    """
    import sys
    from .data import empty_data1D, empty_data2D
    from .core import load_model_info, build_model

    if len(sys.argv) < 3:
        print("usage: python -m sasmodels.direct_model modelname (q|qx,qy) par=val ...")
        sys.exit(1)
    model_name = sys.argv[1]
    call = sys.argv[2].upper()
    if call != "ER_VR":
        try:
            values = [float(v) for v in call.split(',')]
        except ValueError:
            values = []
        if len(values) == 1:
            q, = values
            data = empty_data1D([q])
        elif len(values) == 2:
            qx, qy = values
            data = empty_data2D([qx], [qy])
        else:
            print("use q or qx,qy or ER or VR")
            sys.exit(1)
    else:
        data = empty_data1D([0.001])  # Data not used in ER/VR

    model_info = load_model_info(model_name)
    model = build_model(model_info)
    calculator = DirectModel(data, model)
    pars = dict((k, (float(v) if not k.endswith("_pd_type") else v))
                for pair in sys.argv[3:]
                for k, v in [pair.split('=')])
    if call == "ER_VR":
        ER = call_ER(model_info, pars)
        VR = call_VR(model_info, pars)
        print(ER, VR)
    else:
        Iq = calculator(**pars)
        print(Iq[0])

if __name__ == "__main__":
    main()
