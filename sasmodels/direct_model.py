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

import numpy as np

from .core import make_kernel
from .core import call_kernel, call_ER_VR
from . import sesans
from . import resolution
from . import resolution2d

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
        # pylint: disable=attribute-defined-outside-init

        self._data = data
        self._model = model

        # interpret data
        if hasattr(data, 'lam'):
            self.data_type = 'sesans'
        elif hasattr(data, 'qx_data'):
            self.data_type = 'Iqxy'
        elif getattr(data, 'oriented', False):
            self.data_type = 'Iq-oriented'
        else:
            self.data_type = 'Iq'

        partype = model.info['partype']

        if self.data_type == 'sesans':
            
            q = sesans.make_q(data.sample.zacceptance, data.Rmax)
            index = slice(None, None)
            res = None
            if data.y is not None:
                Iq, dIq = data.y, data.dy
            else:
                Iq, dIq = None, None
            #self._theory = np.zeros_like(q)
            q_vectors = [q]            
            q_mono = sesans.make_all_q(data)
        elif self.data_type == 'Iqxy':
            if not partype['orientation'] and not partype['magnetic']:
                raise ValueError("not 2D without orientation or magnetic parameters")
            q = np.sqrt(data.qx_data**2 + data.qy_data**2)
            qmin = getattr(data, 'qmin', 1e-16)
            qmax = getattr(data, 'qmax', np.inf)
            accuracy = getattr(data, 'accuracy', 'Low')
            index = ~data.mask & (q >= qmin) & (q <= qmax)
            if data.data is not None:
                index &= ~np.isnan(data.data)
                Iq = data.data[index]
                dIq = data.err_data[index]
            else:
                Iq, dIq = None, None
            res = resolution2d.Pinhole2D(data=data, index=index,
                                         nsigma=3.0, accuracy=accuracy)
            #self._theory = np.zeros_like(self.Iq)
            q_vectors = res.q_calc
            q_mono = []
        elif self.data_type == 'Iq':
            index = (data.x >= data.qmin) & (data.x <= data.qmax)
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
                                        qx_width=data.dxw[index],
                                        qy_width=data.dxl[index])
            else:
                res = resolution.Perfect1D(data.x[index])

            #self._theory = np.zeros_like(self.Iq)
            q_vectors = [res.q_calc]
            q_mono = []
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
            q_vectors = res.q_calc
            q_mono = []
        else:
            raise ValueError("Unknown data type") # never gets here

        # Remember function inputs so we can delay loading the function and
        # so we can save/restore state
        self._kernel_inputs = q_vectors
        self._kernel_mono_inputs = q_mono
        self._kernel = None
        self.Iq, self.dIq, self.index = Iq, dIq, index
        self.resolution = res

    def _set_data(self, Iq, noise=None):
        # pylint: disable=attribute-defined-outside-init
        if noise is not None:
            self.dIq = Iq*noise*0.01
        dy = self.dIq
        y = Iq + np.random.randn(*dy.shape) * dy
        self.Iq = y
        if self.data_type in ('Iq', 'Iq-oriented'):
            self._data.dy[self.index] = dy
            self._data.y[self.index] = y
        elif self.data_type == 'Iqxy':
            self._data.data[self.index] = y
        elif self.data_type == 'sesans':
            self._data.y[self.index] = y
        else:
            raise ValueError("Unknown model")

    def _calc_theory(self, pars, cutoff=0.0):
        if self._kernel is None:
            self._kernel = make_kernel(self._model, self._kernel_inputs)  # pylint: disable=attribute-dedata_type
            self._kernel_mono = (make_kernel(self._model, self._kernel_mono_inputs)
                                 if self._kernel_mono_inputs else None)

        Iq_calc = call_kernel(self._kernel, pars, cutoff=cutoff)
        # TODO: may want to plot the raw Iq for other than oriented usans
        self.Iq_calc = None
        if self.data_type == 'sesans':
            Iq_mono = (call_kernel(self._kernel_mono, pars, mono=True)
                       if self._kernel_mono_inputs else None)
            result = sesans.transform(self._data,
                                   self._kernel_inputs[0], Iq_calc, 
                                   self._kernel_mono_inputs, Iq_mono)
        else:
            result = self.resolution.apply(Iq_calc)
            if hasattr(self.resolution, 'nx'):
                self.Iq_calc = (
                    self.resolution.qx_calc, self.resolution.qy_calc,
                    np.reshape(Iq_calc, (self.resolution.ny, self.resolution.nx))
                )
        return result        


class DirectModel(DataMixin):
    """
    Create a calculator object for a model.

    *data* is 1D SAS, 2D SAS or SESANS data

    *model* is a model calculator return from :func:`generate.load_model`

    *cutoff* is the polydispersity weight cutoff.
    """
    def __init__(self, data, model, cutoff=1e-5):
        self.model = model
        self.cutoff = cutoff
        # Note: _interpret_data defines the model attributes
        self._interpret_data(data, model)

    def __call__(self, **pars):
        return self._calc_theory(pars, cutoff=self.cutoff)

    def ER_VR(self, **pars):
        """
        Compute the equivalent radius and volume ratio for the model.
        """
        return call_ER_VR(self.model.info, pars)

    def simulate_data(self, noise=None, **pars):
        """
        Generate simulated data for the model.
        """
        Iq = self.__call__(**pars)
        self._set_data(Iq, noise=noise)

def main():
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
        except:
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
    pars = dict((k, float(v))
                for pair in sys.argv[3:]
                for k, v in [pair.split('=')])
    if call == "ER_VR":
        print(calculator.ER_VR(**pars))
    else:
        Iq = calculator(**pars)
        print(Iq[0])

if __name__ == "__main__":
    main()
