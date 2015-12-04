import warnings

import numpy as np

from .core import load_model_definition, load_model, make_kernel
from .core import call_kernel, call_ER, call_VR
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
    """
    def _interpret_data(self, data, model):
        self._data = data
        self._model = model

        # interpret data
        if hasattr(data, 'lam'):
            self.data_type = 'sesans'
        elif hasattr(data, 'qx_data'):
            self.data_type = 'Iqxy'
        else:
            self.data_type = 'Iq'

        partype = model.info['partype']

        if self.data_type == 'sesans':
            q = sesans.make_q(data.sample.zacceptance, data.Rmax)
            self.index = slice(None, None)
            if data.y is not None:
                self.Iq = data.y
                self.dIq = data.dy
            #self._theory = np.zeros_like(q)
            q_vectors = [q]
        elif self.data_type == 'Iqxy':
            if not partype['orientation'] and not partype['magnetic']:
                raise ValueError("not 2D without orientation or magnetic parameters")
            q = np.sqrt(data.qx_data**2 + data.qy_data**2)
            qmin = getattr(data, 'qmin', 1e-16)
            qmax = getattr(data, 'qmax', np.inf)
            accuracy = getattr(data, 'accuracy', 'Low')
            self.index = ~data.mask & (q >= qmin) & (q <= qmax)
            if data.data is not None:
                self.index &= ~np.isnan(data.data)
                self.Iq = data.data[self.index]
                self.dIq = data.err_data[self.index]
            self.resolution = resolution2d.Pinhole2D(data=data, index=self.index,
                                                     nsigma=3.0, accuracy=accuracy)
            #self._theory = np.zeros_like(self.Iq)
            q_vectors = self.resolution.q_calc
        elif self.data_type == 'Iq':
            self.index = (data.x >= data.qmin) & (data.x <= data.qmax)
            if data.y is not None:
                self.index &= ~np.isnan(data.y)
                self.Iq = data.y[self.index]
                self.dIq = data.dy[self.index]
            if getattr(data, 'dx', None) is not None:
                q, dq = data.x[self.index], data.dx[self.index]
                if (dq>0).any():
                    self.resolution = resolution.Pinhole1D(q, dq)
                else:
                    self.resolution = resolution.Perfect1D(q)
            elif (getattr(data, 'dxl', None) is not None and
                          getattr(data, 'dxw', None) is not None):
                self.resolution = resolution.Slit1D(data.x[self.index],
                                                    width=data.dxh[self.index],
                                                    height=data.dxw[self.index])
            else:
                self.resolution = resolution.Perfect1D(data.x[self.index])

            #self._theory = np.zeros_like(self.Iq)
            q_vectors = [self.resolution.q_calc]
        else:
            raise ValueError("Unknown data type") # never gets here

        # Remember function inputs so we can delay loading the function and
        # so we can save/restore state
        self._kernel_inputs = [v for v in q_vectors]
        self._kernel = None

    def _set_data(self, Iq, noise=None):
        if noise is not None:
            self.dIq = Iq*noise*0.01
        dy = self.dIq
        y = Iq + np.random.randn(*dy.shape) * dy
        self.Iq = y
        if self.data_type == 'Iq':
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
            q_input = self._model.make_input(self._kernel_inputs)
            self._kernel = self._model(q_input)

        Iq_calc = call_kernel(self._kernel, pars, cutoff=cutoff)
        if self.data_type == 'sesans':
            result = sesans.hankel(self._data.x, self._data.lam * 1e-9,
                                   self._data.sample.thickness / 10,
                                   self._kernel_inputs[0], Iq_calc)
        else:
            result = self.resolution.apply(Iq_calc)
        return result


class DirectModel(DataMixin):
    def __init__(self, data, model, cutoff=1e-5):
        self.model = model
        self.cutoff = cutoff
        self._interpret_data(data, model)
        self.kernel = make_kernel(self.model, self._kernel_inputs)
    def __call__(self, **pars):
        return self._calc_theory(pars, cutoff=self.cutoff)
    def ER(self, **pars):
        return call_ER(self.model.info, pars)
    def VR(self, **pars):
        return call_VR(self.model.info, pars)
    def simulate_data(self, noise=None, **pars):
        Iq = self.__call__(**pars)
        self._set_data(Iq, noise=noise)

def demo():
    import sys
    from .data import empty_data1D, empty_data2D

    if len(sys.argv) < 3:
        print("usage: python -m sasmodels.direct_model modelname (q|qx,qy) par=val ...")
        sys.exit(1)
    model_name = sys.argv[1]
    call = sys.argv[2].upper()
    if call not in ("ER","VR"):
        try:
            values = [float(v) for v in call.split(',')]
        except:
            values = []
        if len(values) == 1:
            q, = values
            data = empty_data1D([q])
        elif len(values) == 2:
            qx,qy = values
            data = empty_data2D([qx],[qy])
        else:
            print("use q or qx,qy or ER or VR")
            sys.exit(1)
    else:
        data = empty_data1D([0.001])  # Data not used in ER/VR

    model_definition = load_model_definition(model_name)
    model = load_model(model_definition, dtype='single')
    calculator = DirectModel(data, model)
    pars = dict((k,float(v))
                for pair in sys.argv[3:]
                for k,v in [pair.split('=')])
    if call == "ER":
        print(calculator.ER(**pars))
    elif call == "VR":
        print(calculator.VR(**pars))
    else:
        Iq = calculator(**pars)
        print(Iq[0])

if __name__ == "__main__":
    demo()
