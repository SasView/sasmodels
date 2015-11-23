"""
Wrap sasmodels for direct use by bumps.

:class:`Model` is a wrapper for the sasmodels kernel which defines a
bumps *Parameter* box for each kernel parameter.  *Model* accepts keyword
arguments to set the initial value for each parameter.

:class:`Experiment` combines the *Model* function with a data file loaded by the
sasview data loader.  *Experiment* takes a *cutoff* parameter controlling
how far the polydispersity integral extends.

"""

import datetime
import warnings

import numpy as np

from . import sesans
from . import weights
from .data import plot_theory
from .direct_model import DataMixin

# CRUFT: old style bumps wrapper which doesn't separate data and model
def BumpsModel(data, model, cutoff=1e-5, **kw):
    warnings.warn("Use of BumpsModel is deprecated.  Use bumps_model.Experiment instead.")
    model = Model(model, **kw)
    experiment = Experiment(data=data, model=model, cutoff=cutoff)
    for k in model._parameter_names:
        setattr(experiment, k, getattr(model, k))
    return experiment


class Model(object):
    def __init__(self, model, **kw):
        # lazy import; this allows the doc builder and nosetests to run even
        # when bumps is not on the path.
        from bumps.names import Parameter

        self._sasmodel = model
        partype = model.info['partype']

        pars = []
        for p in model.info['parameters']:
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


class Experiment(DataMixin):
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
        self.model = model
        self.cutoff = cutoff
        self._interpret_data(data, model._sasmodel)
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
            pars = dict((k, v.value) for k,v in self.model.parameters().items())
            self._cache['theory'] = self._calc_theory(pars, cutoff=self.cutoff)
            """
            if self._fn is None:
                q_input = self.model.kernel.make_input(self._kernel_inputs)
                self._fn = self.model.kernel(q_input)

            fixed_pars = [getattr(self.model, p).value for p in self._fn.fixed_pars]
            pd_pars = [self._get_weights(p) for p in self._fn.pd_pars]
            #print fixed_pars,pd_pars
            Iq_calc = self._fn(fixed_pars, pd_pars, self.cutoff)
            #self._theory[:] = self._fn.eval(pars, pd_pars)
            if self.data_type == 'sesans':
                result = sesans.hankel(self.data.x, self.data.lam * 1e-9,
                                       self.data.sample.thickness / 10,
                                       self._kernel_inputs[0], Iq_calc)
                self._cache['theory'] = result
            else:
                Iq = self.resolution.apply(Iq_calc)
                self._cache['theory'] = Iq
            """
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
        data, theory, resid = self._data, self.theory(), self.residuals()
        plot_theory(data, theory, resid, view)

    def simulate_data(self, noise=None):
        Iq = self.theory()
        self._set_data(Iq, noise)

    def save(self, basename):
        pass

    def remove_get_weights(self, name):
        """
        Get parameter dispersion weights
        """
        info = self.model.kernel.info
        relative = name in info['partype']['pd-rel']
        limits = info['limits'][name]
        disperser, value, npts, width, nsigma = [
            getattr(self.model, name + ext)
            for ext in ('_pd_type', '', '_pd_n', '_pd', '_pd_nsigma')]
        value, weight = weights.get_weights(
            disperser, int(npts.value), width.value, nsigma.value,
            value.value, limits, relative)
        return value, weight / np.sum(weight)

    def __getstate__(self):
        # Can't pickle gpu functions, so instead make them lazy
        state = self.__dict__.copy()
        state['_kernel'] = None
        return state

    def __setstate__(self, state):
        # pylint: disable=attribute-defined-outside-init
        self.__dict__ = state
