"""
Wrap sasmodels for direct use by bumps.

:class:`Model` is a wrapper for the sasmodels kernel which defines a
bumps *Parameter* box for each kernel parameter.  *Model* accepts keyword
arguments to set the initial value for each parameter.

:class:`Experiment` combines the *Model* function with a data file loaded by the
sasview data loader.  *Experiment* takes a *cutoff* parameter controlling
how far the polydispersity integral extends.

"""

import warnings

import numpy as np

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

def create_parameters(model_info, **kwargs):
    # lazy import; this allows the doc builder and nosetests to run even
    # when bumps is not on the path.
    from bumps.names import Parameter

    partype = model_info['partype']

    pars = {}
    for p in model_info['parameters']:
        name, default, limits = p[0], p[2], p[3]
        value = kwargs.pop(name, default)
        pars[name] = Parameter.default(value, name=name, limits=limits)
    for name in partype['pd-2d']:
        for xpart, xdefault, xlimits in [
            ('_pd', 0., limits),
            ('_pd_n', 35., (0, 1000)),
            ('_pd_nsigma', 3., (0, 10)),
        ]:
            xname = name + xpart
            xvalue = kwargs.pop(xname, xdefault)
            pars[xname] = Parameter.default(xvalue, name=xname, limits=xlimits)

    pd_types = {}
    for name in partype['pd-2d']:
        xname = name + '_pd_type'
        xvalue = kwargs.pop(xname, 'gaussian')
        pd_types[xname] = xvalue

    if kwargs:  # args not corresponding to parameters
        raise TypeError("unexpected parameters: %s"
                        % (", ".join(sorted(kwargs.keys()))))

    return pars, pd_types

class Model(object):
    def __init__(self, model, **kwargs):
        self._sasmodel = model
        pars, pd_types = create_parameters(model.info, **kwargs)
        for k,v in pars.items():
            setattr(self, k, v)
        for k,v in pd_types.items():
            setattr(self, k, v)
        self._parameter_names = list(pars.keys())
        self._pd_type_names = list(pd_types.keys())

    def parameters(self):
        """
        Return a dictionary of parameters
        """
        return dict((k, getattr(self, k)) for k in self._parameter_names)

    def state(self):
        pars = dict((k, getattr(self, k).value) for k in self._parameter_names)
        pars.update((k, getattr(self, k)) for k in self._pd_type_names)
        return pars

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
            pars = self.model.state()
            self._cache['theory'] = self._calc_theory(pars, cutoff=self.cutoff)
        return self._cache['theory']

    def residuals(self):
        #if np.any(self.err ==0): print("zeros in err")
        return (self.theory() - self.Iq) / self.dIq

    def nllf(self):
        delta = self.residuals()
        #if np.any(np.isnan(R)): print("NaN in residuals")
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

    def __getstate__(self):
        # Can't pickle gpu functions, so instead make them lazy
        state = self.__dict__.copy()
        state['_kernel'] = None
        return state

    def __setstate__(self, state):
        # pylint: disable=attribute-defined-outside-init
        self.__dict__ = state
