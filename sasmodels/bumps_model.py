"""
Wrap sasmodels for direct use by bumps.

:class:`Model` is a wrapper for the sasmodels kernel which defines a
bumps *Parameter* box for each kernel parameter.  *Model* accepts keyword
arguments to set the initial value for each parameter.

:class:`Experiment` combines the *Model* function with a data file loaded by
the sasview data loader.  *Experiment* takes a *cutoff* parameter controlling
how far the polydispersity integral extends.

"""

import warnings

import numpy as np

from .data import plot_theory
from .direct_model import DataMixin

__all__ = [
    "Model", "Experiment",
    ]

# CRUFT: old style bumps wrapper which doesn't separate data and model
# pylint: disable=invalid-name
def BumpsModel(data, model, cutoff=1e-5, **kw):
    r"""
    Bind a model to data, along with a polydispersity cutoff.

    *data* is a :class:`data.Data1D`, :class:`data.Data2D` or
    :class:`data.Sesans` object.  Use :func:`data.empty_data1D` or
    :func:`data.empty_data2D` to define $q, \Delta q$ calculation
    points for displaying the SANS curve when there is no measured data.

    *model* is a runnable module as returned from :func:`core.load_model`.

    *cutoff* is the polydispersity weight cutoff.

    Any additional *key=value* pairs are model dependent parameters.

    Returns an :class:`Experiment` object.

    Note that the usual Bumps semantics is not fully supported, since
    assigning *M.name = parameter* on the returned experiment object
    does not set that parameter in the model.  Range setting will still
    work as expected though.

    .. deprecated:: 0.1
        Use :class:`Experiment` instead.
    """
    warnings.warn("Use of BumpsModel is deprecated.  Use bumps_model.Experiment instead.")

    # Create the model and experiment
    model = Model(model, **kw)
    experiment = Experiment(data=data, model=model, cutoff=cutoff)

    # Copy the model parameters up to the experiment object.
    for k, v in model.parameters().items():
        setattr(experiment, k, v)
    return experiment


def create_parameters(model_info, **kwargs):
    """
    Generate Bumps parameters from the model info.

    *model_info* is returned from :func:`generate.model_info` on the
    model definition module.

    Any additional *key=value* pairs are initial values for the parameters
    to the models.  Uninitialized parameters will use the model default
    value.

    Returns a dictionary of *{name: Parameter}* containing the bumps
    parameters for each model parameter, and a dictionary of
    *{name: str}* containing the polydispersity distribution types.
    """
    # lazy import; this allows the doc builder and nosetests to run even
    # when bumps is not on the path.
    from bumps.names import Parameter

    pars = {}
    for p in model_info['parameters']:
        value = kwargs.pop([p.name, p.default])
        pars[p.name] = Parameter.default(value, name=p.name, limits=p.limits)
    for name in model_info['partype']['pd-2d']:
        for xpart, xdefault, xlimits in [
                ('_pd', 0., pars[name].limits),
                ('_pd_n', 35., (0, 1000)),
                ('_pd_nsigma', 3., (0, 10)),
            ]:
            xname = name + xpart
            xvalue = kwargs.pop(xname, xdefault)
            pars[xname] = Parameter.default(xvalue, name=xname, limits=xlimits)

    pd_types = {}
    for name in model_info['partype']['pd-2d']:
        xname = name + '_pd_type'
        xvalue = kwargs.pop(xname, 'gaussian')
        pd_types[xname] = xvalue

    if kwargs:  # args not corresponding to parameters
        raise TypeError("unexpected parameters: %s"
                        % (", ".join(sorted(kwargs.keys()))))

    return pars, pd_types

class Model(object):
    """
    Bumps wrapper for a SAS model.

    *model* is a runnable module as returned from :func:`core.load_model`.

    *cutoff* is the polydispersity weight cutoff.

    Any additional *key=value* pairs are model dependent parameters.
    """
    def __init__(self, model, **kwargs):
        self._sasmodel = model
        pars, pd_types = create_parameters(model.info, **kwargs)
        for k, v in pars.items():
            setattr(self, k, v)
        for k, v in pd_types.items():
            setattr(self, k, v)
        self._parameter_names = list(pars.keys())
        self._pd_type_names = list(pd_types.keys())

    def parameters(self):
        """
        Return a dictionary of parameters objects for the parameters,
        excluding polydispersity distribution type.
        """
        return dict((k, getattr(self, k)) for k in self._parameter_names)

    def state(self):
        """
        Return a dictionary of current values for all the parameters,
        including polydispersity distribution type.
        """
        pars = dict((k, getattr(self, k).value) for k in self._parameter_names)
        pars.update((k, getattr(self, k)) for k in self._pd_type_names)
        return pars

class Experiment(DataMixin):
    r"""
    Bumps wrapper for a SAS experiment.

    *data* is a :class:`data.Data1D`, :class:`data.Data2D` or
    :class:`data.Sesans` object.  Use :func:`data.empty_data1D` or
    :func:`data.empty_data2D` to define $q, \Delta q$ calculation
    points for displaying the SANS curve when there is no measured data.

    *model* is a :class:`Model` object.

    *cutoff* is the integration cutoff, which avoids computing the
    the SAS model where the polydispersity weight is low.

    The resulting model can be used directly in a Bumps FitProblem call.
    """
    def __init__(self, data, model, cutoff=1e-5):

        # remember inputs so we can inspect from outside
        self.model = model
        self.cutoff = cutoff
        self._interpret_data(data, model._sasmodel)
        self.update()

    def update(self):
        """
        Call when model parameters have changed and theory needs to be
        recalculated.
        """
        self._cache = {}

    def numpoints(self):
        """
        Return the number of data points
        """
        return len(self.Iq)

    def parameters(self):
        """
        Return a dictionary of parameters
        """
        return self.model.parameters()

    def theory(self):
        """
        Return the theory corresponding to the model parameters.

        This method uses lazy evaluation, and requires model.update() to be
        called when the parameters have changed.
        """
        if 'theory' not in self._cache:
            pars = self.model.state()
            self._cache['theory'] = self._calc_theory(pars, cutoff=self.cutoff)
        return self._cache['theory']

    def residuals(self):
        """
        Return theory minus data normalized by uncertainty.
        """
        #if np.any(self.err ==0): print("zeros in err")
        return (self.theory() - self.Iq) / self.dIq

    def nllf(self):
        """
        Return the negative log likelihood of seeing data given the model
        parameters, up to a normalizing constant which depends on the data
        uncertainty.
        """
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
        """
        Generate simulated data.
        """
        Iq = self.theory()
        self._set_data(Iq, noise)

    def save(self, basename):
        """
        Save the model parameters and data into a file.

        Not Implemented.
        """
        pass

    def __getstate__(self):
        # Can't pickle gpu functions, so instead make them lazy
        state = self.__dict__.copy()
        state['_kernel'] = None
        return state

    def __setstate__(self, state):
        # pylint: disable=attribute-defined-outside-init
        self.__dict__ = state
