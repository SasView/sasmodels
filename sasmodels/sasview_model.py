"""
Sasview model constructor.

Given a module defining an OpenCL kernel such as sasmodels.models.cylinder,
create a sasview model class to run that kernel as follows::

    from sasmodels.sasview_model import make_class
    from sasmodels.models import cylinder
    CylinderModel = make_class(cylinder, dtype='single')

The model parameters for sasmodels are different from those in sasview.
When reloading previously saved models, the parameters should be converted
using :func:`sasmodels.convert.convert`.
"""
from __future__ import print_function

import math
from copy import deepcopy
import collections
import traceback
import logging

import numpy as np

from . import core
from . import custom
from . import generate
from . import weights

def load_standard_models():
    """
    Load and return the list of predefined models.

    If there is an error loading a model, then a traceback is logged and the
    model is not returned.
    """
    models = []
    for name in core.list_models():
        try:
            models.append(_make_standard_model(name))
        except:
            logging.error(traceback.format_exc())
    return models


def load_custom_model(path):
    """
    Load a custom model given the model path.
    """
    kernel_module = custom.load_custom_kernel_module(path)
    model_info = generate.make_model_info(kernel_module)
    return _make_model_from_info(model_info)


def _make_standard_model(name):
    """
    Load the sasview model defined by *name*.

    *name* can be a standard model name or a path to a custom model.

    Returns a class that can be used directly as a sasview model.
    """
    kernel_module = generate.load_kernel_module(name)
    model_info = generate.make_model_info(kernel_module)
    return _make_model_from_info(model_info)


def _make_model_from_info(model_info):
    """
    Convert *model_info* into a SasView model wrapper.
    """
    def __init__(self, multfactor=1):
        SasviewModel.__init__(self)
    attrs = dict(__init__=__init__, _model_info=model_info)
    ConstructedModel = type(model_info['name'], (SasviewModel,), attrs)
    return ConstructedModel


class SasviewModel(object):
    """
    Sasview wrapper for opencl/ctypes model.
    """
    _model_info = {}
    def __init__(self):
        self._model = None
        model_info = self._model_info
        parameters = model_info['parameters']

        self.name = model_info['name']
        self.description = model_info['description']
        self.category = None
        #self.is_multifunc = False
        for p in parameters.kernel_parameters:
            if p.is_control:
                profile_axes = model_info['profile_axes']
                self.multiplicity_info = [
                    p.limits[1], p.name, p.choices, profile_axes[0]
                    ]
                break
        else:
            self.multiplicity_info = []

        ## interpret the parameters
        ## TODO: reorganize parameter handling
        self.details = dict()
        self.params = collections.OrderedDict()
        self.dispersion = dict()

        self.orientation_params = []
        self.magnetic_params = []
        self.fixed = []
        for p in parameters.user_parameters():
            self.params[p.name] = p.default
            self.details[p.name] = [p.units] + p.limits
            if p.polydisperse:
                self.dispersion[p.name] = {
                    'width': 0,
                    'npts': 35,
                    'nsigmas': 3,
                    'type': 'gaussian',
                }
            if p.type == 'orientation':
                self.orientation_params.append(p.name)
                self.orientation_params.append(p.name+".width")
                self.fixed.append(p.name+".width")
            if p.type == 'magnetic':
                self.orientation_params.append(p.name)
                self.magnetic_params.append(p.name)
                self.fixed.append(p.name+".width")

        self.non_fittable = []

        ## independent parameter name and unit [string]
        self.input_name = model_info.get("input_name", "Q")
        self.input_unit = model_info.get("input_unit", "A^{-1}")
        self.output_name = model_info.get("output_name", "Intensity")
        self.output_unit = model_info.get("output_unit", "cm^{-1}")

        ## _persistency_dict is used by sas.perspectives.fitting.basepage
        ## to store dispersity reference.
        ## TODO: _persistency_dict to persistency_dict throughout sasview
        self._persistency_dict = {}

        ## New fields introduced for opencl rewrite
        self.cutoff = 1e-5

    def __get_state__(self):
        state = self.__dict__.copy()
        state.pop('_model')
        # May need to reload model info on set state since it has pointers
        # to python implementations of Iq, etc.
        #state.pop('_model_info')
        return state

    def __set_state__(self, state):
        self.__dict__ = state
        self._model = None

    def __str__(self):
        """
        :return: string representation
        """
        return self.name

    def is_fittable(self, par_name):
        """
        Check if a given parameter is fittable or not

        :param par_name: the parameter name to check
        """
        return par_name.lower() in self.fixed
        #For the future
        #return self.params[str(par_name)].is_fittable()


    # pylint: disable=no-self-use
    def getProfile(self):
        """
        Get SLD profile

        : return: (z, beta) where z is a list of depth of the transition points
                beta is a list of the corresponding SLD values
        """
        return None, None

    def setParam(self, name, value):
        """
        Set the value of a model parameter

        :param name: name of the parameter
        :param value: value of the parameter

        """
        # Look for dispersion parameters
        toks = name.split('.')
        if len(toks) == 2:
            for item in self.dispersion.keys():
                if item.lower() == toks[0].lower():
                    for par in self.dispersion[item]:
                        if par.lower() == toks[1].lower():
                            self.dispersion[item][par] = value
                            return
        else:
            # Look for standard parameter
            for item in self.params.keys():
                if item.lower() == name.lower():
                    self.params[item] = value
                    return

        raise ValueError("Model does not contain parameter %s" % name)

    def getParam(self, name):
        """
        Set the value of a model parameter

        :param name: name of the parameter

        """
        # Look for dispersion parameters
        toks = name.split('.')
        if len(toks) == 2:
            for item in self.dispersion.keys():
                if item.lower() == toks[0].lower():
                    for par in self.dispersion[item]:
                        if par.lower() == toks[1].lower():
                            return self.dispersion[item][par]
        else:
            # Look for standard parameter
            for item in self.params.keys():
                if item.lower() == name.lower():
                    return self.params[item]

        raise ValueError("Model does not contain parameter %s" % name)

    def getParamList(self):
        """
        Return a list of all available parameters for the model
        """
        param_list = self.params.keys()
        # WARNING: Extending the list with the dispersion parameters
        param_list.extend(self.getDispParamList())
        return param_list

    def getDispParamList(self):
        """
        Return a list of polydispersity parameters for the model
        """
        # TODO: fix test so that parameter order doesn't matter
        ret = ['%s.%s' % (p.name.lower(), ext)
               for p in self._model_info['parameters'].user_parameters()
               for ext in ('npts', 'nsigmas', 'width')
               if p.polydisperse]
        #print(ret)
        return ret

    def clone(self):
        """ Return a identical copy of self """
        return deepcopy(self)

    def run(self, x=0.0):
        """
        Evaluate the model

        :param x: input q, or [q,phi]

        :return: scattering function P(q)

        **DEPRECATED**: use calculate_Iq instead
        """
        if isinstance(x, (list, tuple)):
            # pylint: disable=unpacking-non-sequence
            q, phi = x
            return self.calculate_Iq([q * math.cos(phi)],
                                     [q * math.sin(phi)])[0]
        else:
            return self.calculate_Iq([float(x)])[0]


    def runXY(self, x=0.0):
        """
        Evaluate the model in cartesian coordinates

        :param x: input q, or [qx, qy]

        :return: scattering function P(q)

        **DEPRECATED**: use calculate_Iq instead
        """
        if isinstance(x, (list, tuple)):
            return self.calculate_Iq([float(x[0])], [float(x[1])])[0]
        else:
            return self.calculate_Iq([float(x)])[0]

    def evalDistribution(self, qdist):
        r"""
        Evaluate a distribution of q-values.

        :param qdist: array of q or a list of arrays [qx,qy]

        * For 1D, a numpy array is expected as input

        ::

            evalDistribution(q)

          where *q* is a numpy array.

        * For 2D, a list of *[qx,qy]* is expected with 1D arrays as input

        ::

              qx = [ qx[0], qx[1], qx[2], ....]
              qy = [ qy[0], qy[1], qy[2], ....]

        If the model is 1D only, then

        .. math::

            q = \sqrt{q_x^2+q_y^2}

        """
        if isinstance(qdist, (list, tuple)):
            # Check whether we have a list of ndarrays [qx,qy]
            qx, qy = qdist
            if not self._model_info['parameters'].has_2d:
                return self.calculate_Iq(np.sqrt(qx ** 2 + qy ** 2))
            else:
                return self.calculate_Iq(qx, qy)

        elif isinstance(qdist, np.ndarray):
            # We have a simple 1D distribution of q-values
            return self.calculate_Iq(qdist)

        else:
            raise TypeError("evalDistribution expects q or [qx, qy], not %r"
                            % type(qdist))

    def calculate_Iq(self, *args):
        """
        Calculate Iq for one set of q with the current parameters.

        If the model is 1D, use *q*.  If 2D, use *qx*, *qy*.

        This should NOT be used for fitting since it copies the *q* vectors
        to the card for each evaluation.
        """
        if self._model is None:
            self._model = core.build_model(self._model_info)
        q_vectors = [np.asarray(q) for q in args]
        kernel = self._model.make_kernel(q_vectors)
        pairs = [self._get_weights(p)
                 for p in self._model_info['parameters'].call_parameters]
        details, weights, values = core.build_details(kernel, pairs)
        result = kernel(details, weights, values, cutoff=self.cutoff)
        kernel.q_input.release()
        kernel.release()
        return result

    def calculate_ER(self):
        """
        Calculate the effective radius for P(q)*S(q)

        :return: the value of the effective radius
        """
        ER = self._model_info.get('ER', None)
        if ER is None:
            return 1.0
        else:
            values, weights = self._dispersion_mesh()
            fv = ER(*values)
            #print(values[0].shape, weights.shape, fv.shape)
            return np.sum(weights * fv) / np.sum(weights)

    def calculate_VR(self):
        """
        Calculate the volf ratio for P(q)*S(q)

        :return: the value of the volf ratio
        """
        VR = self._model_info.get('VR', None)
        if VR is None:
            return 1.0
        else:
            values, weights = self._dispersion_mesh()
            whole, part = VR(*values)
            return np.sum(weights * part) / np.sum(weights * whole)

    def set_dispersion(self, parameter, dispersion):
        """
        Set the dispersion object for a model parameter

        :param parameter: name of the parameter [string]
        :param dispersion: dispersion object of type Dispersion
        """
        if parameter.lower() in (s.lower() for s in self.params.keys()):
            # TODO: Store the disperser object directly in the model.
            # The current method of creating one on the fly whenever it is
            # needed is kind of funky.
            # Note: can't seem to get disperser parameters from sasview
            # (1) Could create a sasview model that has not yet # been
            # converted, assign the disperser to one of its polydisperse
            # parameters, then retrieve the disperser parameters from the
            # sasview model.  (2) Could write a disperser parameter retriever
            # in sasview.  (3) Could modify sasview to use sasmodels.weights
            # dispersers.
            # For now, rely on the fact that the sasview only ever uses
            # new dispersers in the set_dispersion call and create a new
            # one instead of trying to assign parameters.
            from . import weights
            disperser = weights.dispersers[dispersion.__class__.__name__]
            dispersion = weights.models[disperser]()
            self.dispersion[parameter] = dispersion.get_pars()
        else:
            raise ValueError("%r is not a dispersity or orientation parameter")

    def _dispersion_mesh(self):
        """
        Create a mesh grid of dispersion parameters and weights.

        Returns [p1,p2,...],w where pj is a vector of values for parameter j
        and w is a vector containing the products for weights for each
        parameter set in the vector.
        """
        pars = self._model_info['partype']['volume']
        return core.dispersion_mesh([self._get_weights(p) for p in pars])

    def _get_weights(self, par):
        """
        Return dispersion weights for parameter
        """
        if par.polydisperse:
            dis = self.dispersion[par.name]
            value, weight = weights.get_weights(
                dis['type'], dis['npts'], dis['width'], dis['nsigmas'],
                self.params[par.name], par.limits, par.relative_pd)
            return value, weight / np.sum(weight)
        else:
            return [self.params[par.name]], []

def test_model():
    """
    Test that a sasview model (cylinder) can be run.
    """
    Cylinder = _make_standard_model('cylinder')
    cylinder = Cylinder()
    return cylinder.evalDistribution([0.1,0.1])


def test_model_list():
    """
    Make sure that all models build as sasview models.
    """
    from .exception import annotate_exception
    for name in core.list_models():
        try:
            _make_standard_model(name)
        except:
            annotate_exception("when loading "+name)
            raise

if __name__ == "__main__":
    print("cylinder(0.1,0.1)=%g"%test_model())
