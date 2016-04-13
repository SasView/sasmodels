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

import numpy as np  # type: ignore

from . import core
from . import custom
from . import generate
from . import weights
from . import details
from . import modelinfo

try:
    from typing import Dict, Mapping, Any, Sequence, Tuple, NamedTuple, List, Optional
    from .modelinfo import ModelInfo, Parameter
    from .kernel import KernelModel
    MultiplicityInfoType = NamedTuple(
        'MuliplicityInfo',
        [("number", int), ("control", str), ("choices", List[str]),
         ("x_axis_label", str)])
except ImportError:
    pass

# TODO: separate x_axis_label from multiplicity info
# The x-axis label belongs with the profile generating function
MultiplicityInfo = collections.namedtuple(
    'MultiplicityInfo',
    ["number", "control", "choices", "x_axis_label"],
)

# TODO: figure out how to say that the return type is a subclass
def load_standard_models():
    # type: () -> List[type]
    """
    Load and return the list of predefined models.

    If there is an error loading a model, then a traceback is logged and the
    model is not returned.
    """
    models = []
    for name in core.list_models():
        try:
            models.append(_make_standard_model(name))
        except Exception:
            logging.error(traceback.format_exc())
    return models


def load_custom_model(path):
    # type: (str) -> type
    """
    Load a custom model given the model path.
    """
    kernel_module = custom.load_custom_kernel_module(path)
    model_info = modelinfo.make_model_info(kernel_module)
    return _make_model_from_info(model_info)


def _make_standard_model(name):
    # type: (str) -> type
    """
    Load the sasview model defined by *name*.

    *name* can be a standard model name or a path to a custom model.

    Returns a class that can be used directly as a sasview model.
    """
    kernel_module = generate.load_kernel_module(name)
    model_info = modelinfo.make_model_info(kernel_module)
    return _make_model_from_info(model_info)


def _make_model_from_info(model_info):
    # type: (ModelInfo) -> type
    """
    Convert *model_info* into a SasView model wrapper.
    """
    def __init__(self, multiplicity=None):
        SasviewModel.__init__(self, multiplicity=multiplicity)
    attrs = _generate_model_attributes(model_info)
    attrs['__init__'] = __init__
    ConstructedModel = type(model_info.name, (SasviewModel,), attrs)
    return ConstructedModel

def _generate_model_attributes(model_info):
    # type: (ModelInfo) -> Dict[str, Any]
    """
    Generate the class attributes for the model.

    This should include all the information necessary to query the model
    details so that you do not need to instantiate a model to query it.

    All the attributes should be immutable to avoid accidents.
    """

    # TODO: allow model to override axis labels input/output name/unit

    # Process multiplicity
    non_fittable = []  # type: List[str]
    xlabel = model_info.profile_axes[0] if model_info.profile is not None else ""
    variants = MultiplicityInfo(0, "", [], xlabel)
    for p in model_info.parameters.kernel_parameters:
        if p.name == model_info.control:
            non_fittable.append(p.name)
            variants = MultiplicityInfo(
                len(p.choices), p.name, p.choices, xlabel
            )
            break
        elif p.is_control:
            non_fittable.append(p.name)
            variants = MultiplicityInfo(
                int(p.limits[1]), p.name, p.choices, xlabel
            )
            break

    # Organize parameter sets
    orientation_params = []
    magnetic_params = []
    fixed = []
    for p in model_info.parameters.user_parameters():
        if p.type == 'orientation':
            orientation_params.append(p.name)
            orientation_params.append(p.name+".width")
            fixed.append(p.name+".width")
        if p.type == 'magnetic':
            orientation_params.append(p.name)
            magnetic_params.append(p.name)
            fixed.append(p.name+".width")

    # Build class dictionary
    attrs = {}  # type: Dict[str, Any]
    attrs['_model_info'] = model_info
    attrs['name'] = model_info.name
    attrs['id'] = model_info.id
    attrs['description'] = model_info.description
    attrs['category'] = model_info.category
    attrs['is_structure_factor'] = model_info.structure_factor
    attrs['is_form_factor'] = model_info.ER is not None
    attrs['is_multiplicity_model'] = variants[0] > 1
    attrs['multiplicity_info'] = variants
    attrs['orientation_params'] = tuple(orientation_params)
    attrs['magnetic_params'] = tuple(magnetic_params)
    attrs['fixed'] = tuple(fixed)
    attrs['non_fittable'] = tuple(non_fittable)

    return attrs

class SasviewModel(object):
    """
    Sasview wrapper for opencl/ctypes model.
    """
    # Model parameters for the specific model are set in the class constructor
    # via the _generate_model_attributes function, which subclasses
    # SasviewModel.  They are included here for typing and documentation
    # purposes.
    _model = None       # type: KernelModel
    _model_info = None  # type: ModelInfo
    #: load/save name for the model
    id = None           # type: str
    #: display name for the model
    name = None         # type: str
    #: short model description
    description = None  # type: str
    #: default model category
    category = None     # type: str

    #: names of the orientation parameters in the order they appear
    orientation_params = None # type: Sequence[str]
    #: names of the magnetic parameters in the order they appear
    magnetic_params = None    # type: Sequence[str]
    #: names of the fittable parameters
    fixed = None              # type: Sequence[str]
    # TODO: the attribute fixed is ill-named

    # Axis labels
    input_name = "Q"
    input_unit = "A^{-1}"
    output_name = "Intensity"
    output_unit = "cm^{-1}"

    #: default cutoff for polydispersity
    cutoff = 1e-5

    # Note: Use non-mutable values for class attributes to avoid errors
    #: parameters that are not fitted
    non_fittable = ()        # type: Sequence[str]

    #: True if model should appear as a structure factor
    is_structure_factor = False
    #: True if model should appear as a form factor
    is_form_factor = False
    #: True if model has multiplicity
    is_multiplicity_model = False
    #: Mulitplicity information
    multiplicity_info = None # type: MultiplicityInfoType

    # Per-instance variables
    #: parameter {name: value} mapping
    params = None      # type: Dict[str, float]
    #: values for dispersion width, npts, nsigmas and type
    dispersion = None  # type: Dict[str, Any]
    #: units and limits for each parameter
    details = None     # type: Mapping[str, Tuple(str, float, float)]
    #: multiplicity used, or None if no multiplicity controls
    multiplicity = None     # type: Optional[int]

    def __init__(self, multiplicity=None):
        # type: () -> None

        ## _persistency_dict is used by sas.perspectives.fitting.basepage
        ## to store dispersity reference.
        self._persistency_dict = {}

        # TODO: _persistency_dict to persistency_dict throughout sasview
        # TODO: refactor multiplicity to encompass variants
        # TODO: dispersion should be a class
        # TODO: refactor multiplicity info
        # TODO: separate profile view from multiplicity
        # The button label, x and y axis labels and scale need to be under
        # the control of the model, not the fit page.  Maximum flexibility,
        # the fit page would supply the canvas and the profile could plot
        # how it wants, but this assumes matplotlib.  Next level is that
        # we provide some sort of data description including title, labels
        # and lines to plot.

        # Get the list of hidden parameters given the mulitplicity
        # Don't include multiplicity in the list of parameters
        self.multiplicity = multiplicity
        if multiplicity is not None:
            hidden = self._model_info.get_hidden_parameters(multiplicity)
            hidden |= set([self.multiplicity_info.control])
        else:
            hidden = set()

        self.params = collections.OrderedDict()
        self.dispersion = {}
        self.details = {}
        for p in self._model_info.parameters.user_parameters():
            if p.name in hidden:
                continue
            self.params[p.name] = p.default
            self.details[p.id] = [p.units, p.limits[0], p.limits[1]]
            if p.polydisperse:
                self.details[p.id+".width"] = [
                    "", 0.0, 1.0 if p.relative_pd else np.inf
                ]
                self.dispersion[p.name] = {
                    'width': 0,
                    'npts': 35,
                    'nsigmas': 3,
                    'type': 'gaussian',
                }

    def __get_state__(self):
        # type: () -> Dict[str, Any]
        state = self.__dict__.copy()
        state.pop('_model')
        # May need to reload model info on set state since it has pointers
        # to python implementations of Iq, etc.
        #state.pop('_model_info')
        return state

    def __set_state__(self, state):
        # type: (Dict[str, Any]) -> None
        self.__dict__ = state
        self._model = None

    def __str__(self):
        # type: () -> str
        """
        :return: string representation
        """
        return self.name

    def is_fittable(self, par_name):
        # type: (str) -> bool
        """
        Check if a given parameter is fittable or not

        :param par_name: the parameter name to check
        """
        return par_name.lower() in self.fixed
        #For the future
        #return self.params[str(par_name)].is_fittable()


    def getProfile(self):
        # type: () -> (np.ndarray, np.ndarray)
        """
        Get SLD profile

        : return: (z, beta) where z is a list of depth of the transition points
                beta is a list of the corresponding SLD values
        """
        args = []
        for p in self._model_info.parameters.kernel_parameters:
            if p.id == self.multiplicity_info.control:
                args.append(self.multiplicity)
            elif p.length == 1:
                args.append(self.params.get(p.id, np.NaN))
            else:
                args.append([self.params.get(p.id+str(k), np.NaN)
                             for k in range(1,p.length+1)])
        return self._model_info.profile(*args)

    def setParam(self, name, value):
        # type: (str, float) -> None
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
        # type: (str) -> float
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
        # type: () - > Sequence[str]
        """
        Return a list of all available parameters for the model
        """
        param_list = self.params.keys()
        # WARNING: Extending the list with the dispersion parameters
        param_list.extend(self.getDispParamList())
        return param_list

    def getDispParamList(self):
        # type: () - > Sequence[str]
        """
        Return a list of polydispersity parameters for the model
        """
        # TODO: fix test so that parameter order doesn't matter
        ret = ['%s.%s' % (p.name.lower(), ext)
               for p in self._model_info.parameters.user_parameters()
               for ext in ('npts', 'nsigmas', 'width')
               if p.polydisperse]
        #print(ret)
        return ret

    def clone(self):
        # type: () - > "SasviewModel"
        """ Return a identical copy of self """
        return deepcopy(self)

    def run(self, x=0.0):
        # type: (Union[float, (float, float), List[float]]) -> float
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
        # type: (Union[float, (float, float), List[float]]) -> float
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
        # type: (Union[np.ndarray, Tuple[np.ndarray, np.ndarray], List[np.ndarray]) -> np.ndarray
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
            if not self._model_info.parameters.has_2d:
                return self.calculate_Iq(np.sqrt(qx ** 2 + qy ** 2))
            else:
                return self.calculate_Iq(qx, qy)

        elif isinstance(qdist, np.ndarray):
            # We have a simple 1D distribution of q-values
            return self.calculate_Iq(qdist)

        else:
            raise TypeError("evalDistribution expects q or [qx, qy], not %r"
                            % type(qdist))

    def calculate_Iq(self, qx, qy=None):
        # type: (Sequence[float], Optional[Sequence[float]]) -> np.ndarray
        """
        Calculate Iq for one set of q with the current parameters.

        If the model is 1D, use *q*.  If 2D, use *qx*, *qy*.

        This should NOT be used for fitting since it copies the *q* vectors
        to the card for each evaluation.
        """
        if self._model is None:
            self._model = core.build_model(self._model_info)
        if qy is not None:
            q_vectors = [np.asarray(qx), np.asarray(qy)]
        else:
            q_vectors = [np.asarray(qx)]
        kernel = self._model.make_kernel(q_vectors)
        pairs = [self._get_weights(p)
                 for p in self._model_info.parameters.call_parameters]
        call_details, weight, value = details.build_details(kernel, pairs)
        result = kernel(call_details, weight, value, cutoff=self.cutoff)
        kernel.release()
        return result

    def calculate_ER(self):
        # type: () -> float
        """
        Calculate the effective radius for P(q)*S(q)

        :return: the value of the effective radius
        """
        if self._model_info.ER is None:
            return 1.0
        else:
            value, weight = self._dispersion_mesh()
            fv = self._model_info.ER(*value)
            #print(values[0].shape, weights.shape, fv.shape)
            return np.sum(weight * fv) / np.sum(weight)

    def calculate_VR(self):
        # type: () -> float
        """
        Calculate the volf ratio for P(q)*S(q)

        :return: the value of the volf ratio
        """
        if self._model_info.VR is None:
            return 1.0
        else:
            value, weight = self._dispersion_mesh()
            whole, part = self._model_info.VR(*value)
            return np.sum(weight * part) / np.sum(weight * whole)

    def set_dispersion(self, parameter, dispersion):
        # type: (str, weights.Dispersion) -> Dict[str, Any]
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
            dispersion = weights.MODELS[disperser]()
            self.dispersion[parameter] = dispersion.get_pars()
        else:
            raise ValueError("%r is not a dispersity or orientation parameter")

    def _dispersion_mesh(self):
        # type: () -> List[Tuple[np.ndarray, np.ndarray]]
        """
        Create a mesh grid of dispersion parameters and weights.

        Returns [p1,p2,...],w where pj is a vector of values for parameter j
        and w is a vector containing the products for weights for each
        parameter set in the vector.
        """
        pars = [self._get_weights(p)
                for p in self._model_info.parameters.call_parameters
                if p.type == 'volume']
        return details.dispersion_mesh(self._model_info, pars)

    def _get_weights(self, par):
        # type: (Parameter) -> Tuple[np.ndarray, np.ndarray]
        """
        Return dispersion weights for parameter
        """
        if par.name not in self.params:
            if par.name == self.multiplicity_info.control:
                return [self.multiplicity], []
            else:
                return [np.NaN], []
        elif par.polydisperse:
            dis = self.dispersion[par.name]
            value, weight = weights.get_weights(
                dis['type'], dis['npts'], dis['width'], dis['nsigmas'],
                self.params[par.name], par.limits, par.relative_pd)
            return value, weight / np.sum(weight)
        else:
            return [self.params[par.name]], []

def test_model():
    # type: () -> float
    """
    Test that a sasview model (cylinder) can be run.
    """
    Cylinder = _make_standard_model('cylinder')
    cylinder = Cylinder()
    return cylinder.evalDistribution([0.1,0.1])

def test_rpa():
    # type: () -> float
    """
    Test that a sasview model (cylinder) can be run.
    """
    RPA = _make_standard_model('rpa')
    rpa = RPA(3)
    return rpa.evalDistribution([0.1,0.1])


def test_model_list():
    # type: () -> None
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
