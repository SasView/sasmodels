"""
Sasview model constructor.

Given a module defining an OpenCL kernel such as sasmodels.models.cylinder,
create a sasview model class to run that kernel as follows::

    from sasmodels.sasview_model import load_custom_model
    CylinderModel = load_custom_model('sasmodels/models/cylinder.py')
"""
from __future__ import print_function

import math
from copy import deepcopy
import collections
import traceback
import logging
from os.path import basename, splitext

import numpy as np  # type: ignore

from . import core
from . import custom
from . import generate
from . import weights
from . import modelinfo
from .details import make_kernel_args, dispersion_mesh

try:
    from typing import Dict, Mapping, Any, Sequence, Tuple, NamedTuple, List, Optional, Union, Callable
    from .modelinfo import ModelInfo, Parameter
    from .kernel import KernelModel
    MultiplicityInfoType = NamedTuple(
        'MuliplicityInfo',
        [("number", int), ("control", str), ("choices", List[str]),
         ("x_axis_label", str)])
    SasviewModelType = Callable[[int], "SasviewModel"]
except ImportError:
    pass

SUPPORT_OLD_STYLE_PLUGINS = True

def _register_old_models():
    # type: () -> None
    """
    Place the new models into sasview under the old names.

    Monkey patch sas.sascalc.fit as sas.models so that sas.models.pluginmodel
    is available to the plugin modules.
    """
    import sys
    import sas
    import sas.sascalc.fit
    sys.modules['sas.models'] = sas.sascalc.fit
    sas.models = sas.sascalc.fit

    import sas.models
    from sasmodels.conversion_table import CONVERSION_TABLE
    for new_name, conversion in CONVERSION_TABLE.items():
        old_name = conversion[0]
        module_attrs = {old_name: find_model(new_name)}
        ConstructedModule = type(old_name, (), module_attrs)
        old_path = 'sas.models.' + old_name
        setattr(sas.models, old_path, ConstructedModule)
        sys.modules[old_path] = ConstructedModule


# TODO: separate x_axis_label from multiplicity info
MultiplicityInfo = collections.namedtuple(
    'MultiplicityInfo',
    ["number", "control", "choices", "x_axis_label"],
)

MODELS = {}
def find_model(modelname):
    # type: (str) -> SasviewModelType
    """
    Find a model by name.  If the model name ends in py, try loading it from
    custom models, otherwise look for it in the list of builtin models.
    """
    # TODO: used by sum/product model to load an existing model
    # TODO: doesn't handle custom models properly
    if modelname.endswith('.py'):
        return load_custom_model(modelname)
    elif modelname in MODELS:
        return MODELS[modelname]
    else:
        raise ValueError("unknown model %r"%modelname)


# TODO: figure out how to say that the return type is a subclass
def load_standard_models():
    # type: () -> List[SasviewModelType]
    """
    Load and return the list of predefined models.

    If there is an error loading a model, then a traceback is logged and the
    model is not returned.
    """
    models = []
    for name in core.list_models():
        try:
            MODELS[name] = _make_standard_model(name)
            models.append(MODELS[name])
        except Exception:
            logging.error(traceback.format_exc())
    if SUPPORT_OLD_STYLE_PLUGINS:
        _register_old_models()

    return models


def load_custom_model(path):
    # type: (str) -> SasviewModelType
    """
    Load a custom model given the model path.
    """
    kernel_module = custom.load_custom_kernel_module(path)
    try:
        model = kernel_module.Model
        # Old style models do not set the name in the class attributes, so
        # set it here; this name will be overridden when the object is created
        # with an instance variable that has the same value.
        if model.name == "":
            model.name = splitext(basename(path))[0]
        if not hasattr(model, 'filename'):
            model.filename = kernel_module.__file__
    except AttributeError:
        model_info = modelinfo.make_model_info(kernel_module)
        model = _make_model_from_info(model_info)

    # If a model name already exists and we are loading a different model,
    # use the model file name as the model name.
    if model.name in MODELS and not model.filename == MODELS[model.name].filename:
        _previous_name = model.name
        model.name = model.id
        
        # If the new model name is still in the model list (for instance,
        # if we put a cylinder.py in our plug-in directory), then append
        # an identifier.
        if model.name in MODELS and not model.filename == MODELS[model.name].filename:
            model.name = model.id + '_user'
        logging.info("Model %s already exists: using %s [%s]", _previous_name, model.name, model.filename)

    MODELS[model.name] = model
    return model


def _make_standard_model(name):
    # type: (str) -> SasviewModelType
    """
    Load the sasview model defined by *name*.

    *name* can be a standard model name or a path to a custom model.

    Returns a class that can be used directly as a sasview model.
    """
    kernel_module = generate.load_kernel_module(name)
    model_info = modelinfo.make_model_info(kernel_module)
    return _make_model_from_info(model_info)


def _make_model_from_info(model_info):
    # type: (ModelInfo) -> SasviewModelType
    """
    Convert *model_info* into a SasView model wrapper.
    """
    def __init__(self, multiplicity=None):
        SasviewModel.__init__(self, multiplicity=multiplicity)
    attrs = _generate_model_attributes(model_info)
    attrs['__init__'] = __init__
    attrs['filename'] = model_info.filename
    ConstructedModel = type(model_info.name, (SasviewModel,), attrs) # type: SasviewModelType
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
                len(p.choices) if p.choices else int(p.limits[1]),
                p.name, p.choices, xlabel
            )
            break

    # Only a single drop-down list parameter available
    fun_list = []
    for p in model_info.parameters.kernel_parameters:
        if p.choices:
            fun_list = p.choices
            if p.length > 1:
                non_fittable.extend(p.id+str(k) for k in range(1, p.length+1))
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
        elif p.type == 'magnetic':
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
    attrs['fun_list'] = tuple(fun_list)

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
    details = None     # type: Dict[str, Sequence[Any]]
    #                  # actual type is Dict[str, List[str, float, float]]
    #: multiplicity value, or None if no multiplicity on the model
    multiplicity = None     # type: Optional[int]
    #: memory for polydispersity array if using ArrayDispersion (used by sasview).
    _persistency_dict = None # type: Dict[str, Tuple[np.ndarray, np.ndarray]]

    def __init__(self, multiplicity=None):
        # type: (Optional[int]) -> None

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
        if self._model_info.structure_factor:
            hidden.add('scale')
            hidden.add('background')
            self._model_info.parameters.defaults['background'] = 0.

        self._persistency_dict = {}
        self.params = collections.OrderedDict()
        self.dispersion = collections.OrderedDict()
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
        return par_name in self.fixed
        #For the future
        #return self.params[str(par_name)].is_fittable()


    def getProfile(self):
        # type: () -> (np.ndarray, np.ndarray)
        """
        Get SLD profile

        : return: (z, beta) where z is a list of depth of the transition points
                beta is a list of the corresponding SLD values
        """
        args = {} # type: Dict[str, Any]
        for p in self._model_info.parameters.kernel_parameters:
            if p.id == self.multiplicity_info.control:
                value = float(self.multiplicity)
            elif p.length == 1:
                value = self.params.get(p.id, np.NaN)
            else:
                value = np.array([self.params.get(p.id+str(k), np.NaN)
                                  for k in range(1, p.length+1)])
            args[p.id] = value

        x, y = self._model_info.profile(**args)
        return x, 1e-6*y

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
                if item == toks[0]:
                    for par in self.dispersion[item]:
                        if par == toks[1]:
                            self.dispersion[item][par] = value
                            return
        else:
            # Look for standard parameter
            for item in self.params.keys():
                if item == name:
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
                if item == toks[0]:
                    for par in self.dispersion[item]:
                        if par == toks[1]:
                            return self.dispersion[item][par]
        else:
            # Look for standard parameter
            for item in self.params.keys():
                if item == name:
                    return self.params[item]

        raise ValueError("Model does not contain parameter %s" % name)

    def getParamList(self):
        # type: () -> Sequence[str]
        """
        Return a list of all available parameters for the model
        """
        param_list = list(self.params.keys())
        # WARNING: Extending the list with the dispersion parameters
        param_list.extend(self.getDispParamList())
        return param_list

    def getDispParamList(self):
        # type: () -> Sequence[str]
        """
        Return a list of polydispersity parameters for the model
        """
        # TODO: fix test so that parameter order doesn't matter
        ret = ['%s.%s' % (p_name, ext)
               for p_name in self.dispersion.keys()
               for ext in ('npts', 'nsigmas', 'width')]
        #print(ret)
        return ret

    def clone(self):
        # type: () -> "SasviewModel"
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
            return self.calculate_Iq([q*math.cos(phi)], [q*math.sin(phi)])[0]
        else:
            return self.calculate_Iq([x])[0]


    def runXY(self, x=0.0):
        # type: (Union[float, (float, float), List[float]]) -> float
        """
        Evaluate the model in cartesian coordinates

        :param x: input q, or [qx, qy]

        :return: scattering function P(q)

        **DEPRECATED**: use calculate_Iq instead
        """
        if isinstance(x, (list, tuple)):
            return self.calculate_Iq([x[0]], [x[1]])[0]
        else:
            return self.calculate_Iq([x])[0]

    def evalDistribution(self, qdist):
        # type: (Union[np.ndarray, Tuple[np.ndarray, np.ndarray], List[np.ndarray]]) -> np.ndarray
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
        #core.HAVE_OPENCL = False
        if self._model is None:
            self._model = core.build_model(self._model_info)
        if qy is not None:
            q_vectors = [np.asarray(qx), np.asarray(qy)]
        else:
            q_vectors = [np.asarray(qx)]
        calculator = self._model.make_kernel(q_vectors)
        parameters = self._model_info.parameters
        pairs = [self._get_weights(p) for p in parameters.call_parameters]
        #weights.plot_weights(self._model_info, pairs)
        call_details, values, is_magnetic = make_kernel_args(calculator, pairs)
        #call_details.show()
        #print("pairs", pairs)
        #print("params", self.params)
        #print("values", values)
        #print("is_mag", is_magnetic)
        result = calculator(call_details, values, cutoff=self.cutoff,
                            magnetic=is_magnetic)
        calculator.release()
        try:
            self._model.release()
        except:
            pass
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
        if parameter in self.params:
            # TODO: Store the disperser object directly in the model.
            # The current method of relying on the sasview GUI to
            # remember them is kind of funky.
            # Note: can't seem to get disperser parameters from sasview
            # (1) Could create a sasview model that has not yet been
            # converted, assign the disperser to one of its polydisperse
            # parameters, then retrieve the disperser parameters from the
            # sasview model.
            # (2) Could write a disperser parameter retriever in sasview.
            # (3) Could modify sasview to use sasmodels.weights dispersers.
            # For now, rely on the fact that the sasview only ever uses
            # new dispersers in the set_dispersion call and create a new
            # one instead of trying to assign parameters.
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
        return dispersion_mesh(self._model_info, pars)

    def _get_weights(self, par):
        # type: (Parameter) -> Tuple[np.ndarray, np.ndarray]
        """
        Return dispersion weights for parameter
        """
        if par.name not in self.params:
            if par.name == self.multiplicity_info.control:
                return [self.multiplicity], [1.0]
            else:
                # For hidden parameters use the default value.
                value = self._model_info.parameters.defaults.get(par.name, np.NaN)
                return [value], [1.0]
        elif par.polydisperse:
            dis = self.dispersion[par.name]
            if dis['type'] == 'array':
                value, weight = dis['values'], dis['weights']
            else:
                value, weight = weights.get_weights(
                    dis['type'], dis['npts'], dis['width'], dis['nsigmas'],
                    self.params[par.name], par.limits, par.relative_pd)
            return value, weight / np.sum(weight)
        else:
            return [self.params[par.name]], [1.0]

def test_model():
    # type: () -> float
    """
    Test that a sasview model (cylinder) can be run.
    """
    Cylinder = _make_standard_model('cylinder')
    cylinder = Cylinder()
    return cylinder.evalDistribution([0.1, 0.1])

def test_structure_factor():
    # type: () -> float
    """
    Test that a sasview model (cylinder) can be run.
    """
    Model = _make_standard_model('hardsphere')
    model = Model()
    value = model.evalDistribution([0.1, 0.1])
    if np.isnan(value):
        raise ValueError("hardsphere returns null")

def test_rpa():
    # type: () -> float
    """
    Test that a sasview model (cylinder) can be run.
    """
    RPA = _make_standard_model('rpa')
    rpa = RPA(3)
    return rpa.evalDistribution([0.1, 0.1])


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

def test_old_name():
    # type: () -> None
    """
    Load and run cylinder model from sas.models.CylinderModel
    """
    if not SUPPORT_OLD_STYLE_PLUGINS:
        return
    try:
        # if sasview is not on the path then don't try to test it
        import sas
    except ImportError:
        return
    load_standard_models()
    from sas.models.CylinderModel import CylinderModel
    CylinderModel().evalDistribution([0.1, 0.1])

if __name__ == "__main__":
    print("cylinder(0.1,0.1)=%g"%test_model())
