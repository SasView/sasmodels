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
from os.path import basename, splitext, abspath, getmtime
try:
    import _thread as thread
except ImportError:
    import thread

import numpy as np  # type: ignore

from . import core
from . import custom
from . import kernelcl
from . import product
from . import generate
from . import weights
from . import modelinfo
from .details import make_kernel_args, dispersion_mesh
from .kernelcl import reset_environment

# pylint: disable=unused-import
try:
    from typing import (Dict, Mapping, Any, Sequence, Tuple, NamedTuple,
                        List, Optional, Union, Callable)
    from .modelinfo import ModelInfo, Parameter
    from .kernel import KernelModel
    MultiplicityInfoType = NamedTuple(
        'MultiplicityInfo',
        [("number", int), ("control", str), ("choices", List[str]),
         ("x_axis_label", str)])
    SasviewModelType = Callable[[int], "SasviewModel"]
except ImportError:
    pass
# pylint: enable=unused-import

logger = logging.getLogger(__name__)

calculation_lock = thread.allocate_lock()

#: True if pre-existing plugins, with the old names and parameters, should
#: continue to be supported.
SUPPORT_OLD_STYLE_PLUGINS = True

# TODO: separate x_axis_label from multiplicity info
MultiplicityInfo = collections.namedtuple(
    'MultiplicityInfo',
    ["number", "control", "choices", "x_axis_label"],
)

#: set of defined models (standard and custom)
MODELS = {}  # type: Dict[str, SasviewModelType]
# TODO: remove unused MODEL_BY_PATH cache once sasview no longer references it
#: custom model {path: model} mapping so we can check timestamps
MODEL_BY_PATH = {}  # type: Dict[str, SasviewModelType]
#: Track modules that we have loaded so we can determine whether the model
#: has changed since we last reloaded.
_CACHED_MODULE = {}  # type: Dict[str, "module"]

def reset_environment():
    # type: () -> None
    """
    Clear the compute engine context so that the GUI can change devices.

    This removes all compiled kernels, even those that are active on fit
    pages, but they will be restored the next time they are needed.
    """
    kernelcl.reset_environment()
    for model in MODELS.values():
        model._model = None

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
    for name in core.list_models():
        try:
            MODELS[name] = _make_standard_model(name)
        except Exception:
            logger.error(traceback.format_exc())
    if SUPPORT_OLD_STYLE_PLUGINS:
        _register_old_models()

    return list(MODELS.values())


def load_custom_model(path):
    # type: (str) -> SasviewModelType
    """
    Load a custom model given the model path.
    """
    #logger.info("Loading model %s", path)

    # Load the kernel module.  This may already be cached by the loader, so
    # only requires checking the timestamps of the dependents.
    kernel_module = custom.load_custom_kernel_module(path)

    # Check if the module has changed since we last looked.
    reloaded = kernel_module != _CACHED_MODULE.get(path, None)
    _CACHED_MODULE[path] = kernel_module

    # Turn the module into a model.  We need to do this in even if the
    # model has already been loaded so that we can determine the model
    # name and retrieve it from the MODELS cache.
    model = getattr(kernel_module, 'Model', None)
    if model is not None:
        # Old style models do not set the name in the class attributes, so
        # set it here; this name will be overridden when the object is created
        # with an instance variable that has the same value.
        if model.name == "":
            model.name = splitext(basename(path))[0]
        if not hasattr(model, 'filename'):
            model.filename = abspath(kernel_module.__file__).replace('.pyc', '.py')
        if not hasattr(model, 'id'):
            model.id = splitext(basename(model.filename))[0]
    else:
        model_info = modelinfo.make_model_info(kernel_module)
        model = make_model_from_info(model_info)

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
        logger.info("Model %s already exists: using %s [%s]",
                    _previous_name, model.name, model.filename)

    # Only update the model if the module has changed
    if reloaded or model.name not in MODELS:
        MODELS[model.name] = model

    return MODELS[model.name]


def make_model_from_info(model_info):
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


def _make_standard_model(name):
    # type: (str) -> SasviewModelType
    """
    Load the sasview model defined by *name*.

    *name* can be a standard model name or a path to a custom model.

    Returns a class that can be used directly as a sasview model.
    """
    kernel_module = generate.load_kernel_module(name)
    model_info = modelinfo.make_model_info(kernel_module)
    return make_model_from_info(model_info)


def _register_old_models():
    # type: () -> None
    """
    Place the new models into sasview under the old names.

    Monkey patch sas.sascalc.fit as sas.models so that sas.models.pluginmodel
    is available to the plugin modules.
    """
    import sys
    import sas   # needed in order to set sas.models
    import sas.sascalc.fit
    sys.modules['sas.models'] = sas.sascalc.fit
    sas.models = sas.sascalc.fit
    import sas.models
    from sasmodels.conversion_table import CONVERSION_TABLE

    for new_name, conversion in CONVERSION_TABLE.get((3, 1, 2), {}).items():
        # CoreShellEllipsoidModel => core_shell_ellipsoid:1
        new_name = new_name.split(':')[0]
        old_name = conversion[0] if len(conversion) < 3 else conversion[2]
        module_attrs = {old_name: find_model(new_name)}
        ConstructedModule = type(old_name, (), module_attrs)
        old_path = 'sas.models.' + old_name
        setattr(sas.models, old_path, ConstructedModule)
        sys.modules[old_path] = ConstructedModule


def MultiplicationModel(form_factor, structure_factor):
    # type: ("SasviewModel", "SasviewModel") -> "SasviewModel"
    """
    Returns a constructed product model from form_factor and structure_factor.
    """
    model_info = product.make_product_info(form_factor._model_info,
                                           structure_factor._model_info)
    ConstructedModel = make_model_from_info(model_info)
    return ConstructedModel(form_factor.multiplicity)


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
    for p in model_info.parameters.user_parameters({}, is2d=True):
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
    orientation_params = None # type: List[str]
    #: names of the magnetic parameters in the order they appear
    magnetic_params = None    # type: List[str]
    #: names of the fittable parameters
    fixed = None              # type: List[str]
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
    #: Multiplicity information
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

        # Get the list of hidden parameters given the multiplicity
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

        # Update the parameter lists to exclude any hidden parameters
        self.magnetic_params = tuple(pname for pname in self.magnetic_params
                                     if pname not in hidden)
        self.orientation_params = tuple(pname for pname in self.orientation_params
                                        if pname not in hidden)

        self._persistency_dict = {}
        self.params = collections.OrderedDict()
        self.dispersion = collections.OrderedDict()
        self.details = {}
        for p in self._model_info.parameters.user_parameters({}, is2d=True):
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
            return self.calculate_Iq(qx, qy)

        elif isinstance(qdist, np.ndarray):
            # We have a simple 1D distribution of q-values
            return self.calculate_Iq(qdist)

        else:
            raise TypeError("evalDistribution expects q or [qx, qy], not %r"
                            % type(qdist))

    def calc_composition_models(self, qx):
        """
        returns parts of the composition model or None if not a composition
        model.
        """
        # TODO: have calculate_Iq return the intermediates.
        #
        # The current interface causes calculate_Iq() to be called twice,
        # once to get the combined result and again to get the intermediate
        # results.  This is necessary for now.
        # Long term, the solution is to change the interface to calculate_Iq
        # so that it returns a results object containing all the bits:
        #     the A, B, C, ... of the composition model (and any subcomponents?)
        #     the P and S of the product model,
        #     the combined model before resolution smearing,
        #     the sasmodel before sesans conversion,
        #     the oriented 2D model used to fit oriented usans data,
        #     the final I(q),
        #     ...
        #
        # Have the model calculator add all of these blindly to the data
        # tree, and update the graphs which contain them.  The fitter
        # needs to be updated to use the I(q) value only, ignoring the rest.
        #
        # The simple fix of returning the existing intermediate results
        # will not work for a couple of reasons: (1) another thread may
        # sneak in to compute its own results before calc_composition_models
        # is called, and (2) calculate_Iq is currently called three times:
        # once with q, once with q values before qmin and once with q values
        # after q max.  Both of these should be addressed before
        # replacing this code.
        composition = self._model_info.composition
        if composition and composition[0] == 'product': # only P*S for now
            with calculation_lock:
                self._calculate_Iq(qx)
                return self._intermediate_results
        else:
            return None

    def calculate_Iq(self, qx, qy=None):
        # type: (Sequence[float], Optional[Sequence[float]]) -> np.ndarray
        """
        Calculate Iq for one set of q with the current parameters.

        If the model is 1D, use *q*.  If 2D, use *qx*, *qy*.

        This should NOT be used for fitting since it copies the *q* vectors
        to the card for each evaluation.
        """
        ## uncomment the following when trying to debug the uncoordinated calls
        ## to calculate_Iq
        #if calculation_lock.locked():
        #    logger.info("calculation waiting for another thread to complete")
        #    logger.info("\n".join(traceback.format_stack()))

        with calculation_lock:
            return self._calculate_Iq(qx, qy)

    def _calculate_Iq(self, qx, qy=None):
        if self._model is None:
            # Only need one copy of the compiled kernel regardless of how many
            # times it is used, so store it in the class.  Also, to reset the
            # compute engine, need to clear out all existing compiled kernels,
            # which is much easier to do if we store them in the class.
            self.__class__._model = core.build_model(self._model_info)
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
        #print("================ parameters ==================")
        #for p, v in zip(parameters.call_parameters, pairs): print(p.name, v[0])
        #for k, p in enumerate(self._model_info.parameters.call_parameters):
        #    print(k, p.name, *pairs[k])
        #print("params", self.params)
        #print("values", values)
        #print("is_mag", is_magnetic)
        result = calculator(call_details, values, cutoff=self.cutoff,
                            magnetic=is_magnetic)
        #print("result", result)
        self._intermediate_results = getattr(calculator, 'results', None)
        calculator.release()
        #self._model.release()
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
        # type: (str, weights.Dispersion) -> None
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
            raise ValueError("%r is not a dispersity or orientation parameter"
                             % parameter)

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
                return self.multiplicity, [self.multiplicity], [1.0]
            else:
                # For hidden parameters use default values.  This sets
                # scale=1 and background=0 for structure factors
                default = self._model_info.parameters.defaults.get(par.name, np.NaN)
                return default, [default], [1.0]
        elif par.polydisperse:
            value = self.params[par.name]
            dis = self.dispersion[par.name]
            if dis['type'] == 'array':
                dispersity, weight = dis['values'], dis['weights']
            else:
                dispersity, weight = weights.get_weights(
                    dis['type'], dis['npts'], dis['width'], dis['nsigmas'],
                    value, par.limits, par.relative_pd)
            return value, dispersity, weight
        else:
            value = self.params[par.name]
            return value, [value], [1.0]

    @classmethod
    def runTests(cls):
        """
        Run any tests built into the model and captures the test output.

        Returns success flag and output
        """
        from .model_test import check_model
        return check_model(cls._model_info)

def test_cylinder():
    # type: () -> float
    """
    Test that the cylinder model runs, returning the value at [0.1,0.1].
    """
    Cylinder = _make_standard_model('cylinder')
    cylinder = Cylinder()
    return cylinder.evalDistribution([0.1, 0.1])

def test_structure_factor():
    # type: () -> float
    """
    Test that 2-D hardsphere model runs and doesn't produce NaN.
    """
    Model = _make_standard_model('hardsphere')
    model = Model()
    value2d = model.evalDistribution([0.1, 0.1])
    value1d = model.evalDistribution(np.array([0.1*np.sqrt(2)]))
    #print("hardsphere", value1d, value2d)
    if np.isnan(value1d) or np.isnan(value2d):
        raise ValueError("hardsphere returns nan")

def test_product():
    # type: () -> float
    """
    Test that 2-D hardsphere model runs and doesn't produce NaN.
    """
    S = _make_standard_model('hayter_msa')()
    P = _make_standard_model('cylinder')()
    model = MultiplicationModel(P, S)
    value = model.evalDistribution([0.1, 0.1])
    if np.isnan(value):
        raise ValueError("cylinder*hatyer_msa returns null")

def test_rpa():
    # type: () -> float
    """
    Test that the 2-D RPA model runs
    """
    RPA = _make_standard_model('rpa')
    rpa = RPA(3)
    return rpa.evalDistribution([0.1, 0.1])

def test_empty_distribution():
    # type: () -> None
    """
    Make sure that sasmodels returns NaN when there are no polydispersity points
    """
    Cylinder = _make_standard_model('cylinder')
    cylinder = Cylinder()
    cylinder.setParam('radius', -1.0)
    cylinder.setParam('background', 0.)
    Iq = cylinder.evalDistribution(np.asarray([0.1]))
    assert Iq[0] == 0., "empty distribution fails"

def test_model_list():
    # type: () -> None
    """
    Make sure that all models build as sasview models
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
    Load and run cylinder model as sas-models-CylinderModel
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

def magnetic_demo():
    Model = _make_standard_model('sphere')
    model = Model()
    model.setParam('sld_M0', 8)
    q = np.linspace(-0.35, 0.35, 500)
    qx, qy = np.meshgrid(q, q)
    result = model.calculate_Iq(qx.flatten(), qy.flatten())
    result = result.reshape(qx.shape)

    import pylab
    pylab.imshow(np.log(result + 0.001))
    pylab.show()

if __name__ == "__main__":
    print("cylinder(0.1,0.1)=%g"%test_cylinder())
    #magnetic_demo()
    #test_product()
    #test_structure_factor()
    #print("rpa:", test_rpa())
    #test_empty_distribution()
