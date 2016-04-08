"""
Model Info and Parameter Tables
===============================

Defines :class:`ModelInfo` and :class:`ParameterTable` and the routines for
manipulating them.  In particular, :func:`make_model_info` converts a kernel
module into the model info block as seen by the rest of the sasmodels library.
"""
from copy import copy
from os.path import abspath, basename, splitext

import numpy as np

from .details import mono_details

# Optional typing
try:
    from typing import Tuple, List, Union, Dict, Optional, Any, Callable
except ImportError:
    pass
else:
    from .details import CallDetails
    Limits = Tuple[float, float]
    LimitsOrChoice = Union[Limits, Tuple[str]]
    ParameterDef = Tuple[str, str, float, LimitsOrChoice, str, str]
    ParameterSetUser = Dict[str, Union[float, List[float]]]
    ParameterSet = Dict[str, float]
    TestInput = Union[str, float, List[float], Tuple[float, float], List[Tuple[float, float]]]
    TestValue = Union[float, List[float]]
    TestCondition = Tuple[ParameterSetUser, TestInput, TestValue]

MAX_PD = 4 #: Maximum number of simultaneously polydisperse parameters

# assumptions about common parameters exist throughout the code, such as:
# (1) kernel functions Iq, Iqxy, form_volume, ... don't see them
# (2) kernel drivers assume scale is par[0] and background is par[1]
# (3) mixture models drop the background on components and replace the scale
#     with a scale that varies from [-inf, inf]
# (4) product models drop the background and reassign scale
# and maybe other places.
# Note that scale and background cannot be coordinated parameters whose value
# depends on the some polydisperse parameter with the current implementation
COMMON_PARAMETERS = [
    ["scale", "", 1, [0, np.inf], "", "Source intensity"],
    ["background", "1/cm", 1e-3, [0, np.inf], "", "Source background"],
]
assert (len(COMMON_PARAMETERS) == 2
        and COMMON_PARAMETERS[0][0]=="scale"
        and COMMON_PARAMETERS[1][0]=="background"), "don't change common parameters"


def make_parameter_table(pars):
    # type: (List[ParameterDefinition) -> ParameterTable
    """
    Construct a parameter table from a list of parameter definitions.

    This is used by the module processor to convert the parameter block into
    the parameter table seen in the :class:`ModelInfo` for the module.
    """
    processed = []
    for p in pars:
        if not isinstance(p, (list, tuple)) or len(p) != 6:
            raise ValueError("Parameter should be [name, units, default, limits, type, desc], but got %r"
                             %str(p))
        processed.append(parse_parameter(*p))
    partable = ParameterTable(processed)
    return partable

def parse_parameter(name, units='', default=np.NaN,
                    limits=(-np.inf, np.inf), ptype='', description=''):
    """
    Parse an individual parameter from the parameter definition block.

    This does type and value checking on the definition, leading
    to early failure in the model loading process and easier debugging.
    """
    # type: (str, str, float, LimitsOrChoice, str, str) -> Parameter
    # Parameter is a user facing class.  Do robust type checking.
    if not isstr(name):
        raise ValueError("expected string for parameter name %r"%name)
    if not isstr(units):
        raise ValueError("expected units to be a string for %s"%name)
    # if limits is a list of strings, then this is a choice list
    # field, and limits are 1 to length of string list
    if isinstance(limits, list) and all(isstr(k) for k in limits):
        choices = limits
        limits = [0, len(choices)-1]
    else:
        choices = []
    # TODO: maybe allow limits of None for (-inf, inf)
    try:
        low, high = limits
        if not isinstance(low, (int, float)):
            raise TypeError("low is not numeric")
        if not isinstance(high, (int, float)):
            raise TypeError("high is not numeric")
        if low >= high:
            raise ValueError("require low < high")
    except Exception:
        raise ValueError("invalid limits %s for %s"%(limits, name))

    if not isinstance(default, (int, float)):
        raise ValueError("expected default %r to be a number for %s"
                         % (default, name))
    if default < low or default > high:
        raise ValueError("default value %r not in range for %s"
                         % (default, name))

    if ptype not in ("volume", "orientation", "sld", "magnetic", ""):
        raise ValueError("unexpected type %r for %s" % (ptype, name))

    if not isstr(description):
        raise ValueError("expected description to be a string")


    # Parameter id for name[n] does not include [n]
    if "[" in name:
        if not name.endswith(']'):
            raise ValueError("Expected name[len] for vector parameter %s"%name)
        pid, ref = name[:-1].split('[', 1)
        ref = ref.strip()
    else:
        pid, ref = name, None


    # automatically identify sld types
    if ptype== '' and (pid.startswith('sld') or pid.endswith('sld')):
        ptype = 'sld'

    # Check if using a vector definition, name[k], as the parameter name
    if ref:
        if ref == '':
            raise ValueError("Need to specify vector length for %s"%name)
        try:
            length = int(ref)
            control = None
        except ValueError:
            length = None
            control = ref
    else:
        length = 1
        control = None

    # Build the parameter
    parameter = Parameter(name=name, units=units, default=default,
                          limits=limits, ptype=ptype, description=description)

    # TODO: need better control over whether a parameter is polydisperse
    parameter.polydisperse = ptype in ('orientation', 'volume')
    parameter.relative_pd = ptype == 'volume'
    parameter.choices = choices
    parameter.length = length
    parameter.length_control = control

    return parameter


def expand_pars(partable, pars):
    # type: (ParameterTable, ParameterSetUser) ->  ParameterSet
    """
    Create demo parameter set from key-value pairs.

    *pars* are the key-value pairs to use for the parameters.  Any
    parameters not specified in *pars* are set from the *partable* defaults.

    If *pars* references vector fields, such as thickness[n], then support
    different ways of assigning the demo values, including assigning a
    specific value (e.g., thickness3=50.0), assigning a new value to all
    (e.g., thickness=50.0) or assigning values using list notation.
    """
    if pars is None:
        result = partable.defaults
    else:
        lookup = dict((p.id, p) for p in partable.kernel_parameters)
        result = partable.defaults.copy()
        vectors = dict((name,value) for name,value in pars.items()
                       if name in lookup and lookup[name].length > 1)
        if vectors:
            scalars = dict((name, value) for name, value in pars.items()
                           if name not in lookup or lookup[name].length == 1)
            for name, value in vectors.items():
                if np.isscalar(value):
                    # support for the form
                    #    dict(thickness=0, thickness2=50)
                    for k in range(1, lookup[name].length+1):
                        key = name+str(k)
                        if key not in scalars:
                            scalars[key] = vectors
                else:
                    # supoprt for the form
                    #    dict(thickness=[20,10,3])
                    for (k,v) in enumerate(value):
                        scalars[name+str(k)] = v
            result.update(scalars)
        else:
            result.update(pars)

    return result

def prefix_parameter(par, prefix):
    # type: (Parameter, str) -> Parameter
    """
    Return a copy of the parameter with its name prefixed.
    """
    new_par = copy(par)
    new_par.name = prefix + par.name
    new_par.id = prefix + par.id

def suffix_parameter(par, suffix):
    # type: (Parameter, str) -> Parameter
    """
    Return a copy of the parameter with its name prefixed.
    """
    new_par = copy(par)
    # If name has the form x[n], replace with x_suffix[n]
    new_par.name = par.id + suffix + par.name[len(par.id):]
    new_par.id = par.id + suffix

class Parameter(object):
    """
    The available kernel parameters are defined as a list, with each parameter
    defined as a sublist with the following elements:

    *name* is the name that will be used in the call to the kernel
    function and the name that will be displayed to the user.  Names
    should be lower case, with words separated by underscore.  If
    acronyms are used, the whole acronym should be upper case.

    *units* should be one of *degrees* for angles, *Ang* for lengths,
    *1e-6/Ang^2* for SLDs.

    *default value* will be the initial value for  the model when it
    is selected, or when an initial value is not otherwise specified.

    *limits = [lb, ub]* are the hard limits on the parameter value, used to
    limit the polydispersity density function.  In the fit, the parameter limits
    given to the fit are the limits  on the central value of the parameter.
    If there is polydispersity, it will evaluate parameter values outside
    the fit limits, but not outside the hard limits specified in the model.
    If there are no limits, use +/-inf imported from numpy.

    *type* indicates how the parameter will be used.  "volume" parameters
    will be used in all functions.  "orientation" parameters will be used
    in *Iqxy* and *Imagnetic*.  "magnetic* parameters will be used in
    *Imagnetic* only.  If *type* is the empty string, the parameter will
    be used in all of *Iq*, *Iqxy* and *Imagnetic*.  "sld" parameters
    can automatically be promoted to magnetic parameters, each of which
    will have a magnitude and a direction, which may be different from
    other sld parameters. The volume parameters are used for calls
    to form_volume within the kernel (required for volume normalization)
    and for calls to ER and VR for effective radius and volume ratio
    respectively.

    *description* is a short description of the parameter.  This will
    be displayed in the parameter table and used as a tool tip for the
    parameter value in the user interface.

    Additional values can be set after the parameter is created:

    * *length* is the length of the field if it is a vector field

    * *length_control* is the parameter which sets the vector length

    * *is_control* is True if the parameter is a control parameter for a vector

    * *polydisperse* is true if the parameter accepts a polydispersity

    * *relative_pd* is true if that polydispersity is a portion of the
    value (so a 10% length dipsersity would use a polydispersity value of 0.1)
    rather than absolute dispersisity (such as an angle plus or minus
    15 degrees).

    In the usual process these values are set by :func:`make_parameter_table`
    and :func:`parse_parameter` therein.
    """
    def __init__(self, name, units='', default=None, limits=(-np.inf, np.inf),
                 ptype='', description=''):
        # type: (str, str, float, Limits, str, str)
        self.id = name.split('[')[0].strip() # type: str
        self.name = name                     # type: str
        self.units = units                   # type: str
        self.default = default               # type: float
        self.limits = limits                 # type: Limits
        self.type = ptype                    # type: str
        self.description = description       # type: str

        # Length and length_control will be filled in once the complete
        # parameter table is available.
        self.length = 1                      # type: int
        self.length_control = None           # type: Optional[str]
        self.is_control = False              # type: bool

        # TODO: need better control over whether a parameter is polydisperse
        self.polydisperse = False            # type: bool
        self.relative_pd = False             # type: bool

        # choices are also set externally.
        self.choices = []                    # type: List[str]

    def as_definition(self):
        # type: () -> str
        """
        Declare space for the variable in a parameter structure.

        For example, the parameter thickness with length 3 will
        return "double thickness[3];", with no spaces before and
        no new line character afterward.
        """
        if self.length == 1:
            return "double %s;"%self.id
        else:
            return "double %s[%d];"%(self.id, self.length)

    def as_function_argument(self):
        # type: () -> str
        """
        Declare the variable as a function argument.

        For example, the parameter thickness with length 3 will
        return "double *thickness", with no spaces before and
        no comma afterward.
        """
        if self.length == 1:
            return "double %s"%self.id
        else:
            return "double *%s"%self.id

    def as_call_reference(self, prefix=""):
        # type: () -> str
        # Note: if the parameter is a struct type, then we will need to use
        # &prefix+id.  For scalars and vectors we can just use prefix+id.
        return prefix + self.id

    def __str__(self):
        # type: () -> str
        return "<%s>"%self.name

    def __repr__(self):
        # type: () -> str
        return "P<%s>"%self.name


class ParameterTable(object):
    """
    ParameterTable manages the list of available parameters.

    There are a couple of complications which mean that the list of parameters
    for the kernel differs from the list of parameters that the user sees.

    (1) Common parameters.  Scale and background are implicit to every model,
    but are not passed to the kernel.

    (2) Vector parameters.  Vector parameters are passed to the kernel as a
    pointer to an array, e.g., thick[], but they are seen by the user as n
    separate parameters thick1, thick2, ...

    Therefore, the parameter table is organized by how it is expected to be
    used. The following information is needed to set up the kernel functions:

    * *kernel_parameters* is the list of parameters in the kernel parameter
    table, with vector parameter p declared as p[].

    * *iq_parameters* is the list of parameters to the Iq(q, ...) function,
    with vector parameter p sent as p[].

    * *iqxy_parameters* is the list of parameters to the Iqxy(qx, qy, ...)
    function, with vector parameter p sent as p[].

    * *form_volume_parameters* is the list of parameters to the form_volume(...)
    function, with vector parameter p sent as p[].

    Problem details, which sets up the polydispersity loops, requires the
    following:

    * *theta_offset* is the offset of the theta parameter in the kernel parameter
    table, with vector parameters counted as n individual parameters
    p1, p2, ..., or offset is -1 if there is no theta parameter.

    * *max_pd* is the maximum number of polydisperse parameters, with vector
    parameters counted as n individual parameters p1, p2, ...  Note that
    this number is limited to sasmodels.modelinfo.MAX_PD.

    * *npars* is the total number of parameters to the kernel, with vector
    parameters counted as n individual parameters p1, p2, ...

    * *call_parameters* is the complete list of parameters to the kernel,
    including scale and background, with vector parameters recorded as
    individual parameters p1, p2, ...

    * *active_1d* is the set of names that may be polydisperse for 1d data

    * *active_2d* is the set of names that may be polydisperse for 2d data

    User parameters are the set of parameters visible to the user, including
    the scale and background parameters that the kernel does not see.  User
    parameters don't use vector notation, and instead use p1, p2, ...

    * *control_parameters* is the

    """
    # scale and background are implicit parameters
    COMMON = [Parameter(*p) for p in COMMON_PARAMETERS]

    def __init__(self, parameters):
        # type: (List[Parameter]) -> None
        self.kernel_parameters = parameters
        self._set_vector_lengths()
        self.call_parameters = self._get_call_parameters()
        self.defaults = self._get_defaults()
        #self._name_table= dict((p.id, p) for p in parameters)

        # Set the kernel parameters.  Assumes background and scale are the
        # first two parameters in the parameter list, but these are not sent
        # to the underlying kernel functions.
        self.iq_parameters = [p for p in self.kernel_parameters
                              if p.type not in ('orientation', 'magnetic')]
        self.iqxy_parameters = [p for p in self.kernel_parameters
                                if p.type != 'magnetic']
        self.form_volume_parameters = [p for p in self.kernel_parameters
                                       if p.type == 'volume']

        # Theta offset
        offset = 0
        for p in self.kernel_parameters:
            if p.name == 'theta':
                self.theta_offset = offset
                break
            offset += p.length
        else:
            self.theta_offset = -1

        # number of polydisperse parameters
        num_pd = sum(p.length for p in self.kernel_parameters if p.polydisperse)
        # Don't use more polydisperse parameters than are available in the model
        # Note: we can do polydispersity on arbitrary parameters, so it is not
        # clear that this is a good idea; it does however make the poly_details
        # code easier to write, so we will leave it in for now.
        self.max_pd = min(num_pd, MAX_PD)

        self.npars = sum(p.length for p in self.kernel_parameters)

        # true if has 2D parameters
        self.has_2d = any(p.type in ('orientation', 'magnetic')
                          for p in self.kernel_parameters)

        self.pd_1d = set(p.name for p in self.call_parameters
                         if p.polydisperse and p.type not in ('orientation', 'magnetic'))
        self.pd_2d = set(p.name for p in self.call_parameters
                         if p.polydisperse and p.type != 'magnetic')


    def _set_vector_lengths(self):
        # type: () -> None
        """
        Walk the list of kernel parameters, setting the length field of the
        vector parameters from the upper limit of the reference parameter.

        This needs to be done once the entire parameter table is available
        since the reference may still be undefined when the parameter is
        initially created.

        Note: This modifies the underlying parameter object.
        """
        # Sort out the length of the vector parameters such as thickness[n]
        for p in self.kernel_parameters:
            if p.length_control:
                for ref in self.kernel_parameters:
                    if ref.id == p.length_control:
                        break
                else:
                    raise ValueError("no reference variable %r for %s"
                                     % (p.length_control, p.name))
                ref.is_control = True
                low, high = ref.limits
                if int(low) != low or int(high) != high or low < 0 or high > 20:
                    raise ValueError("expected limits on %s to be within [0, 20]"
                                     % ref.name)
                # TODO: may want to make a copy of the parameter before updating
                # this introduces other potential problems, since the same
                # parameter may be referenced elsewhere
                p.length = high

    def _get_defaults(self):
        # type: () -> ParameterSet
        """
        Get a list of parameter defaults from the parameters.

        Expands vector parameters into parameter id+number.
        """
        # Construct default values, including vector defaults
        defaults = {}
        for p in self.call_parameters:
            if p.length == 1:
                defaults[p.id] = p.default
            else:
                for k in range(1, p.length+1):
                    defaults["%s%d"%(p.id, k)] = p.default
        return defaults

    def _get_call_parameters(self):
        # type: () -> List[Parameter]
        full_list = self.COMMON[:]
        for p in self.kernel_parameters:
            if p.length == 1:
                full_list.append(p)
            else:
                for k in range(1, p.length+1):
                    pk = Parameter(p.id+str(k), p.units, p.default,
                                   p.limits, p.type, p.description)
                    pk.polydisperse = p.polydisperse
                    pk.relative_pd = p.relative_pd
                    full_list.append(pk)
        return full_list

    def user_parameters(self, pars={}, is2d=True):
        # type: (Dict[str, float], bool) -> List[Parameter]
        """
        Return the list of parameters for the given data type.

        Vector parameters are expanded as in place.  If multiple parameters
        share the same vector length, then the parameters will be interleaved
        in the result.  The control parameters come first.  For example,
        if the parameter table is ordered as::

            sld_core
            sld_shell[num_shells]
            sld_solvent
            thickness[num_shells]
            num_shells

        and *pars[num_shells]=2* then the returned list will be::

            num_shells
            scale
            background
            sld_core
            sld_shell1
            thickness1
            sld_shell2
            thickness2
            sld_solvent

        Note that shell/thickness pairs are grouped together in the result
        even though they were not grouped in the incoming table.  The control
        parameter is always returned first since the GUI will want to set it
        early, and rerender the table when it is changed.
        """
        control = [p for p in self.kernel_parameters if p.is_control]

        # Gather entries such as name[n] into groups of the same n
        dependent = dict((p.id, []) for p in control)
        for p in self.kernel_parameters:
            if p.length_control is not None:
                dependent[p.length_control].append(p)

        # Gather entries such as name[4] into groups of the same length
        fixed = {}
        for p in self.kernel_parameters:
            if p.length > 1 and p.length_control is None:
                fixed.setdefault(p.length, []).append(p)

        # Using the call_parameters table, we already have expanded forms
        # for each of the vector parameters; put them in a lookup table
        expanded_pars = dict((p.name, p) for p in self.call_parameters)

        # Gather the user parameters in order
        result = control + self.COMMON
        for p in self.kernel_parameters:
            if not is2d and p.type in ('orientation', 'magnetic'):
                pass
            elif p.is_control:
                pass # already added
            elif p.length_control is not None:
                table = dependent.get(p.length_control, [])
                if table:
                    # look up length from incoming parameters
                    table_length = int(pars.get(p.length_control, p.length))
                    del dependent[p.length_control] # first entry seen
                    for k in range(1, table_length+1):
                        for entry in table:
                            result.append(expanded_pars[entry.id+str(k)])
                else:
                    pass # already processed all entries
            elif p.length > 1:
                table = fixed.get(p.length, [])
                if table:
                    table_length = p.length
                    del fixed[p.length]
                    for k in range(1, table_length+1):
                        for entry in table:
                            result.append(expanded_pars[entry.id+str(k)])
                else:
                    pass # already processed all entries
            else:
                result.append(p)

        return result

def isstr(x):
    # type: (Any) -> bool
    """
    Return True if the object is a string.
    """
    # TODO: 2-3 compatible tests for str, including unicode strings
    return isinstance(x, str)

def make_model_info(kernel_module):
    # type: (module) -> ModelInfo
    """
    Extract the model definition from the loaded kernel module.

    Fill in default values for parts of the module that are not provided.

    Note: vectorized Iq and Iqxy functions will be created for python
    models when the model is first called, not when the model is loaded.
    """
    info = ModelInfo()
    #print("make parameter table", kernel_module.parameters)
    parameters = make_parameter_table(kernel_module.parameters)
    demo = expand_pars(parameters, getattr(kernel_module, 'demo', None))
    filename = abspath(kernel_module.__file__)
    kernel_id = splitext(basename(filename))[0]
    name = getattr(kernel_module, 'name', None)
    if name is None:
        name = " ".join(w.capitalize() for w in kernel_id.split('_'))

    info.id = kernel_id  # string used to load the kernel
    info.filename = abspath(kernel_module.__file__)
    info.name = name
    info.title = getattr(kernel_module, 'title', name+" model")
    info.description = getattr(kernel_module, 'description', 'no description')
    info.parameters = parameters
    info.demo = demo
    info.composition = None
    info.docs = kernel_module.__doc__
    info.category = getattr(kernel_module, 'category', None)
    info.single = getattr(kernel_module, 'single', True)
    info.structure_factor = getattr(kernel_module, 'structure_factor', False)
    info.profile_axes = getattr(kernel_module, 'profile_axes', ['x','y'])
    info.variant_info = getattr(kernel_module, 'variant_info', None)
    info.source = getattr(kernel_module, 'source', [])
    # TODO: check the structure of the tests
    info.tests = getattr(kernel_module, 'tests', [])
    info.ER = getattr(kernel_module, 'ER', None)
    info.VR = getattr(kernel_module, 'VR', None)
    info.form_volume = getattr(kernel_module, 'form_volume', None)
    info.Iq = getattr(kernel_module, 'Iq', None)
    info.Iqxy = getattr(kernel_module, 'Iqxy', None)
    info.profile = getattr(kernel_module, 'profile', None)
    info.sesans = getattr(kernel_module, 'sesans', None)

    # Precalculate the monodisperse parameter details
    info.mono_details = mono_details(info)
    return info

class ModelInfo(object):
    """
    Interpret the model definition file, categorizing the parameters.

    The module can be loaded with a normal python import statement if you
    know which module you need, or with __import__('sasmodels.model.'+name)
    if the name is in a string.

    The *model_info* structure contains the following fields:

    * *id* is the id of the kernel
    * *name* is the display name of the kernel
    * *filename* is the full path to the module defining the file (if any)
    * *title* is a short description of the kernel
    * *description* is a long description of the kernel (this doesn't seem
      very useful since the Help button on the model page brings you directly
      to the documentation page)
    * *docs* is the docstring from the module.  Use :func:`make_doc` to
    * *category* specifies the model location in the docs
    * *parameters* is the model parameter table
    * *single* is True if the model allows single precision
    * *structure_factor* is True if the model is useable in a product
    * *variant_info* contains the information required to select between
      model variants (e.g., the list of cases) or is None if there are no
      model variants
    * *par_type* categorizes the model parameters. See
      :func:`categorize_parameters` for details.
    * *demo* contains the *{parameter: value}* map used in compare (and maybe
      for the demo plot, if plots aren't set up to use the default values).
      If *demo* is not given in the file, then the default values will be used.
    * *tests* is a set of tests that must pass
    * *source* is the list of library files to include in the C model build
    * *Iq*, *Iqxy*, *form_volume*, *ER*, *VR* and *sesans* are python functions
      implementing the kernel for the module, or None if they are not
      defined in python
    * *composition* is None if the model is independent, otherwise it is a
      tuple with composition type ('product' or 'mixture') and a list of
      *model_info* blocks for the composition objects.  This allows us to
      build complete product and mixture models from just the info.

    The structure should be mostly static, other than the delayed definition
    of *Iq* and *Iqxy* if they need to be defined.
    """
    id = None               # type: str
    filename = None         # type: str
    name = None             # type: str
    title = None            # type: str
    description = None      # type: str
    parameters = None       # type: ParameterTable
    demo = None             # type: Dict[str, float]
    composition = None      # type: Optional[Tuple[str, List[ModelInfo]]]
    docs = None             # type: str
    category = None         # type: Optional[str]
    single = None           # type: bool
    structure_factor = None # type: bool
    profile_axes = None     # type: Tuple[str, str]
    variant_info = None     # type: Optional[List[str]]
    source = None           # type: List[str]
    tests = None            # type: List[TestCondition]
    ER = None               # type: Optional[Callable[[np.ndarray, ...], np.ndarray]]
    VR = None               # type: Optional[Callable[[np.ndarray, ...], Tuple[np.ndarray, np.ndarray]]]
    form_volume = None      # type: Optional[Callable[[np.ndarray, ...], float]]
    Iq = None               # type: Optional[Callable[[np.ndarray, ...], np.ndarray]]
    Iqxy = None             # type: Optional[Callable[[np.ndarray, ...], np.ndarray]]
    profile = None
    sesans = None
    mono_details = None     # type: CallDetails

    def __init__(self):
        # type: () -> None
        pass


