"""
Model Info and Parameter Tables
===============================

Defines :class:`ModelInfo` and :class:`ParameterTable` and the routines for
manipulating them.  In particular, :func:`make_model_info` converts a kernel
module into the model info block as seen by the rest of the sasmodels library.
"""
from __future__ import print_function

from copy import copy
from os.path import abspath, basename, splitext
import inspect

import numpy as np  # type: ignore

# Optional typing
# pylint: disable=unused-import
try:
    from typing import Tuple, List, Union, Dict, Optional, Any, Callable, Sequence, Set
    from types import ModuleType
except ImportError:
    pass
else:
    Limits = Tuple[float, float]
    #LimitsOrChoice = Union[Limits, Tuple[Sequence[str]]]
    ParameterDef = Tuple[str, str, float, Limits, str, str]
    ParameterSetUser = Dict[str, Union[float, List[float]]]
    ParameterSet = Dict[str, float]
    TestInput = Union[str, float, List[float], Tuple[float, float], List[Tuple[float, float]]]
    TestValue = Union[float, List[float]]
    TestCondition = Tuple[ParameterSetUser, TestInput, TestValue]
# pylint: enable=unused-import

# If MAX_PD changes, need to change the loop macros in kernel_iq.c
MAX_PD = 5 #: Maximum number of simultaneously polydisperse parameters

# assumptions about common parameters exist throughout the code, such as:
# (1) kernel functions Iq, Iqxy, Iqac, Iqabc, form_volume, ... don't see them
# (2) kernel drivers assume scale is par[0] and background is par[1]
# (3) mixture models drop the background on components and replace the scale
#     with a scale that varies from [-inf, inf]
# (4) product models drop the background and reassign scale
# and maybe other places.
# Note that scale and background cannot be coordinated parameters whose value
# depends on the some polydisperse parameter with the current implementation
DEFAULT_BACKGROUND = 1e-3
COMMON_PARAMETERS = [
    ("scale", "", 1, (0.0, np.inf), "", "Source intensity"),
    ("background", "1/cm", DEFAULT_BACKGROUND, (-np.inf, np.inf), "", "Source background"),
]
assert (len(COMMON_PARAMETERS) == 2
        and COMMON_PARAMETERS[0][0] == "scale"
        and COMMON_PARAMETERS[1][0] == "background"), "don't change common parameters"


def make_parameter_table(pars):
    # type: (List[ParameterDef]) -> ParameterTable
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
    partable.check_angles()
    return partable

def parse_parameter(name, units='', default=np.NaN,
                    user_limits=None, ptype='', description=''):
    # type: (str, str, float, Sequence[Any], str, str) -> Parameter
    """
    Parse an individual parameter from the parameter definition block.

    This does type and value checking on the definition, leading
    to early failure in the model loading process and easier debugging.
    """
    # Parameter is a user facing class.  Do robust type checking.
    if not isstr(name):
        raise ValueError("expected string for parameter name %r"%name)
    if not isstr(units):
        raise ValueError("expected units to be a string for %s"%name)

    # Process limits as [float, float] or [[str, str, ...]]
    choices = []  # type: List[str]
    if user_limits is None:
        limits = (-np.inf, np.inf)
    elif not isinstance(user_limits, (tuple, list)):
        raise ValueError("invalid limits for %s"%name)
    else:
        # if limits is [[str,...]], then this is a choice list field,
        # and limits are 1 to length of string list
        if isinstance(user_limits[0], (tuple, list)):
            choices = user_limits[0]
            limits = (0., len(choices)-1.)
            if not all(isstr(k) for k in choices):
                raise ValueError("choices must be strings for %s"%name)
        else:
            try:
                low, high = user_limits
                limits = (float(low), float(high))
            except Exception:
                raise ValueError("invalid limits for %s: %r"%(name, user_limits))
            if low >= high:
                raise ValueError("require lower limit < upper limit")

    # Process default value as float, making sure it is in range
    if not isinstance(default, (int, float)):
        raise ValueError("expected default %r to be a number for %s"
                         % (default, name))
    if default < limits[0] or default > limits[1]:
        raise ValueError("default value %r not in range for %s"
                         % (default, name))

    # Check for valid parameter type
    if ptype not in ("volume", "orientation", "sld", "magnetic", ""):
        raise ValueError("unexpected type %r for %s" % (ptype, name))

    # Check for valid parameter description
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
    if ptype == '' and (pid.startswith('sld') or pid.endswith('sld')):
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
        scalars = dict((name, value) for name, value in pars.items()
                       if name not in lookup or lookup[name].length == 1)
        vectors = dict((name, value) for name, value in pars.items()
                       if name in lookup and lookup[name].length > 1)
        #print("lookup", lookup)
        #print("scalars", scalars)
        #print("vectors", vectors)
        if vectors:
            for name, value in vectors.items():
                if np.isscalar(value):
                    # support for the form
                    #    dict(thickness=0, thickness2=50)
                    for k in range(1, lookup[name].length+1):
                        key = name+str(k)
                        if key not in scalars:
                            scalars[key] = value
                else:
                    # supoprt for the form
                    #    dict(thickness=[20,10,3])
                    for (k, v) in enumerate(value):
                        scalars[name+str(k+1)] = v
        result.update(scalars)
        #print("expanded", result)

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

    *name* is the name that will be displayed to the user.  Names
    should be lower case, with words separated by underscore.  If
    acronyms are used, the whole acronym should be upper case. For vector
    parameters, the name will be followed by *[len]* where *len* is an
    integer length of the vector, or the name of the parameter which
    controls the length.  The attribute *id* will be created from name
    without the length.

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
    will be used in all functions.  "orientation" parameters are not passed,
    but will be used to convert from *qx*, *qy* to *qa*, *qb*, *qc* in calls to
    *Iqxy* and *Imagnetic*.  If *type* is the empty string, the parameter will
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
      value (so a 10% length dipsersity would use a polydispersity value
      of 0.1) rather than absolute dispersisity (such as an angle plus or
      minus 15 degrees).

    *choices* is the option names for a drop down list of options, as for
    example, might be used to set the value of a shape parameter.

    These values are set by :func:`make_parameter_table` and
    :func:`parse_parameter` therein.
    """
    def __init__(self, name, units='', default=None, limits=(-np.inf, np.inf),
                 ptype='', description=''):
        # type: (str, str, float, Limits, str, str) -> None
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
        r"""
        Declare the variable as a function argument.

        For example, the parameter thickness with length 3 will
        return "double \*thickness", with no spaces before and
        no comma afterward.
        """
        if self.length == 1:
            return "double %s"%self.id
        else:
            return "double *%s"%self.id

    def as_call_reference(self, prefix=""):
        # type: (str) -> str
        """
        Return *prefix* + parameter name.  For struct references, use "v."
        as the prefix.
        """
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
    """
    # scale and background are implicit parameters
    COMMON = [Parameter(*p) for p in COMMON_PARAMETERS]

    def __init__(self, parameters):
        # type: (List[Parameter]) -> None
        self.kernel_parameters = parameters
        self._set_vector_lengths()

        self.npars = sum(p.length for p in self.kernel_parameters)
        self.nmagnetic = sum(p.length for p in self.kernel_parameters
                             if p.type == 'sld')
        self.nvalues = 2 + self.npars
        if self.nmagnetic:
            self.nvalues += 3 + 3*self.nmagnetic

        self.call_parameters = self._get_call_parameters()
        self.defaults = self._get_defaults()
        #self._name_table= dict((p.id, p) for p in parameters)

        # Set the kernel parameters.  Assumes background and scale are the
        # first two parameters in the parameter list, but these are not sent
        # to the underlying kernel functions.
        self.iq_parameters = [p for p in self.kernel_parameters
                              if p.type not in ('orientation', 'magnetic')]
        self.orientation_parameters = [p for p in self.kernel_parameters
                                       if p.type == 'orientation']
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
        self.max_pd = min(num_pd, MAX_PD)

        # true if has 2D parameters
        self.has_2d = any(p.type in ('orientation', 'magnetic')
                          for p in self.kernel_parameters)
        self.is_asymmetric = any(p.name == 'psi' for p in self.kernel_parameters)
        self.magnetism_index = [k for k, p in enumerate(self.call_parameters)
                                if p.id.endswith('_M0')]

        self.pd_1d = set(p.name for p in self.call_parameters
                         if p.polydisperse and p.type not in ('orientation', 'magnetic'))
        self.pd_2d = set(p.name for p in self.call_parameters if p.polydisperse)

    def check_angles(self):
        """
        Check that orientation angles are theta, phi and possibly psi.
        """
        theta = phi = psi = -1
        for k, p in enumerate(self.kernel_parameters):
            if p.name == 'theta':
                theta = k
                if p.type != 'orientation':
                    raise TypeError("theta must be an orientation parameter")
            elif p.name == 'phi':
                phi = k
                if p.type != 'orientation':
                    raise TypeError("phi must be an orientation parameter")
            elif p.name == 'psi':
                psi = k
                if p.type != 'orientation':
                    raise TypeError("psi must be an orientation parameter")
            elif p.type == 'orientation':
                raise TypeError("only theta, phi and psi can be orientation parameters")
        if theta >= 0 and phi >= 0:
            last_par = len(self.kernel_parameters) - 1
            if phi != theta+1:
                raise TypeError("phi must follow theta")
            if psi >= 0 and psi != phi+1:
                raise TypeError("psi must follow phi")
            if (psi >= 0 and psi != last_par) or (psi < 0 and phi != last_par):
                raise TypeError("orientation parameters must appear at the "
                                "end of the parameter table")
        elif theta >= 0 or phi >= 0 or psi >= 0:
            raise TypeError("oriented shapes must have both theta and phi and maybe psi")

    def __getitem__(self, key):
        # Find the parameter definition
        for par in self.call_parameters:
            if par.name == key:
                return par
        raise KeyError("unknown parameter %r"%key)

    def __contains__(self, key):
        for par in self.call_parameters:
            if par.name == key:
                return True
        return False

    def _set_vector_lengths(self):
        # type: () -> List[str]
        """
        Walk the list of kernel parameters, setting the length field of the
        vector parameters from the upper limit of the reference parameter.

        This needs to be done once the entire parameter table is available
        since the reference may still be undefined when the parameter is
        initially created.

        Returns the list of control parameter names.

        Note: This modifies the underlying parameter object.
        """
        # Sort out the length of the vector parameters such as thickness[n]
        for p in self.kernel_parameters:
            if p.length_control:
                ref = self._get_ref(p)
                ref.is_control = True
                ref.polydisperse = False
                low, high = ref.limits
                if int(low) != low or int(high) != high or low < 0 or high > 20:
                    raise ValueError("expected limits on %s to be within [0, 20]"
                                     % ref.name)
                p.length = int(high)

    def _get_ref(self, p):
        # type: (Parameter) -> Parameter
        for ref in self.kernel_parameters:
            if ref.id == p.length_control:
                return ref
        raise ValueError("no reference variable %r for %s"
                         % (p.length_control, p.name))

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
                    pk.choices = p.choices
                    full_list.append(pk)

        # Add the magnetic parameters to the end of the call parameter list.
        if self.nmagnetic > 0:
            full_list.extend([
                Parameter('up_frac_i', '', 0., [0., 1.],
                          'magnetic', 'fraction of spin up incident'),
                Parameter('up_frac_f', '', 0., [0., 1.],
                          'magnetic', 'fraction of spin up final'),
                Parameter('up_angle', 'degrees', 0., [0., 360.],
                          'magnetic', 'spin up angle'),
            ])
            slds = [p for p in full_list if p.type == 'sld']
            for p in slds:
                full_list.extend([
                    Parameter(p.id+'_M0', '1e-6/Ang^2', 0., [-np.inf, np.inf],
                              'magnetic', 'magnetic amplitude for '+p.description),
                    Parameter(p.id+'_mtheta', 'degrees', 0., [-90., 90.],
                              'magnetic', 'magnetic latitude for '+p.description),
                    Parameter(p.id+'_mphi', 'degrees', 0., [-180., 180.],
                              'magnetic', 'magnetic longitude for '+p.description),
                ])
        #print("call parameters", full_list)
        return full_list

    def user_parameters(self, pars, is2d=True):
        # type: (Dict[str, float], bool) -> List[Parameter]
        """
        Return the list of parameters for the given data type.

        Vector parameters are expanded in place.  If multiple parameters
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

        Parameters marked as sld will automatically have a set of associated
        magnetic parameters (p_M0, p_mtheta, p_mphi), as well as polarization
        information (up_theta, up_frac_i, up_frac_f).
        """
        # control parameters go first
        control = [p for p in self.kernel_parameters if p.is_control]

        # Gather entries such as name[n] into groups of the same n
        dependent = {} # type: Dict[str, List[Parameter]]
        dependent.update((p.id, []) for p in control)
        for p in self.kernel_parameters:
            if p.length_control is not None:
                dependent[p.length_control].append(p)

        # Gather entries such as name[4] into groups of the same length
        fixed_length = {}  # type: Dict[int, List[Parameter]]
        for p in self.kernel_parameters:
            if p.length > 1 and p.length_control is None:
                fixed_length.setdefault(p.length, []).append(p)

        # Using the call_parameters table, we already have expanded forms
        # for each of the vector parameters; put them in a lookup table
        # Note: p.id and p.name are currently identical for the call parameters
        expanded_pars = dict((p.id, p) for p in self.call_parameters)

        def append_group(name):
            """add the named parameter, and related magnetic parameters if any"""
            result.append(expanded_pars[name])
            if is2d:
                for tag in '_M0', '_mtheta', '_mphi':
                    if name+tag in expanded_pars:
                        result.append(expanded_pars[name+tag])

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
                            append_group(entry.id+str(k))
                else:
                    pass # already processed all entries
            elif p.length > 1:
                table = fixed_length.get(p.length, [])
                if table:
                    table_length = p.length
                    del fixed_length[p.length]
                    for k in range(1, table_length+1):
                        for entry in table:
                            append_group(entry.id+str(k))
                else:
                    pass # already processed all entries
            else:
                append_group(p.id)

        if is2d and 'up_angle' in expanded_pars:
            result.extend([
                expanded_pars['up_frac_i'],
                expanded_pars['up_frac_f'],
                expanded_pars['up_angle'],
            ])

        return result

def isstr(x):
    # type: (Any) -> bool
    """
    Return True if the object is a string.
    """
    # TODO: 2-3 compatible tests for str, including unicode strings
    return isinstance(x, str)


#: Set of variables defined in the model that might contain C code
C_SYMBOLS = ['Imagnetic', 'Iq', 'Iqxy', 'Iqac', 'Iqabc', 'form_volume', 'c_code']

def _find_source_lines(model_info, kernel_module):
    # type: (ModelInfo, ModuleType) -> None
    """
    Identify the location of the C source inside the model definition file.

    This code runs through the source of the kernel module looking for lines
    that contain C code (because it is a c function definition).  Clearly
    there are all sorts of reasons why this might not work (e.g., code
    commented out in a triple-quoted line block, code built using string
    concatenation, code defined in the branch of an 'if' block, code imported
    from another file), but it should work properly in the 95% case, and for
    the remainder, getting the incorrect line number will merely be
    inconvenient.
    """
    # Only need line numbers if we are creating a C module and the C symbols
    # are defined.
    if (callable(model_info.Iq)
            or not any(hasattr(model_info, s) for s in C_SYMBOLS)):
        return

    # load the module source if we can
    try:
        source = inspect.getsource(kernel_module)
    except IOError:
        return

    # look for symbol at the start of the line
    for lineno, line in enumerate(source.split('\n')):
        for name in C_SYMBOLS:
            if line.startswith(name):
                # Add 1 since some compilers complain about "#line 0"
                model_info.lineno[name] = lineno + 1
                break

def make_model_info(kernel_module):
    # type: (module) -> ModelInfo
    """
    Extract the model definition from the loaded kernel module.

    Fill in default values for parts of the module that are not provided.

    Note: vectorized Iq and Iqac/Iqabc functions will be created for python
    models when the model is first called, not when the model is loaded.
    """
    if hasattr(kernel_module, "model_info"):
        # Custom sum/multi models
        return kernel_module.model_info
    info = ModelInfo()
    #print("make parameter table", kernel_module.parameters)
    parameters = make_parameter_table(getattr(kernel_module, 'parameters', []))
    demo = expand_pars(parameters, getattr(kernel_module, 'demo', None))
    filename = abspath(kernel_module.__file__).replace('.pyc', '.py')
    kernel_id = splitext(basename(filename))[0]
    name = getattr(kernel_module, 'name', None)
    if name is None:
        name = " ".join(w.capitalize() for w in kernel_id.split('_'))

    info.id = kernel_id  # string used to load the kernel
    info.filename = filename
    info.name = name
    info.title = getattr(kernel_module, 'title', name+" model")
    info.description = getattr(kernel_module, 'description', 'no description')
    info.parameters = parameters
    info.demo = demo
    info.composition = None
    info.docs = kernel_module.__doc__
    info.category = getattr(kernel_module, 'category', None)
    info.structure_factor = getattr(kernel_module, 'structure_factor', False)
    info.profile_axes = getattr(kernel_module, 'profile_axes', ['x', 'y'])
    # Note: custom.load_custom_kernel_module assumes the C sources are defined
    # by this attribute.
    info.source = getattr(kernel_module, 'source', [])
    info.c_code = getattr(kernel_module, 'c_code', None)
    # TODO: check the structure of the tests
    info.tests = getattr(kernel_module, 'tests', [])
    info.ER = getattr(kernel_module, 'ER', None) # type: ignore
    info.VR = getattr(kernel_module, 'VR', None) # type: ignore
    info.form_volume = getattr(kernel_module, 'form_volume', None) # type: ignore
    info.Iq = getattr(kernel_module, 'Iq', None) # type: ignore
    info.Iqxy = getattr(kernel_module, 'Iqxy', None) # type: ignore
    info.Iqac = getattr(kernel_module, 'Iqac', None) # type: ignore
    info.Iqabc = getattr(kernel_module, 'Iqabc', None) # type: ignore
    info.Imagnetic = getattr(kernel_module, 'Imagnetic', None) # type: ignore
    info.profile = getattr(kernel_module, 'profile', None) # type: ignore
    info.sesans = getattr(kernel_module, 'sesans', None) # type: ignore
    # Default single and opencl to True for C models.  Python models have callable Iq.
    info.opencl = getattr(kernel_module, 'opencl', not callable(info.Iq))
    info.single = getattr(kernel_module, 'single', not callable(info.Iq))
    info.random = getattr(kernel_module, 'random', None)

    # multiplicity info
    control_pars = [p.id for p in parameters.kernel_parameters if p.is_control]
    default_control = control_pars[0] if control_pars else None
    info.control = getattr(kernel_module, 'control', default_control)
    info.hidden = getattr(kernel_module, 'hidden', None) # type: ignore

    if callable(info.Iq) and parameters.has_2d:
        raise ValueError("oriented python models not supported")

    info.lineno = {}
    _find_source_lines(info, kernel_module)

    return info

class ModelInfo(object):
    """
    Interpret the model definition file, categorizing the parameters.

    The module can be loaded with a normal python import statement if you
    know which module you need, or with __import__('sasmodels.model.'+name)
    if the name is in a string.

    The structure should be mostly static, other than the delayed definition
    of *Iq*, *Iqac* and *Iqabc* if they need to be defined.
    """
    #: Full path to the file defining the kernel, if any.
    filename = None         # type: Optional[str]
    #: Id of the kernel used to load it from the filesystem.
    id = None               # type: str
    #: Display name of the model, which defaults to the model id but with
    #: capitalization of the parts so for example core_shell defaults to
    #: "Core Shell".
    name = None             # type: str
    #: Short description of the model.
    title = None            # type: str
    #: Long description of the model.
    description = None      # type: str
    #: Model parameter table. Parameters are defined using a list of parameter
    #: definitions, each of which is contains parameter name, units,
    #: default value, limits, type and description.  See :class:`Parameter`
    #: for details on the individual parameters.  The parameters are gathered
    #: into a :class:`ParameterTable`, which provides various views into the
    #: parameter list.
    parameters = None       # type: ParameterTable
    #: Demo parameters as a *parameter:value* map used as the default values
    #: for :mod:`compare`.  Any parameters not set in *demo* will use the
    #: defaults from the parameter table.  That means no polydispersity, and
    #: in the case of multiplicity models, a minimal model with no interesting
    #: scattering.
    demo = None             # type: Dict[str, float]
    #: Composition is None if this is an independent model, or it is a
    #: tuple with comoposition type ('product' or 'misture') and a list of
    #: :class:`ModelInfo` blocks for the composed objects.  This allows us
    #: to rebuild a complete mixture or product model from the info block.
    #: *composition* is not given in the model definition file, but instead
    #: arises when the model is constructed using names such as
    #: *sphere*hardsphere* or *cylinder+sphere*.
    composition = None      # type: Optional[Tuple[str, List[ModelInfo]]]
    #: Name of the control parameter for a variant model such as :ref:`rpa`.
    #: The *control* parameter should appear in the parameter table, with
    #: limits defined as *[CASES]*, for case names such as
    #: *CASES = ["diblock copolymer", "triblock copolymer", ...]*.
    #: This should give *limits=[[case1, case2, ...]]*, but the
    #: model loader translates this to *limits=[0, len(CASES)-1]*, and adds
    #: *choices=CASES* to the :class:`Parameter` definition. Note that
    #: models can use a list of cases as a parameter without it being a
    #: control parameter.  Either way, the parameter is sent to the model
    #: evaluator as *float(choice_num)*, where choices are numbered from 0.
    #: See also :attr:`hidden`.
    control = None          # type: str
    #: Different variants require different parameters.  In order to show
    #: just the parameters needed for the variant selected by :attr:`control`,
    #: you should provide a function *hidden(control) -> set(['a', 'b', ...])*
    #: indicating which parameters need to be hidden.  For multiplicity
    #: models, you need to use the complete name of the parameter, including
    #: its number.  So for example, if variant "a" uses only *sld1* and *sld2*,
    #: then *sld3*, *sld4* and *sld5* of multiplicity parameter *sld[5]*
    #: should be in the hidden set.
    hidden = None           # type: Optional[Callable[[int], Set[str]]]
    #: Doc string from the top of the model file.  This should be formatted
    #: using ReStructuredText format, with latex markup in ".. math"
    #: environments, or in dollar signs.  This will be automatically
    #: extracted to a .rst file by :func:`generate.make_docs`, then
    #: converted to HTML or PDF by Sphinx.
    docs = None             # type: str
    #: Location of the model description in the documentation.  This takes the
    #: form of "section" or "section:subsection".  So for example,
    #: :ref:`porod` uses *category="shape-independent"* so it is in the
    #: :ref:`shape-independent` section whereas
    #: :ref:`capped-cylinder` uses: *category="shape:cylinder"*, which puts
    #: it in the :ref:`shape-cylinder` section.
    category = None         # type: Optional[str]
    #: True if the model can be computed accurately with single precision.
    #: This is True by default, but models such as :ref:`bcc-paracrystal` set
    #: it to False because they require double precision calculations.
    single = None           # type: bool
    #: True if the model can be run as an opencl model.  If for some reason
    #: the model cannot be run in opencl (e.g., because the model passes
    #: functions by reference), then set this to false.
    opencl = None           # type: bool
    #: True if the model is a structure factor used to model the interaction
    #: between form factor models.  This will default to False if it is not
    #: provided in the file.
    structure_factor = None # type: bool
    #: List of C source files used to define the model.  The source files
    #: should define the *Iq* function, and possibly *Iqac* or *Iqabc* if the
    #: model defines orientation parameters. Files containing the most basic
    #: functions must appear first in the list, followed by the files that
    #: use those functions.  Form factors are indicated by providing
    #: an :attr:`ER` function.
    source = None           # type: List[str]
    #: The set of tests that must pass.  The format of the tests is described
    #: in :mod:`model_test`.
    tests = None            # type: List[TestCondition]
    #: Returns the effective radius of the model given its volume parameters.
    #: The presence of *ER* indicates that the model is a form factor model
    #: that may be used together with a structure factor to form an implicit
    #: multiplication model.
    #:
    #: The parameters to the *ER* function must be marked with type *volume*.
    #: in the parameter table.  They will appear in the same order as they
    #: do in the table.  The values passed to *ER* will be vectors, with one
    #: value for each polydispersity condition.  For example, if the model
    #: is polydisperse over both length and radius, then both length and
    #: radius will have the same number of values in the vector, with one
    #: value for each *length X radius*.  If only *radius* is polydisperse,
    #: then the value for *length* will be repeated once for each value of
    #: *radius*.  The *ER* function should return one effective radius for
    #: each parameter set.  Multiplicity parameters will be received as
    #: arrays, with one row per polydispersity condition.
    ER = None               # type: Optional[Callable[[np.ndarray], np.ndarray]]
    #: Returns the occupied volume and the total volume for each parameter set.
    #: See :attr:`ER` for details on the parameters.
    VR = None               # type: Optional[Callable[[np.ndarray], Tuple[np.ndarray, np.ndarray]]]
    #: Arbitrary C code containing supporting functions, etc., to be inserted
    #: after everything in source.  This can include Iq and Iqxy functions with
    #: the full function signature, including all parameters.
    c_code = None
    #: Returns the form volume for python-based models.  Form volume is needed
    #: for volume normalization in the polydispersity integral.  If no
    #: parameters are *volume* parameters, then form volume is not needed.
    #: For C-based models, (with :attr:`sources` defined, or with :attr:`Iq`
    #: defined using a string containing C code), form_volume must also be
    #: C code, either defined as a string, or in the sources.
    form_volume = None      # type: Union[None, str, Callable[[np.ndarray], float]]
    #: Returns *I(q, a, b, ...)* for parameters *a*, *b*, etc. defined
    #: by the parameter table.  *Iq* can be defined as a python function, or
    #: as a C function.  If it is defined in C, then set *Iq* to the body of
    #: the C function, including the return statement.  This function takes
    #: values for *q* and each of the parameters as separate *double* values
    #: (which may be converted to float or long double by sasmodels).  All
    #: source code files listed in :attr:`sources` will be loaded before the
    #: *Iq* function is defined.  If *Iq* is not present, then sources should
    #: define *static double Iq(double q, double a, double b, ...)* which
    #: will return *I(q, a, b, ...)*.  Multiplicity parameters are sent as
    #: pointers to doubles.  Constants in floating point expressions should
    #: include the decimal point. See :mod:`generate` for more details.
    Iq = None               # type: Union[None, str, Callable[[np.ndarray], np.ndarray]]
    #: Returns *I(qab, qc, a, b, ...)*.  The interface follows :attr:`Iq`.
    Iqac = None             # type: Union[None, str, Callable[[np.ndarray], np.ndarray]]
    #: Returns *I(qa, qb, qc, a, b, ...)*.  The interface follows :attr:`Iq`.
    Iqabc = None            # type: Union[None, str, Callable[[np.ndarray], np.ndarray]]
    #: Returns *I(qx, qy, a, b, ...)*.  The interface follows :attr:`Iq`.
    Imagnetic = None        # type: Union[None, str, Callable[[np.ndarray], np.ndarray]]
    #: Returns a model profile curve *x, y*.  If *profile* is defined, this
    #: curve will appear in response to the *Show* button in SasView.  Use
    #: :attr:`profile_axes` to set the axis labels.  Note that *y* values
    #: will be scaled by 1e6 before plotting.
    profile = None          # type: Optional[Callable[[np.ndarray], None]]
    #: Axis labels for the :attr:`profile` plot.  The default is *['x', 'y']*.
    #: Only the *x* component is used for now.
    profile_axes = None     # type: Tuple[str, str]
    #: Returns *sesans(z, a, b, ...)* for models which can directly compute
    #: the SESANS correlation function.  Note: not currently implemented.
    sesans = None           # type: Optional[Callable[[np.ndarray], np.ndarray]]
    #: Returns a random parameter set for the model
    random = None           # type: Optional[Callable[[], Dict[str, float]]]
    #: Line numbers for symbols defining C code
    lineno = None           # type: Dict[str, int]

    def __init__(self):
        # type: () -> None
        pass

    def get_hidden_parameters(self, control):
        """
        Returns the set of hidden parameters for the model.  *control* is the
        value of the control parameter.  Note that multiplicity models have
        an implicit control parameter, which is the parameter that controls
        the multiplicity.
        """
        if self.hidden is not None:
            hidden = self.hidden(control)
        else:
            controls = [p for p in self.parameters.kernel_parameters
                        if p.is_control]
            if len(controls) != 1:
                raise ValueError("more than one control parameter")
            hidden = set(p.id+str(k)
                         for p in self.parameters.kernel_parameters
                         for k in range(control+1, p.length+1)
                         if p.length > 1)
            for p in self.parameters.kernel_parameters:
                if p.length > 1 and p.type == "sld":
                    for k in range(control+1, p.length+1):
                        base = p.id+str(k)
                        hidden.update((base+"_M0", base+"_mtheta", base+"_mphi"))
        return hidden
