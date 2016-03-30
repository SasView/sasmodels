
import numpy as np

# TODO: turn ModelInfo into a proper class
ModelInfo = dict

MAX_PD = 4

COMMON_PARAMETERS = [
    ["scale", "", 1, [0, np.inf], "", "Source intensity"],
    ["background", "1/cm", 1e-3, [0, np.inf], "", "Source background"],
]
assert (len(COMMON_PARAMETERS) == 2
        and COMMON_PARAMETERS[0][0]=="scale"
        and COMMON_PARAMETERS[1][0]=="background"), "don't change common parameters"
# assumptions about common parameters exist throughout the code, such as:
# (1) kernel functions Iq, Iqxy, form_volume, ... don't see them
# (2) kernel drivers assume scale is par[0] and background is par[1]
# (3) mixture models drop the background on components and replace the scale
#     with a scale that varies from [-inf, inf]
# (4) product models drop the background and reassign scale
# and maybe other places.
# Note that scale and background cannot be coordinated parameters whose value
# depends on the some polydisperse parameter with the current implementation

def make_parameter_table(pars):
    processed = []
    for p in pars:
        if not isinstance(p, list) or len(p) != 6:
            raise ValueError("Parameter should be [name, units, default, limits, type, desc], but got %r"
                             %str(p))
        processed.append(parse_parameter(*p))
    partable = ParameterTable(processed)
    return partable

def parse_parameter(name, units='', default=None,
                    limits=(-np.inf, np.inf), type='', description=''):
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
    except:
        raise ValueError("invalid limits %s for %s"%(limits, name))

    if not isinstance(default, (int, float)):
        raise ValueError("expected default %r to be a number for %s"
                         % (default, name))
    if default < low or default > high:
        raise ValueError("default value %r not in range for %s"
                         % (default, name))

    if type not in ("volume", "orientation", "sld", "magnetic", ""):
        raise ValueError("unexpected type %r for %s" % (type, name))

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
    if type=='' and (pid.startswith('sld') or pid.endswith('sld')):
        type = 'sld'

    # Check if using a vector definition, name[k], as the parameter name
    if ref:
        if ref == '':
            raise ValueError("Need to specify vector length for %s"%name)
        try:
            length = int(ref)
            control = None
        except:
            length = None
            control = ref
    else:
        length = 1
        control = None

    # Build the parameter
    parameter = Parameter(name=name, units=units, default=default,
                          limits=limits, type=type, description=description)

    # TODO: need better control over whether a parameter is polydisperse
    parameter.polydisperse = type in ('orientation', 'volume')
    parameter.relative_pd = type in ('volume')
    parameter.choices = choices
    parameter.length = length
    parameter.length_control = control

    return parameter


def set_demo(model_info, demo):
    """
    Assign demo parameters to model_info['demo']

    If demo is not defined, then use defaults.

    If demo is defined, override defaults with value from demo.

    If demo references vector fields, such as thickness[n], then support
    different ways of assigning the demo values, including assigning a
    specific value (e.g., thickness3=50.0), assigning a new value to all
    (e.g., thickness=50.0) or assigning values using list notation.
    """
    partable = model_info['parameters']
    if demo is None:
        result = partable.defaults
    else:
        pars = dict((p.id, p) for p in partable.kernel_parameters)
        result = partable.defaults.copy()
        vectors = dict((name,value) for name,value in demo.items()
                       if name in pars and pars[name].length > 1)
        if vectors:
            scalars = dict((name, value) for name, value in demo.items()
                           if name not in pars or pars[name].length == 1)
            for name, value in vectors.items():
                if np.isscalar(value):
                    # support for the form
                    #    demo(thickness=0, thickness2=50)
                    for k in pars[name].length:
                        key = name+str(k)
                        if key not in scalars:
                            scalars[key] = vectors
                else:
                    # supoprt for the form
                    #    demo(thickness=[20,10,3])
                    for (k,v) in enumerate(value):
                        scalars[name+str(k)] = v
            result.update(scalars)
        else:
            result.update(demo)

    model_info['demo'] = result

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
    other sld parameters.

    *description* is a short description of the parameter.  This will
    be displayed in the parameter table and used as a tool tip for the
    parameter value in the user interface.

    Additional values can be set after the parameter is created:

    *length* is the length of the field if it is a vector field
    *length_control* is the parameter which sets the vector length
    *is_control* is True if the parameter is a control parameter for a vector
    *polydisperse* is true if the parameter accepts a polydispersity
    *relative_pd* is true if that polydispersity is relative

    In the usual process these values are set by :func:`make_parameter_table`
    and :func:`parse_parameter` therein.
    """
    def __init__(self, name, units='', default=None, limits=(-np.inf, np.inf),
                 type='', description=''):
        self.id = name.split('[')[0].strip()
        self.name = name
        self.units = units
        self.default = default
        self.limits = limits
        self.type = type
        self.description = description

        # Length and length_control will be filled in once the complete
        # parameter table is available.
        self.length = 1
        self.length_control = None
        self.is_control = False

        # TODO: need better control over whether a parameter is polydisperse
        self.polydisperse = False
        self.relative_pd = False

        # choices are also set externally.
        self.choices = []

    def as_definition(self):
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
        # Note: if the parameter is a struct type, then we will need to use
        # &prefix+id.  For scalars and vectors we can just use prefix+id.
        return prefix + self.id

    def __str__(self):
        return "<%s>"%self.name

    def __repr__(self):
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
        self.kernel_parameters = parameters
        self._set_vector_lengths()
        self._make_call_parameter_list()
        self._categorize_parameters()
        self._set_defaults()
        #self._name_table= dict((p.id, p) for p in parameters)

    def _set_vector_lengths(self):
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
                if int(low) != low or int(high) != high or low<0 or high>20:
                    raise ValueError("expected limits on %s to be within [0, 20]"
                                     % ref.name)
                p.length = high

    def _set_defaults(self):
        # Construct default values, including vector defaults
        defaults = {}
        for p in self.call_parameters:
            if p.length == 1:
                defaults[p.id] = p.default
            else:
                for k in range(1, p.length+1):
                    defaults["%s%d"%(p.id, k)] = p.default
        self.defaults = defaults

    def _make_call_parameter_list(self):
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
        self.call_parameters = full_list

    """ # Suppress these for now until we see how they are used
    def __getitem__(self, k):
        if isinstance(k, (int, slice)):
            return self.parameters[k]
        else:
            return self._name_table[k]

    def __contains__(self, key):
        return key in self._name_table

    def __iter__(self):
        return iter(self.expanded_parameters)
    """

    def _categorize_parameters(self):
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

    def user_parameters(self, pars={}, is2d=True):
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
    # TODO: 2-3 compatible tests for str, including unicode strings
    return isinstance(x, str)

