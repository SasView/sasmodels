
import numpy as np

# TODO: turn ModelInfo into a proper class
ModelInfo = dict

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
    set_vector_length_from_reference(partable)
    return partable

def set_vector_length_from_reference(partable):
    # Sort out the length of the vector parameters such as thickness[n]
    for p in partable:
        if p.length_control:
           ref = partable[p.length_control]
           low, high = ref.limits
           if int(low) != low or int(high) != high or low<0 or high>20:
               raise ValueError("expected limits on %s to be within [0, 20]"%ref.name)
           p.length = low

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
        limits = [1, len(choices)]
    else:
        choices = None
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
    *polydisperse* is true if the parameter accepts a polydispersity
    *relative_pd* is true if that polydispersity is relative

    In the usual process these values are set by :func:`make_parameter_table`
    and :func:`parse_parameter` therein.
    """
    def __init__(self, name, units='', default=None, limits=(-np.inf, np.inf),
                 type='', description=''):
        self.id = name
        self.name = name
        self.default = default
        self.limits = limits
        self.type = type
        self.description = description
        self.choices = None

        # Length and length_control will be filled in by
        # set_vector_length_from_reference(partable) once the complete
        # parameter table is available.
        self.length = 1
        self.length_control = None

        # TODO: need better control over whether a parameter is polydisperse
        self.polydisperse = False
        self.relative_pd = False

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
    # scale and background are implicit parameters
    COMMON = [Parameter(*p) for p in COMMON_PARAMETERS]

    def __init__(self, parameters):
        self.parameters = self.COMMON + parameters
        self._name_table= dict((p.name, p) for p in parameters)
        self._categorize_parameters()

    def __getitem__(self, k):
        if isinstance(k, (int, slice)):
            return self.parameters[k]
        else:
            return self._name_table[k]

    def __contains__(self, key):
        return key in self._name_table

    def __iter__(self):
        return iter(self.parameters)

    def kernel_pars(self, ptype=None):
        """
        Return the parameters to the user kernel which match the given type.

        Types include '1d' for Iq kernels, '2d' for Iqxy kernels and
        'volume' for form_volume kernels.
        """
        # Assumes background and scale are the first two parameters
        if ptype is None:
            return self.parameters[2:]
        else:
            return [p for p in self.parameters[2:] if p in self.type[ptype]]

    def _categorize_parameters(self):
        """
        Build parameter categories out of the the parameter definitions.

        Returns a dictionary of categories.

        Note: these categories are subject to change, depending on the needs of
        the UI and the needs of the kernel calling function.

        The categories are as follows:

        * *volume* list of volume parameter names
        * *orientation* list of orientation parameters
        * *magnetic* list of magnetic parameters
        * *sld* list of parameters that have no type info
        * *other* list of parameters that have no type info

        Each parameter is in one and only one category.
        """
        pars = self.parameters

        par_type = {
            'volume': [], 'orientation': [], 'magnetic': [], 'sld': [], 'other': [],
        }
        for p in self.parameters:
            par_type[p.type if p.type else 'other'].append(p)
        par_type['1d'] = [p for p in pars if p.type not in ('orientation', 'magnetic')]
        par_type['2d'] = [p for p in pars if p.type != 'magnetic']
        par_type['magnetic'] = [p for p in pars]
        par_type['pd'] = [p for p in pars if p.polydisperse]
        par_type['pd_relative'] = [p for p in pars if p.relative_pd]
        self.type = par_type

        # find index of theta (or whatever variable is used for spherical
        # normalization during polydispersity...
        if 'theta' in par_type['2d']:
            # TODO: may be an off-by 2 bug due to background and scale
            # TODO: is theta always the polar coordinate?
            self.theta_par = [k for k,p in enumerate(pars) if p.name=='theta'][0]
        else:
            self.theta_par = -1

    @property
    def defaults(self):
        return dict((p.name, p.default) for p in self.parameters)

    @property
    def num_pd(self):
        """
        Number of distributional parameters in the model (polydispersity in
        shape dimensions and orientational distributions).
        """
        return len(self.type['pd'])

    @property
    def has_2d(self):
        return self.type['orientation'] or self.type['magnetic']


def isstr(x):
    # TODO: 2-3 compatible tests for str, including unicode strings
    return isinstance(x, str)

