"""
SAS model constructor.

Small angle scattering models are defined by a set of kernel functions:

    *Iq(q, p1, p2, ...)* returns the scattering at q for a form with
    particular dimensions averaged over all orientations.

    *Iqxy(qx, qy, p1, p2, ...)* returns the scattering at qx, qy for a form
    with particular dimensions for a single orientation.

    *Imagnetic(qx, qy, result[], p1, p2, ...)* returns the scattering for the
    polarized neutron spin states (up-up, up-down, down-up, down-down) for
    a form with particular dimensions for a single orientation.

    *form_volume(p1, p2, ...)* returns the volume of the form with particular
    dimension.

    *ER(p1, p2, ...)* returns the effective radius of the form with
    particular dimensions.

    *VR(p1, p2, ...)* returns the volume ratio for core-shell style forms.

    #define INVALID(v) (expr)  returns False if v.parameter is invalid
    for some parameter or other (e.g., v.bell_radius < v.radius).  If
    necessary, the expression can call a function.

These functions are defined in a kernel module .py script and an associated
set of .c files.  The model constructor will use them to create models with
polydispersity across volume and orientation parameters, and provide
scale and background parameters for each model.

*Iq*, *Iqxy*, *Imagnetic* and *form_volume* should be stylized C-99
functions written for OpenCL.  All functions need prototype declarations
even if the are defined before they are used.  OpenCL does not support
*#include* preprocessor directives, so instead the list of includes needs
to be given as part of the metadata in the kernel module definition.
The included files should be listed using a path relative to the kernel
module, or if using "lib/file.c" if it is one of the standard includes
provided with the sasmodels source.  The includes need to be listed in
order so that functions are defined before they are used.

Floating point values should be declared as *double*.  For single precision
calculations, *double* will be replaced by *float*.  The single precision
conversion will also tag floating point constants with "f" to make them
single precision constants.  When using integral values in floating point
expressions, they should be expressed as floating point values by including
a decimal point.  This includes 0., 1. and 2.

OpenCL has a *sincos* function which can improve performance when both
the *sin* and *cos* values are needed for a particular argument.  Since
this function does not exist in C99, all use of *sincos* should be
replaced by the macro *SINCOS(value, sn, cn)* where *sn* and *cn* are
previously declared *double* variables.  When compiled for systems without
OpenCL, *SINCOS* will be replaced by *sin* and *cos* calls.   If *value* is
an expression, it will appear twice in this case; whether or not it will be
evaluated twice depends on the quality of the compiler.

If the input parameters are invalid, the scattering calculator should
return a negative number. Particularly with polydispersity, there are
some sets of shape parameters which lead to nonsensical forms, such
as a capped cylinder where the cap radius is smaller than the
cylinder radius.  The polydispersity calculation will ignore these points,
effectively chopping the parameter weight distributions at the boundary
of the infeasible region.  The resulting scattering will be set to
background.  This will work correctly even when polydispersity is off.

*ER* and *VR* are python functions which operate on parameter vectors.
The constructor code will generate the necessary vectors for computing
them with the desired polydispersity.
The kernel module must set variables defining the kernel meta data:

    *id* is an implicit variable formed from the filename.  It will be
    a valid python identifier, and will be used as the reference into
    the html documentation, with '_' replaced by '-'.

    *name* is the model name as displayed to the user.  If it is missing,
    it will be constructed from the id.

    *title* is a short description of the model, suitable for a tool tip,
    or a one line model summary in a table of models.

    *description* is an extended description of the model to be displayed
    while the model parameters are being edited.

    *parameters* is the list of parameters.  Parameters in the kernel
    functions must appear in the same order as they appear in the
    parameters list.  Two additional parameters, *scale* and *background*
    are added to the beginning of the parameter list.  They will show up
    in the documentation as model parameters, but they are never sent to
    the kernel functions.  Note that *effect_radius* and *volfraction*
    must occur first in structure factor calculations.

    *category* is the default category for the model.  The category is
    two level structure, with the form "group:section", indicating where
    in the manual the model will be located.  Models are alphabetical
    within their section.

    *source* is the list of C-99 source files that must be joined to
    create the OpenCL kernel functions.  The files defining the functions
    need to be listed before the files which use the functions.

    *ER* is a python function defining the effective radius.  If it is
    not present, the effective radius is 0.

    *VR* is a python function defining the volume ratio.  If it is not
    present, the volume ratio is 1.

    *form_volume*, *Iq*, *Iqxy*, *Imagnetic* are strings containing the
    C source code for the body of the volume, Iq, and Iqxy functions
    respectively.  These can also be defined in the last source file.

    *Iq* and *Iqxy* also be instead be python functions defining the
    kernel.  If they are marked as *Iq.vectorized = True* then the
    kernel is passed the entire *q* vector at once, otherwise it is
    passed values one *q* at a time.  The performance improvement of
    this step is significant.

    *demo* is a dictionary of parameter=value defining a set of
    parameters to use by default when *compare* is called.  Any
    parameter not set in *demo* gets the initial value from the
    parameter list.  *demo* is mostly needed to set the default
    polydispersity values for tests.

    *oldname* is the name of the model in sasview before sasmodels
    was split into its own package, and *oldpars* is a dictionary
    of *parameter: old_parameter* pairs defining the new names for
    the parameters.  This is used by *compare* to check the values
    of the new model against the values of the old model before
    you are ready to add the new model to sasmodels.


An *model_info* dictionary is constructed from the kernel meta data and
returned to the caller.

The model evaluator, function call sequence consists of q inputs and the return vector,
followed by the loop value/weight vector, followed by the values for
the non-polydisperse parameters, followed by the lengths of the
polydispersity loops.  To construct the call for 1D models, the
categories *fixed-1d* and *pd-1d* list the names of the parameters
of the non-polydisperse and the polydisperse parameters respectively.
Similarly, *fixed-2d* and *pd-2d* provide parameter names for 2D models.
The *pd-rel* category is a set of those parameters which give
polydispersitiy as a portion of the value (so a 10% length dispersity
would use a polydispersity value of 0.1) rather than absolute
dispersity such as an angle plus or minus 15 degrees.

The *volume* category lists the volume parameters in order for calls
to volume within the kernel (used for volume normalization) and for
calls to ER and VR for effective radius and volume ratio respectively.

The *orientation* and *magnetic* categories list the orientation and
magnetic parameters.  These are used by the sasview interface.  The
blank category is for parameters such as scale which don't have any
other marking.

The doc string at the start of the kernel module will be used to
construct the model documentation web pages.  Embedded figures should
appear in the subdirectory "img" beside the model definition, and tagged
with the kernel module name to avoid collision with other models.  Some
file systems are case-sensitive, so only use lower case characters for
file names and extensions.


The function :func:`make` loads the metadata from the module and returns
the kernel source.  The function :func:`make_doc` extracts the doc string
and adds the parameter table to the top.  The function :func:`model_sources`
returns a list of files required by the model.

Code follows the C99 standard with the following extensions and conditions::

    M_PI_180 = pi/180
    M_4PI_3 = 4pi/3
    square(x) = x*x
    cube(x) = x*x*x
    sinc(x) = sin(x)/x, with sin(0)/0 -> 1
    all double precision constants must include the decimal point
    all double declarations may be converted to half, float, or long double
    FLOAT_SIZE is the number of bytes in the converted variables
"""
from __future__ import print_function

#TODO: determine which functions are useful outside of generate
#__all__ = ["model_info", "make_doc", "make_source", "convert_type"]

import sys
from os.path import abspath, dirname, join as joinpath, exists, basename, \
    splitext, getmtime
import re
import string
import warnings

import numpy as np

from .modelinfo import ModelInfo, Parameter, make_parameter_table

# TODO: identify model files which have changed since loading and reload them.

TEMPLATE_ROOT = dirname(__file__)

MAX_PD = 4

F16 = np.dtype('float16')
F32 = np.dtype('float32')
F64 = np.dtype('float64')
try:  # CRUFT: older numpy does not support float128
    F128 = np.dtype('float128')
except TypeError:
    F128 = None

# Conversion from units defined in the parameter table for each model
# to units displayed in the sphinx documentation.
RST_UNITS = {
    "Ang": "|Ang|",
    "1/Ang": "|Ang^-1|",
    "1/Ang^2": "|Ang^-2|",
    "1e-6/Ang^2": "|1e-6Ang^-2|",
    "degrees": "degree",
    "1/cm": "|cm^-1|",
    "Ang/cm": "|Ang*cm^-1|",
    "g/cm3": "|g/cm^3|",
    "mg/m2": "|mg/m^2|",
    "": "None",
    }

# Headers for the parameters tables in th sphinx documentation
PARTABLE_HEADERS = [
    "Parameter",
    "Description",
    "Units",
    "Default value",
    ]

# Minimum width for a default value (this is shorter than the column header
# width, so will be ignored).
PARTABLE_VALUE_WIDTH = 10

# Documentation header for the module, giving the model name, its short
# description and its parameter table.  The remainder of the doc comes
# from the module docstring.
DOC_HEADER = """.. _%(id)s:

%(name)s
=======================================================

%(title)s

%(parameters)s

%(returns)s

%(docs)s
"""

def format_units(units):
    """
    Convert units into ReStructured Text format.
    """
    return "string" if isinstance(units, list) else RST_UNITS.get(units, units)

def make_partable(pars):
    """
    Generate the parameter table to include in the sphinx documentation.
    """
    column_widths = [
        max(len(p.name) for p in pars),
        max(len(p.description) for p in pars),
        max(len(format_units(p.units)) for p in pars),
        PARTABLE_VALUE_WIDTH,
        ]
    column_widths = [max(w, len(h))
                     for w, h in zip(column_widths, PARTABLE_HEADERS)]

    sep = " ".join("="*w for w in column_widths)
    lines = [
        sep,
        " ".join("%-*s" % (w, h)
                 for w, h in zip(column_widths, PARTABLE_HEADERS)),
        sep,
        ]
    for p in pars:
        lines.append(" ".join([
            "%-*s" % (column_widths[0], p.name),
            "%-*s" % (column_widths[1], p.description),
            "%-*s" % (column_widths[2], format_units(p.units)),
            "%*g" % (column_widths[3], p.default),
            ]))
    lines.append(sep)
    return "\n".join(lines)

def _search(search_path, filename):
    """
    Find *filename* in *search_path*.

    Raises ValueError if file does not exist.
    """
    for path in search_path:
        target = joinpath(path, filename)
        if exists(target):
            return target
    raise ValueError("%r not found in %s" % (filename, search_path))


def model_sources(model_info):
    """
    Return a list of the sources file paths for the module.
    """
    search_path = [dirname(model_info['filename']),
                   abspath(joinpath(dirname(__file__), 'models'))]
    return [_search(search_path, f) for f in model_info['source']]

def timestamp(model_info):
    """
    Return a timestamp for the model corresponding to the most recently
    changed file or dependency.
    """
    source_files = (model_sources(model_info)
                    + model_templates()
                    + [model_info['filename']])
    newest = max(getmtime(f) for f in source_files)
    return newest

def convert_type(source, dtype):
    """
    Convert code from double precision to the desired type.

    Floating point constants are tagged with 'f' for single precision or 'L'
    for long double precision.
    """
    if dtype == F16:
        fbytes = 2
        source = _convert_type(source, "float", "f")
    elif dtype == F32:
        fbytes = 4
        source = _convert_type(source, "float", "f")
    elif dtype == F64:
        fbytes = 8
        # no need to convert if it is already double
    elif dtype == F128:
        fbytes = 16
        source = _convert_type(source, "long double", "L")
    else:
        raise ValueError("Unexpected dtype in source conversion: %s"%dtype)
    return ("#define FLOAT_SIZE %d\n"%fbytes)+source


def _convert_type(source, type_name, constant_flag):
    """
    Replace 'double' with *type_name* in *source*, tagging floating point
    constants with *constant_flag*.
    """
    # Convert double keyword to float/long double/half.
    # Accept an 'n' # parameter for vector # values, where n is 2, 4, 8 or 16.
    # Assume complex numbers are represented as cdouble which is typedef'd
    # to double2.
    source = re.sub(r'(^|[^a-zA-Z0-9_]c?)double(([248]|16)?($|[^a-zA-Z0-9_]))',
                    r'\1%s\2'%type_name, source)
    # Convert floating point constants to single by adding 'f' to the end,
    # or long double with an 'L' suffix.  OS/X complains if you don't do this.
    source = re.sub(r'[^a-zA-Z_](\d*[.]\d+|\d+[.]\d*)([eE][+-]?\d+)?',
                    r'\g<0>%s'%constant_flag, source)
    return source


def kernel_name(model_info, is_2d):
    """
    Name of the exported kernel symbol.
    """
    return model_info['name'] + "_" + ("Iqxy" if is_2d else "Iq")


def indent(s, depth):
    """
    Indent a string of text with *depth* additional spaces on each line.
    """
    spaces = " "*depth
    sep = "\n" + spaces
    return spaces + sep.join(s.split("\n"))


_template_cache = {}
def load_template(filename):
    path = joinpath(TEMPLATE_ROOT, filename)
    mtime = getmtime(path)
    if filename not in _template_cache or mtime > _template_cache[filename][0]:
        with open(path) as fid:
            _template_cache[filename] = (mtime, fid.read(), path)
    return _template_cache[filename][1]

def model_templates():
    # TODO: fails DRY; templates are listed in two places.
    # should instead have model_info contain a list of paths
    return [joinpath(TEMPLATE_ROOT, filename)
            for filename in ('kernel_header.c', 'kernel_iq.c')]


_FN_TEMPLATE = """\
double %(name)s(%(pars)s);
double %(name)s(%(pars)s) {
    %(body)s
}


"""

def _gen_fn(name, pars, body):
    """
    Generate a function given pars and body.

    Returns the following string::

         double fn(double a, double b, ...);
         double fn(double a, double b, ...) {
             ....
         }
    """
    par_decl = ', '.join(p.as_argument() for p in pars) if pars else 'void'
    return _FN_TEMPLATE % {'name': name, 'body': body, 'pars': par_decl}

def _call_pars(prefix, pars):
    """
    Return a list of *prefix.parameter* from parameter items.
    """
    return [p.as_call_reference(prefix) for p in pars]

_IQXY_PATTERN = re.compile("^((inline|static) )? *(double )? *Iqxy *([(]|$)",
                           flags=re.MULTILINE)
def _have_Iqxy(sources):
    """
    Return true if any file defines Iqxy.

    Note this is not a C parser, and so can be easily confused by
    non-standard syntax.  Also, it will incorrectly identify the following
    as having Iqxy::

        /*
        double Iqxy(qx, qy, ...) { ... fill this in later ... }
        */

    If you want to comment out an Iqxy function, use // on the front of the
    line instead.
    """
    for code in sources:
        if _IQXY_PATTERN.search(code):
            return True
    else:
        return False

def make_source(model_info):
    """
    Generate the OpenCL/ctypes kernel from the module info.

    Uses source files found in the given search path.
    """
    if callable(model_info['Iq']):
        return None

    # TODO: need something other than volume to indicate dispersion parameters
    # No volume normalization despite having a volume parameter.
    # Thickness is labelled a volume in order to trigger polydispersity.
    # May want a separate dispersion flag, or perhaps a separate category for
    # disperse, but not volume.  Volume parameters also use relative values
    # for the distribution rather than the absolute values used by angular
    # dispersion.  Need to be careful that necessary parameters are available
    # for computing volume even if we allow non-disperse volume parameters.

    partable = model_info['parameters']

    # Identify parameters for Iq, Iqxy, Iq_magnetic and form_volume.
    # Note that scale and volume are not possible types.

    # Load templates and user code
    kernel_header = load_template('kernel_header.c')
    kernel_code = load_template('kernel_iq.c')
    user_code = [open(f).read() for f in model_sources(model_info)]

    # Build initial sources
    source = [kernel_header] + user_code

    vol_parameters = partable.kernel_pars('volume')
    iq_parameters = partable.kernel_pars('1d')
    iqxy_parameters = partable.kernel_pars('2d')

    # Make parameters for q, qx, qy so that we can use them in declarations
    q, qx, qy = [Parameter(name=v) for v in ('q', 'qx', 'qy')]
    # Generate form_volume function, etc. from body only
    if model_info['form_volume'] is not None:
        pars = vol_parameters
        source.append(_gen_fn('form_volume', pars, model_info['form_volume']))
    if model_info['Iq'] is not None:
        pars = [q] + iq_parameters
        source.append(_gen_fn('Iq', pars, model_info['Iq']))
    if model_info['Iqxy'] is not None:
        pars = [qx, qy] + iqxy_parameters
        source.append(_gen_fn('Iqxy', pars, model_info['Iqxy']))

    # Define the parameter table
    source.append("#define PARAMETER_TABLE \\")
    source.append("\\\n".join(p.as_definition()
                                  for p in model_info['parameters'][2:]))

    # Define the function calls
    if vol_parameters:
        refs = _call_pars("v.", vol_parameters)
        call_volume = "#define CALL_VOLUME(v) form_volume(%s)" % (",".join(refs))
    else:
        # Model doesn't have volume.  We could make the kernel run a little
        # faster by not using/transferring the volume normalizations, but
        # the ifdef's reduce readability more than is worthwhile.
        call_volume = "#define CALL_VOLUME(v) 0.0"
    source.append(call_volume)

    refs = ["q[i]"] + _call_pars("v.", iq_parameters)
    call_iq = "#define CALL_IQ(q,i,v) Iq(%s)" % (",".join(refs))
    if _have_Iqxy(user_code):
        # Call 2D model
        refs = ["q[2*i]", "q[2*i+1]"] + _call_pars("v.", iqxy_parameters)
        call_iqxy = "#define CALL_IQ(q,i,v) Iqxy(%s)" % (",".join(refs))
    else:
        # Call 1D model with sqrt(qx^2 + qy^2)
        warnings.warn("Creating Iqxy = Iq(sqrt(qx^2 + qy^2))")
        # still defined:: refs = ["q[i]"] + _call_pars("v", iq_parameters)
        pars_sqrt = ["sqrt(q[2*i]*q[2*i]+q[2*i+1]*q[2*i+1])"] + refs[1:]
        call_iqxy = "#define CALL_IQ(q,i,v) Iq(%s)" % (",".join(pars_sqrt))

    # Fill in definitions for numbers of parameters
    source.append("#define MAX_PD %s"%model_info['max_pd'])
    source.append("#define NPARS %d"%(len(partable.kernel_pars())))

    # TODO: allow mixed python/opencl kernels?

    # define the Iq kernel
    source.append("#define KERNEL_NAME %s_Iq"%model_info['name'])
    source.append(call_iq)
    source.append(kernel_code)
    source.append("#undef CALL_IQ")
    source.append("#undef KERNEL_NAME")

    # define the Iqxy kernel from the same source with different #defines
    source.append("#define KERNEL_NAME %s_Iqxy"%model_info['name'])
    source.append(call_iqxy)
    source.append(kernel_code)
    source.append("#undef CALL_IQ")
    source.append("#undef KERNEL_NAME")

    return '\n'.join(source)

def categorize_parameters(pars):
    """
    Categorize the parameters by use:

    * *pd* list of polydisperse parameters in order; gui should test whether
      they are in *2d* or *magnetic* as appropriate for the data
    * *1d* set of parameters that are used to compute 1D patterns
    * *2d* set of parameters that are used to compute 2D patterns (which
      includes all 1D parameters)
    * *magnetic* set of parameters that are used to compute magnetic
      patterns (which includes all 1D and 2D parameters)
    * *pd_relative* is the set of parameters with relative distribution
      width (e.g., radius +/- 10%) rather than absolute distribution
      width (e.g., theta +/- 6 degrees).
    * *theta_par* is the index of the polar angle polydispersion parameter
      or -1 if no such parameter exists
    """
    par_set = {}

def process_parameters(model_info):
    """
    Process parameter block, precalculating parameter details.
    """
    partable = model_info['parameters']
    if model_info.get('demo', None) is None:
        model_info['demo'] = partable.defaults

    # Don't use more polydisperse parameters than are available in the model
    # Note: we can do polydispersity on arbitrary parameters, so it is not
    # clear that this is a good idea; it does however make the poly_details
    # code easier to write, so we will leave it in for now.
    model_info['max_pd'] = min(partable.num_pd, MAX_PD)

def mono_details(model_info):
    # TODO: move max_pd into ParameterTable?
    max_pd = model_info['max_pd']
    pars = model_info['parameters'].kernel_pars()
    npars = len(pars)
    par_offset = 5*max_pd
    constants_offset = par_offset + 3*npars

    details = np.zeros(constants_offset + 2, 'int32')
    details[0*max_pd:1*max_pd] = range(max_pd)       # pd_par: arbitrary order; use first
    details[1*max_pd:2*max_pd] = [1]*max_pd          # pd_length: only one element
    details[2*max_pd:3*max_pd] = range(max_pd)       # pd_offset: consecutive 1.0 weights
    details[3*max_pd:4*max_pd] = [1]*max_pd          # pd_stride: vectors of length 1
    details[4*max_pd:5*max_pd] = [0]*max_pd          # pd_isvol: doens't matter if no norm
    details[par_offset+0*npars:par_offset+1*npars] = range(2, npars+2) # par_offset: skip scale and background
    details[par_offset+1*npars:par_offset+2*npars] = [0]*npars         # no coordination
    #details[p+npars] = 1 # par_coord[0] is coordinated with the first par?
    details[par_offset+2*npars:par_offset+3*npars] = 0 # fast coord with 0
    details[constants_offset]   = 1     # fast_coord_count: one fast index
    details[constants_offset+1] = -1    # theta_par: None
    return details

def poly_details(model_info, weights):
    weights = weights[2:]

    # TODO: move max_pd into ParameterTable?
    max_pd = model_info['max_pd']
    pars = model_info['parameters'].kernel_pars
    npars = len(pars)
    par_offset = 5*max_pd
    constants_offset = par_offset + 3*npars

    # Decreasing list of polydispersity lengths
    # Note: the reversing view, x[::-1], does not require a copy
    pd_length = np.array([len(w) for w in weights])
    print (pd_length)
    print (weights)
    pd_offset = np.cumsum(np.hstack((0, pd_length)))
    pd_isvol = np.array([p.type=='volume' for p in pars])
    idx = np.argsort(pd_length)[::-1][:max_pd]
    print (idx)
    pd_stride = np.cumprod(np.hstack((1, pd_length[idx][:-1])))
    par_offsets = np.cumsum(np.hstack((2, pd_length)))[:-1]

    theta_par = -1
    if 'theta_par' in model_info:
        theta_par = model_info['theta_par']
        if theta_par >= 0 and pd_length[theta_par] <= 1:
            theta_par = -1

    details = np.empty(constants_offset + 2, 'int32')
    details[0*max_pd:1*max_pd] = idx             # pd_par
    details[1*max_pd:2*max_pd] = pd_length[idx]
    details[2*max_pd:3*max_pd] = pd_offset[idx]
    details[3*max_pd:4*max_pd] = pd_stride
    details[4*max_pd:5*max_pd] = pd_isvol[idx]
    details[par_offset+0*npars:par_offset+1*npars] = par_offsets
    details[par_offset+1*npars:par_offset+2*npars] = 0  # no coordination for most
    details[par_offset+2*npars:par_offset+3*npars] = 0  # no fast coord with 0
    coord_offset = par_offset+1*npars
    for k,parameter_num in enumerate(idx):
        details[coord_offset+parameter_num] = 2**k
    details[constants_offset] = 1   # fast_coord_count: one fast index
    details[constants_offset+1] = theta_par
    print ("details",details)
    return details

def constrained_poly_details(model_info, weights, constraints):
    # Need to find the independently varying pars and sort them
    # Need to build a coordination list for the dependent variables
    # Need to generate a constraints function which takes values
    # and weights, returning par blocks
    raise NotImplementedError("Can't handle constraints yet")


def create_default_functions(model_info):
    """
    Autogenerate missing functions, such as Iqxy from Iq.

    This only works for Iqxy when Iq is written in python. :func:`make_source`
    performs a similar role for Iq written in C.
    """
    if model_info['Iq'] is not None and model_info['Iqxy'] is None:
        partable = model_info['parameters']
        if partable.type['1d'] != partable.type['2d']:
            raise ValueError("Iqxy model is missing")
        Iq = model_info['Iq']
        def Iqxy(qx, qy, **kw):
            return Iq(np.sqrt(qx**2 + qy**2), **kw)
        model_info['Iqxy'] = Iqxy


def make_model_info(kernel_module):
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
    * *oldname* is the model name in pre-4.0 Sasview
    * *oldpars* is the *{new: old}* parameter translation table
      from pre-4.0 Sasview
    * *composition* is None if the model is independent, otherwise it is a
      tuple with composition type ('product' or 'mixture') and a list of
      *model_info* blocks for the composition objects.  This allows us to
      build complete product and mixture models from just the info.
    * *max_pd* is the max polydispersity dimension.  This is constant and
      should not be reset.  You may be able to change it when the program
      starts by setting *sasmodels.generate.MAX_PD*.

    """
    # TODO: maybe turn model_info into a class ModelDefinition
    #print("make parameter table", kernel_module.parameters)
    parameters = make_parameter_table(kernel_module.parameters)
    filename = abspath(kernel_module.__file__)
    kernel_id = splitext(basename(filename))[0]
    name = getattr(kernel_module, 'name', None)
    if name is None:
        name = " ".join(w.capitalize() for w in kernel_id.split('_'))
    model_info = dict(
        id=kernel_id,  # string used to load the kernel
        filename=abspath(kernel_module.__file__),
        name=name,
        title=kernel_module.title,
        description=kernel_module.description,
        parameters=parameters,
        composition=None,
        docs=kernel_module.__doc__,
        category=getattr(kernel_module, 'category', None),
        single=getattr(kernel_module, 'single', True),
        structure_factor=getattr(kernel_module, 'structure_factor', False),
        variant_info=getattr(kernel_module, 'invariant_info', None),
        demo=getattr(kernel_module, 'demo', None),
        source=getattr(kernel_module, 'source', []),
        oldname=getattr(kernel_module, 'oldname', None),
        oldpars=getattr(kernel_module, 'oldpars', {}),
        tests=getattr(kernel_module, 'tests', []),
        )
    process_parameters(model_info)
    # Check for optional functions
    functions = "ER VR form_volume Iq Iqxy shape sesans".split()
    model_info.update((k, getattr(kernel_module, k, None)) for k in functions)
    create_default_functions(model_info)
    # Precalculate the monodisperse parameters
    # TODO: make this a lazy evaluator
    # make_model_info is called for every model on sasview startup
    model_info['mono_details'] = mono_details(model_info)
    return model_info

section_marker = re.compile(r'\A(?P<first>[%s])(?P=first)*\Z'
                            %re.escape(string.punctuation))
def _convert_section_titles_to_boldface(lines):
    """
    Do the actual work of identifying and converting section headings.
    """
    prior = None
    for line in lines:
        if prior is None:
            prior = line
        elif section_marker.match(line):
            if len(line) >= len(prior):
                yield "".join(("**", prior, "**"))
                prior = None
            else:
                yield prior
                prior = line
        else:
            yield prior
            prior = line
    if prior is not None:
        yield prior

def convert_section_titles_to_boldface(s):
    """
    Use explicit bold-face rather than section headings so that the table of
    contents is not polluted with section names from the model documentation.

    Sections are identified as the title line followed by a line of punctuation
    at least as long as the title line.
    """
    return "\n".join(_convert_section_titles_to_boldface(s.split('\n')))

def make_doc(model_info):
    """
    Return the documentation for the model.
    """
    Iq_units = "The returned value is scaled to units of |cm^-1| |sr^-1|, absolute scale."
    Sq_units = "The returned value is a dimensionless structure factor, $S(q)$."
    docs = convert_section_titles_to_boldface(model_info['docs'])
    subst = dict(id=model_info['id'].replace('_', '-'),
                 name=model_info['name'],
                 title=model_info['title'],
                 parameters=make_partable(model_info['parameters']),
                 returns=Sq_units if model_info['structure_factor'] else Iq_units,
                 docs=docs)
    return DOC_HEADER % subst



def demo_time():
    """
    Show how long it takes to process a model.
    """
    from .models import cylinder
    import datetime
    tic = datetime.datetime.now()
    make_source(make_model_info(cylinder))
    toc = (datetime.datetime.now() - tic).total_seconds()
    print("time: %g"%toc)

def main():
    """
    Program which prints the source produced by the model.
    """
    if len(sys.argv) <= 1:
        print("usage: python -m sasmodels.generate modelname")
    else:
        name = sys.argv[1]
        import sasmodels.models
        __import__('sasmodels.models.' + name)
        model = getattr(sasmodels.models, name)
        model_info = make_model_info(model)
        source = make_source(model_info)
        print(source)

if __name__ == "__main__":
    main()
