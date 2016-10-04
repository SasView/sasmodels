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
    dimension, or 1.0 if no volume normalization is required.

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

A :class:`modelinfo.ModelInfo` structure is constructed from the kernel meta
data and returned to the caller.

The doc string at the start of the kernel module will be used to
construct the model documentation web pages.  Embedded figures should
appear in the subdirectory "img" beside the model definition, and tagged
with the kernel module name to avoid collision with other models.  Some
file systems are case-sensitive, so only use lower case characters for
file names and extensions.

Code follows the C99 standard with the following extensions and conditions::

    M_PI_180 = pi/180
    M_4PI_3 = 4pi/3
    square(x) = x*x
    cube(x) = x*x*x
    sinc(x) = sin(x)/x, with sin(0)/0 -> 1
    all double precision constants must include the decimal point
    all double declarations may be converted to half, float, or long double
    FLOAT_SIZE is the number of bytes in the converted variables

:func:`load_kernel_module` loads the model definition file and
:func:`modelinfo.make_model_info` parses it. :func:`make_source`
converts C-based model definitions to C source code, including the
polydispersity integral.  :func:`model_sources` returns the list of
source files the model depends on, and :func:`timestamp` returns
the latest time stamp amongst the source files (so you can check if
the model needs to be rebuilt).

The function :func:`make_doc` extracts the doc string and adds the
parameter table to the top.  *make_figure* in *sasmodels/doc/genmodel*
creates the default figure for the model.  [These two sets of code
should mignrate into docs.py so docs can be updated in one place].
"""
from __future__ import print_function

# TODO: determine which functions are useful outside of generate
#__all__ = ["model_info", "make_doc", "make_source", "convert_type"]

import sys
from os.path import abspath, dirname, join as joinpath, exists, isdir, getmtime
import re
import string

import numpy as np  # type: ignore

from .modelinfo import Parameter
from .custom import load_custom_kernel_module

try:
    from typing import Tuple, Sequence, Iterator, Dict
    from .modelinfo import ModelInfo
except ImportError:
    pass

def get_data_path(external_dir, target_file):
    path = abspath(dirname(__file__))
    if exists(joinpath(path, target_file)):
        return path

    # check next to exe/zip file
    exepath = dirname(sys.executable)
    path = joinpath(exepath, external_dir)
    if exists(joinpath(path, target_file)):
        return path

    # check in py2app Contents/Resources
    path = joinpath(exepath, '..', 'Resources', external_dir)
    if exists(joinpath(path, target_file)):
        return abspath(path)

    raise RuntimeError('Could not find '+joinpath(external_dir, target_file))

EXTERNAL_DIR = 'sasmodels-data'
DATA_PATH = get_data_path(EXTERNAL_DIR, 'kernel_template.c')
MODEL_PATH = joinpath(DATA_PATH, 'models')

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
    "g/cm^3": "|g/cm^3|",
    "mg/m^2": "|mg/m^2|",
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
    # type: (str) -> str
    """
    Convert units into ReStructured Text format.
    """
    return "string" if isinstance(units, list) else RST_UNITS.get(units, units)


def make_partable(pars):
    # type: (List[Parameter]) -> str
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
    # type: (List[str], str) -> str
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
    # type: (ModelInfo) -> List[str]
    """
    Return a list of the sources file paths for the module.
    """
    search_path = [dirname(model_info.filename), MODEL_PATH]
    return [_search(search_path, f) for f in model_info.source]


def dll_timestamp(model_info):
    # type: (ModelInfo) -> int
    """
    Return a timestamp for the model corresponding to the most recently
    changed file or dependency.
    """
    # TODO: fails DRY; templates appear two places.
    model_templates = [joinpath(DATA_PATH, filename)
                       for filename in ('kernel_header.c', 'kernel_iq.c')]
    source_files = (model_sources(model_info)
                    + model_templates
                    + [model_info.filename])
    # Note: file may not exist when it is a standard model from library.zip
    times = [getmtime(f) for f in source_files if exists(f)]
    newest = max(times) if times else 0
    return newest

def ocl_timestamp(model_info):
    # type: (ModelInfo) -> int
    """
    Return a timestamp for the model corresponding to the most recently
    changed file or dependency.

    Note that this does not look at the time stamps for the OpenCL header
    information since that need not trigger a recompile of the DLL.
    """
    # TODO: fails DRY; templates appear two places.
    model_templates = [joinpath(DATA_PATH, filename)
                       for filename in ('kernel_header.c', 'kernel_iq.cl')]
    source_files = (model_sources(model_info)
                    + model_templates
                    + [model_info.filename])
    # Note: file may not exist when it is a standard model from library.zip
    times = [getmtime(f) for f in source_files if exists(f)]
    newest = max(times) if times else 0
    return newest


def convert_type(source, dtype):
    # type: (str, np.dtype) -> str
    """
    Convert code from double precision to the desired type.

    Floating point constants are tagged with 'f' for single precision or 'L'
    for long double precision.
    """
    source = _fix_tgmath_int(source)
    if dtype == F16:
        fbytes = 2
        source = _convert_type(source, "half", "f")
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
        raise ValueError("Unexpected dtype in source conversion: %s" % dtype)
    return ("#define FLOAT_SIZE %d\n" % fbytes)+source


def _convert_type(source, type_name, constant_flag):
    # type: (str, str, str) -> str
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
    source = _tag_float(source, constant_flag)
    return source

TGMATH_INT_RE = re.compile(r"""
(?: # Non-capturing match; not lookbehind since pattern length is variable
  \b              # word boundary
   # various math functions
  (a?(sin|cos|tan)h? | atan2
   | erfc? | tgamma
   | exp(2|10|m1)? | log(2|10|1p)? | pow[nr]? | sqrt | rsqrt | rootn
   | fabs | fmax | fmin
   )
  \s*[(]\s*       # open parenthesis
)
[+-]?(0|[1-9]\d*) # integer
(?=               # lookahead match: don't want to move from end of int
  \s*[,)]         # comma or close parenthesis for end of argument
)                 # end lookahead
""", re.VERBOSE)
def _fix_tgmath_int(source):
    # type: (str) -> str
    """
    Replace f(integer) with f(integer.) for sin, cos, pow, etc.

    OS X OpenCL complains that it can't resolve the type generic calls to
    the standard math functions when they are called with integer constants,
    but this does not happen with the Windows Intel driver for example.
    To avoid confusion on the matrix marketplace, automatically promote
    integers to floats if we recognize them in the source.

    The specific functions we look for are:

        trigonometric: sin, asin, sinh, asinh, etc., and atan2
        exponential:   exp, exp2, exp10, expm1, log, log2, log10, logp1
        power:         pow, pown, powr, sqrt, rsqrt, rootn
        special:       erf, erfc, tgamma
        float:         fabs, fmin, fmax

    Note that we don't convert the second argument of dual argument
    functions: atan2, fmax, fmin, pow, powr.  This could potentially
    be a problem for pow(x, 2), but that case seems to work without change.
    """
    out = TGMATH_INT_RE.sub(r'\g<0>.', source)
    return out


# Floating point regular expression
#
# Define parts:
#
#    E = [eE][+-]?\d+    : Exponent
#    P = [.]             : Decimal separator
#    N = [1-9]\d*        : Natural number, no leading zeros
#    Z = 0               : Zero
#    F = \d+             : Fractional number, maybe leading zeros
#    F? = \d*            : Optional fractional number
#
# We want to reject bare natural numbers and bare decimal points, so we
# need to tediously outline the cases where we have either a fraction or
# an exponent:
#
#   ( ZP | ZPF | ZE | ZPE | ZPFE | NP | NPF | NE | NPE | NPFE | PF | PFE )
#
#
# We can then join cases by making parts optional.  The following are
# some ways to do this:
#
#   ( (Z|N)(P|PF|E|PE|PFE) | PFE? )                   # Split on lead
#     => ( (Z|N)(PF?|(PF?)?E) | PFE? )
#   ( ((Z|N)PF?|PF)E? | (Z|N)E)                       # Split on point
#   ( (ZP|ZPF|NP|NPF|PF) | (Z|ZP|ZPF|N|NP|NPF|PF)E )  # Split on E
#     => ( ((Z|N)PF?|PF) | ((Z|N)(PF?)? | PF) E )
FLOAT_RE = re.compile(r"""
    (?<!\w)  # use negative lookbehind since '.' confuses \b test
    # use split on lead to match float ( (Z|N)(PF?|(PF?)?E) | PFE? )
    ( ( 0 | [1-9]\d* )                     # ( ( Z | N )
      ([.]\d* | ([.]\d*)? [eE][+-]?\d+ )   #   (PF? | (PF?)? E )
    | [.]\d+ ([eE][+-]?\d+)?               # | PF (E)?
    )                                      # )
    (?!\w)  # use negative lookahead since '.' confuses \b test
    """, re.VERBOSE)
def _tag_float(source, constant_flag):
    # Convert floating point constants to single by adding 'f' to the end,
    # or long double with an 'L' suffix.  OS/X complains if you don't do this.
    out = FLOAT_RE.sub(r'\g<0>%s'%constant_flag, source)
    #print("in",repr(source),"out",repr(out), constant_flag)
    return out

def test_tag_float():

    cases="""
ZP  : 0.
ZPF : 0.0,0.01,0.1
Z  E: 0e+001
ZP E: 0.E0
ZPFE: 0.13e-031
NP  : 1., 12.
NPF : 1.0001, 1.1, 1.0
N  E: 1e0, 37E-080
NP E: 1.e0, 37.E-080
NPFE: 845.017e+22
 PF : .1, .0, .0100
 PFE: .6e+9, .82E-004
# isolated cases
0.
1e0
0.13e-013
# untouched
struct3.e3, 03.05.67, 37
# expressions
3.75+-1.6e-7-27+13.2
a3.e2 - 0.
4*atan(1)
4.*atan(1.)
"""

    output="""
ZP  : 0.f
ZPF : 0.0f,0.01f,0.1f
Z  E: 0e+001f
ZP E: 0.E0f
ZPFE: 0.13e-031f
NP  : 1.f, 12.f
NPF : 1.0001f, 1.1f, 1.0f
N  E: 1e0f, 37E-080f
NP E: 1.e0f, 37.E-080f
NPFE: 845.017e+22f
 PF : .1f, .0f, .0100f
 PFE: .6e+9f, .82E-004f
# isolated cases
0.f
1e0f
0.13e-013f
# untouched
struct3.e3, 03.05.67, 37
# expressions
3.75f+-1.6e-7f-27+13.2f
a3.e2 - 0.f
4*atan(1)
4.f*atan(1.f)
"""

    for case_in, case_out in zip(cases.split('\n'), output.split('\n')):
        out = _tag_float(case_in, 'f')
        assert case_out == out, "%r => %r"%(case_in, out)


def kernel_name(model_info, variant):
    # type: (ModelInfo, str) -> str
    """
    Name of the exported kernel symbol.

    *variant* is "Iq", "Iqxy" or "Imagnetic".
    """
    return model_info.name + "_" + variant


def indent(s, depth):
    # type: (str, int) -> str
    """
    Indent a string of text with *depth* additional spaces on each line.
    """
    spaces = " "*depth
    sep = "\n" + spaces
    return spaces + sep.join(s.split("\n"))


_template_cache = {}  # type: Dict[str, Tuple[int, str, str]]
def load_template(filename):
    # type: (str) -> str
    path = joinpath(DATA_PATH, filename)
    mtime = getmtime(path)
    if filename not in _template_cache or mtime > _template_cache[filename][0]:
        with open(path) as fid:
            _template_cache[filename] = (mtime, fid.read(), path)
    return _template_cache[filename][1], path


_FN_TEMPLATE = """\
double %(name)s(%(pars)s);
double %(name)s(%(pars)s) {
#line %(line)d "%(filename)s"
    %(body)s
}

"""
def _gen_fn(name, pars, body, filename, line):
    # type: (str, List[Parameter], str, str, int) -> str
    """
    Generate a function given pars and body.

    Returns the following string::

         double fn(double a, double b, ...);
         double fn(double a, double b, ...) {
             ....
         }
    """
    par_decl = ', '.join(p.as_function_argument() for p in pars) if pars else 'void'
    return _FN_TEMPLATE % {
        'name': name, 'pars': par_decl, 'body': body,
        'filename': filename.replace('\\', '\\\\'), 'line': line,
    }


def _call_pars(prefix, pars):
    # type: (str, List[Parameter]) -> List[str]
    """
    Return a list of *prefix.parameter* from parameter items.
    """
    return [p.as_call_reference(prefix) for p in pars]


# type in IQXY pattern could be single, float, double, long double, ...
_IQXY_PATTERN = re.compile("^((inline|static) )? *([a-z ]+ )? *Iqxy *([(]|$)",
                           flags=re.MULTILINE)
def _have_Iqxy(sources):
    # type: (List[str]) -> bool
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
    for path, code in sources:
        if _IQXY_PATTERN.search(code):
            return True
    else:
        return False


def _add_source(source, code, path):
    """
    Add a file to the list of source code chunks, tagged with path and line.
    """
    path = path.replace('\\', '\\\\')
    source.append('#line 1 "%s"' % path)
    source.append(code)

def make_source(model_info):
    # type: (ModelInfo) -> Dict[str, str]
    """
    Generate the OpenCL/ctypes kernel from the module info.

    Uses source files found in the given search path.  Returns None if this
    is a pure python model, with no C source components.
    """
    if callable(model_info.Iq):
        raise ValueError("can't compile python model")
        #return None

    # TODO: need something other than volume to indicate dispersion parameters
    # No volume normalization despite having a volume parameter.
    # Thickness is labelled a volume in order to trigger polydispersity.
    # May want a separate dispersion flag, or perhaps a separate category for
    # disperse, but not volume.  Volume parameters also use relative values
    # for the distribution rather than the absolute values used by angular
    # dispersion.  Need to be careful that necessary parameters are available
    # for computing volume even if we allow non-disperse volume parameters.

    partable = model_info.parameters

    # Load templates and user code
    kernel_header = load_template('kernel_header.c')
    dll_code = load_template('kernel_iq.c')
    ocl_code = load_template('kernel_iq.cl')
    #ocl_code = load_template('kernel_iq_local.cl')
    user_code = [(f, open(f).read()) for f in model_sources(model_info)]

    # Build initial sources
    source = []
    _add_source(source, *kernel_header)
    for path, code in user_code:
        _add_source(source, code, path)

    # Make parameters for q, qx, qy so that we can use them in declarations
    q, qx, qy = [Parameter(name=v) for v in ('q', 'qx', 'qy')]
    # Generate form_volume function, etc. from body only
    if isinstance(model_info.form_volume, str):
        pars = partable.form_volume_parameters
        source.append(_gen_fn('form_volume', pars, model_info.form_volume,
                              model_info.filename, model_info._form_volume_line))
    if isinstance(model_info.Iq, str):
        pars = [q] + partable.iq_parameters
        source.append(_gen_fn('Iq', pars, model_info.Iq,
                              model_info.filename, model_info._Iq_line))
    if isinstance(model_info.Iqxy, str):
        pars = [qx, qy] + partable.iqxy_parameters
        source.append(_gen_fn('Iqxy', pars, model_info.Iqxy,
                              model_info.filename, model_info._Iqxy_line))

    # Define the parameter table
    # TODO: plug in current line number
    source.append('#line 542 "sasmodels/generate.py"')
    source.append("#define PARAMETER_TABLE \\")
    source.append("\\\n".join(p.as_definition()
                              for p in partable.kernel_parameters))

    # Define the function calls
    if partable.form_volume_parameters:
        refs = _call_pars("_v.", partable.form_volume_parameters)
        call_volume = "#define CALL_VOLUME(_v) form_volume(%s)"%(",".join(refs))
    else:
        # Model doesn't have volume.  We could make the kernel run a little
        # faster by not using/transferring the volume normalizations, but
        # the ifdef's reduce readability more than is worthwhile.
        call_volume = "#define CALL_VOLUME(v) 1.0"
    source.append(call_volume)

    refs = ["_q[_i]"] + _call_pars("_v.", partable.iq_parameters)
    call_iq = "#define CALL_IQ(_q,_i,_v) Iq(%s)" % (",".join(refs))
    if _have_Iqxy(user_code) or isinstance(model_info.Iqxy, str):
        # Call 2D model
        refs = ["_q[2*_i]", "_q[2*_i+1]"] + _call_pars("_v.", partable.iqxy_parameters)
        call_iqxy = "#define CALL_IQ(_q,_i,_v) Iqxy(%s)" % (",".join(refs))
    else:
        # Call 1D model with sqrt(qx^2 + qy^2)
        #warnings.warn("Creating Iqxy = Iq(sqrt(qx^2 + qy^2))")
        # still defined:: refs = ["q[i]"] + _call_pars("v", iq_parameters)
        pars_sqrt = ["sqrt(_q[2*_i]*_q[2*_i]+_q[2*_i+1]*_q[2*_i+1])"] + refs[1:]
        call_iqxy = "#define CALL_IQ(_q,_i,_v) Iq(%s)" % (",".join(pars_sqrt))

    magpars = [k-2 for k,p in enumerate(partable.call_parameters)
               if p.type == 'sld']

    # Fill in definitions for numbers of parameters
    source.append("#define MAX_PD %s"%partable.max_pd)
    source.append("#define NUM_PARS %d"%partable.npars)
    source.append("#define NUM_VALUES %d" % partable.nvalues)
    source.append("#define NUM_MAGNETIC %d" % partable.nmagnetic)
    source.append("#define MAGNETIC_PARS %s"%",".join(str(k) for k in magpars))
    for k,v in enumerate(magpars[:3]):
        source.append("#define MAGNETIC_PAR%d %d"%(k+1, v))

    # TODO: allow mixed python/opencl kernels?

    ocl = kernels(ocl_code, call_iq, call_iqxy, model_info.name)
    dll = kernels(dll_code, call_iq, call_iqxy, model_info.name)
    result = {
        'dll': '\n'.join(source+dll[0]+dll[1]+dll[2]),
        'opencl': '\n'.join(source+ocl[0]+ocl[1]+ocl[2]),
    }

    return result


def kernels(kernel, call_iq, call_iqxy, name):
    # type: ([str,str], str, str, str) -> List[str]
    code = kernel[0]
    path = kernel[1].replace('\\', '\\\\')
    iq = [
        # define the Iq kernel
        "#define KERNEL_NAME %s_Iq" % name,
        call_iq,
        '#line 1 "%s Iq"' % path,
        code,
        "#undef CALL_IQ",
        "#undef KERNEL_NAME",
        ]

    iqxy = [
        # define the Iqxy kernel from the same source with different #defines
        "#define KERNEL_NAME %s_Iqxy" % name,
        call_iqxy,
        '#line 1 "%s Iqxy"' % path,
        code,
        "#undef CALL_IQ",
        "#undef KERNEL_NAME",
         ]

    imagnetic = [
        # define the Imagnetic kernel
        "#define KERNEL_NAME %s_Imagnetic" % name,
        "#define MAGNETIC 1",
        call_iqxy,
        '#line 1 "%s Imagnetic"' % path,
        code,
        "#undef MAGNETIC",
        "#undef CALL_IQ",
        "#undef KERNEL_NAME",
    ]

    return iq, iqxy, imagnetic


def load_kernel_module(model_name):
    # type: (str) -> module
    """
    Return the kernel module named in *model_name*.

    If the name ends in *.py* then load it as a custom model using
    :func:`sasmodels.custom.load_custom_kernel_module`, otherwise
    load it from :mod:`sasmodels.models`.
    """
    if model_name.endswith('.py'):
        kernel_module = load_custom_kernel_module(model_name)
    else:
        from sasmodels import models
        __import__('sasmodels.models.'+model_name)
        kernel_module = getattr(models, model_name, None)
    return kernel_module


section_marker = re.compile(r'\A(?P<first>[%s])(?P=first)*\Z'
                            % re.escape(string.punctuation))
def _convert_section_titles_to_boldface(lines):
    # type: (Sequence[str]) -> Iterator[str]
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
    # type: (str) -> str
    """
    Use explicit bold-face rather than section headings so that the table of
    contents is not polluted with section names from the model documentation.

    Sections are identified as the title line followed by a line of punctuation
    at least as long as the title line.
    """
    return "\n".join(_convert_section_titles_to_boldface(s.split('\n')))


def make_doc(model_info):
    # type: (ModelInfo) -> str
    """
    Return the documentation for the model.
    """
    Iq_units = "The returned value is scaled to units of |cm^-1| |sr^-1|, absolute scale."
    Sq_units = "The returned value is a dimensionless structure factor, $S(q)$."
    docs = model_info.docs if model_info.docs is not None else ""
    docs = convert_section_titles_to_boldface(docs)
    pars = make_partable(model_info.parameters.COMMON
                         + model_info.parameters.kernel_parameters)
    subst = dict(id=model_info.id.replace('_', '-'),
                 name=model_info.name,
                 title=model_info.title,
                 parameters=pars,
                 returns=Sq_units if model_info.structure_factor else Iq_units,
                 docs=docs)
    return DOC_HEADER % subst


# TODO: need a single source for rst_prolog; it is also in doc/rst_prolog
RST_PROLOG = """\
.. |Ang| unicode:: U+212B
.. |Ang^-1| replace:: |Ang|\ :sup:`-1`
.. |Ang^2| replace:: |Ang|\ :sup:`2`
.. |Ang^-2| replace:: |Ang|\ :sup:`-2`
.. |1e-6Ang^-2| replace:: 10\ :sup:`-6`\ |Ang|\ :sup:`-2`
.. |Ang^3| replace:: |Ang|\ :sup:`3`
.. |Ang^-3| replace:: |Ang|\ :sup:`-3`
.. |Ang^-4| replace:: |Ang|\ :sup:`-4`
.. |cm^-1| replace:: cm\ :sup:`-1`
.. |cm^2| replace:: cm\ :sup:`2`
.. |cm^-2| replace:: cm\ :sup:`-2`
.. |cm^3| replace:: cm\ :sup:`3`
.. |1e15cm^3| replace:: 10\ :sup:`15`\ cm\ :sup:`3`
.. |cm^-3| replace:: cm\ :sup:`-3`
.. |sr^-1| replace:: sr\ :sup:`-1`
.. |P0| replace:: P\ :sub:`0`\

.. |equiv| unicode:: U+2261
.. |noteql| unicode:: U+2260
.. |TM| unicode:: U+2122

.. |cdot| unicode:: U+00B7
.. |deg| unicode:: U+00B0
.. |g/cm^3| replace:: g\ |cdot|\ cm\ :sup:`-3`
.. |mg/m^2| replace:: mg\ |cdot|\ m\ :sup:`-2`
.. |fm^2| replace:: fm\ :sup:`2`
.. |Ang*cm^-1| replace:: |Ang|\ |cdot|\ cm\ :sup:`-1`
"""

# TODO: make a better fake reference role
RST_ROLES = """\
.. role:: ref

.. role:: numref

"""

def make_html(model_info):
    """
    Convert model docs directly to html.
    """
    from . import rst2html

    rst = make_doc(model_info)
    return rst2html.rst2html("".join((RST_ROLES, RST_PROLOG, rst)))

def view_html(model_name):
    from . import rst2html
    from . import modelinfo
    kernel_module = load_kernel_module(model_name)
    info = modelinfo.make_model_info(kernel_module)
    url = "file://"+dirname(info.filename)+"/"
    rst2html.wxview(make_html(info), url=url)

def demo_time():
    # type: () -> None
    """
    Show how long it takes to process a model.
    """
    import datetime
    from .modelinfo import make_model_info
    from .models import cylinder

    tic = datetime.datetime.now()
    make_source(make_model_info(cylinder))
    toc = (datetime.datetime.now() - tic).total_seconds()
    print("time: %g"%toc)


def main():
    # type: () -> None
    """
    Program which prints the source produced by the model.
    """
    import sys
    from .modelinfo import make_model_info

    if len(sys.argv) <= 1:
        print("usage: python -m sasmodels.generate modelname")
    else:
        name = sys.argv[1]
        kernel_module = load_kernel_module(name)
        model_info = make_model_info(kernel_module)
        source = make_source(model_info)
        print(source['dll'])


if __name__ == "__main__":
    main()
