"""
SAS model constructor.

Small angle scattering models are defined by a set of kernel functions:

    *Iq(q, p1, p2, ...)* returns the scattering at q for a form with
    particular dimensions averaged over all orientations.

    *Iqac(qab, qc, p1, p2, ...)* returns the scattering at qab, qc
    for a rotationally symmetric form with particular dimensions.
    qab, qc are determined from shape orientation and scattering angles.
    This call is used if the shape has orientation parameters theta and phi.

    *Iqabc(qa, qb, qc, p1, p2, ...)* returns the scattering at qa, qb, qc
    for a form with particular dimensions.  qa, qb, qc are determined from
    shape orientation and scattering angles. This call is used if the shape
    has orientation parameters theta, phi and psi.

    *Iqxy(qx, qy, p1, p2, ...)* returns the scattering at qx, qy.  Use this
    to create an arbitrary 2D theory function, needed for q-dependent
    background functions and for models with non-uniform magnetism.

    *form_volume(p1, p2, ...)* returns the volume of the form with particular
    dimension, or 1.0 if no volume normalization is required.

    *shell_volume(p1, p2, ...)* returns the volume of the shell for forms
    which are hollow.

    *radius_effective(mode, p1, p2, ...)* returns the effective radius of
    the form with particular dimensions.  Mode determines the type of
    effective radius returned, with mode=1 for equivalent volume.

These functions are defined in a kernel module .py script and an associated
set of .c files.  The model constructor will use them to create models with
polydispersity across volume and orientation parameters, and provide
scale and background parameters for each model.

C code should be stylized C-99 functions written for OpenCL. All functions
need prototype declarations even if the are defined before they are used.
Although OpenCL supports *#include* preprocessor directives, the list of
includes should be given as part of the metadata in the kernel module
definition. The included files should be listed using a path relative to the
kernel module, or if using "lib/file.c" if it is one of the standard includes
provided with the sasmodels source. The includes need to be listed in order
so that functions are defined before they are used.

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

    *form_volume*, *Iq*, *Iqac*, *Iqabc* are strings containing
    the C source code for the body of the volume, Iq, and Iqac functions
    respectively.  These can also be defined in the last source file.

    *Iq*, *Iqac*, *Iqabc* also be instead be python functions defining the
    kernel.  If they are marked as *Iq.vectorized = True* then the
    kernel is passed the entire *q* vector at once, otherwise it is
    passed values one *q* at a time.  The performance improvement of
    this step is significant.

    *valid* expression that evaluates to True if the input parameters
    are valid (e.g., "bell_radius >= radius" for the barbell or capped
    cylinder models).  The expression can call C functions, including
    those defined in your model file.

A :class:`.modelinfo.ModelInfo` structure is constructed from the kernel meta
data and returned to the caller.

Valid inputs should be identified by the *valid* expression.
Particularly with polydispersity, there are some sets of shape parameters
which lead to nonsensical forms, such as a capped cylinder where the
cap radius is smaller than the cylinder radius.  The polydispersity
calculation will ignore these points, effectively chopping the
parameter weight distributions at the boundary of the infeasible region.
The resulting scattering will be set to background, even for models with
no polydispersity.  If the *valid* expression misses some parameter
combinations and they reach the kernel, the kernel should probably return
NaN rather than zero.  Even if the volume also evaluates to zero for
these parameters, the distribution weights are still accumulated and the
average volume calculation will be slightly off.

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
    sas_sinx_x(x) = sin(x)/x, with sin(0)/0 -> 1
    all double precision constants must include the decimal point
    all double declarations may be converted to half, float, or long double
    FLOAT_SIZE is the number of bytes in the converted variables

:func:`load_kernel_module` loads the model definition file and
:func:`.modelinfo.make_model_info` parses it. :func:`make_source`
converts C-based model definitions to C source code, including the
polydispersity integral.  :func:`model_sources` returns the list of
source files the model depends on, and :func:`ocl_timestamp` returns
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
from os import environ
from os.path import abspath, dirname, join as joinpath, exists, getmtime, sep
import re
import string
from zlib import crc32
from inspect import currentframe, getframeinfo
import logging

import numpy as np  # type: ignore

from .modelinfo import Parameter
from .custom import load_custom_kernel_module

# pylint: disable=unused-import
try:
    from typing import Tuple, Sequence, Iterator, Dict, List
    from types import ModuleType
    from .modelinfo import ModelInfo
except ImportError:
    pass
# pylint: enable=unused-import

logger = logging.getLogger(__name__)

# jitter projection to use in the kernel code.  See explore/jitter.py
# for details.  To change it from a program, set generate.PROJECTION.
PROJECTION = 1

def get_data_path(external_dir, target_file):
    """
    Search for the target file relative in the installed application.

    Search first in the location of the generate module in case we are
    running directly from the distribution.  Search next to the python
    executable for windows installs.  Search in the ../Resources directory
    next to the executable for Mac OS/X installs.
    """
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
DATA_PATH = get_data_path(EXTERNAL_DIR, 'kernel_iq.c')
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
# This section associates the unit with the macro to use to produce the LaTex
# code.  The macro itself needs to be defined in sasmodels/doc/rst_prolog.
#
# NOTE: there is an RST_PROLOG at the end of this file which is NOT
# used for the bundled documentation. Still as long as we are defining the
# macros in two places any new addition should define the macro in both places.
RST_UNITS = {
    "Ang": "|Ang|",
    "1/Ang": "|Ang^-1|",
    "1/Ang^2": "|Ang^-2|",
    "Ang^3": "|Ang^3|",
    "Ang^2": "|Ang^2|",
    "1e15/cm^3": "|1e15cm^3|",
    "Ang^3/mol": "|Ang^3|/mol",
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


def set_integration_size(info, n):
    # type: (ModelInfo, int) -> None
    """
    Update the model definition, replacing the gaussian integration with
    a gaussian integration of a different size.

    Note: this really ought to be a method in modelinfo, but that leads to
    import loops.
    """
    if info.source and any(lib.startswith('lib/gauss') for lib in info.source):
        from .gengauss import gengauss
        path = joinpath(MODEL_PATH, "lib", "gauss%d.c"%n)
        if not exists(path):
            gengauss(n, path)
        info.source = ["lib/gauss%d.c"%n if lib.startswith('lib/gauss')
                       else lib for lib in info.source]

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

    underbar = " ".join("="*w for w in column_widths)
    lines = [
        underbar,
        " ".join("%-*s" % (w, h)
                 for w, h in zip(column_widths, PARTABLE_HEADERS)),
        underbar,
        ]
    for p in pars:
        lines.append(" ".join([
            "%-*s" % (column_widths[0], p.name),
            "%-*s" % (column_widths[1], p.description),
            "%-*s" % (column_widths[2], format_units(p.units)),
            "%*g" % (column_widths[3], p.default),
            ]))
    lines.append(underbar)
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
    search_path = [dirname(model_info.basefile), MODEL_PATH]
    # If it is a reparameterized model, include both to directory for the
    # original model file and the directory for the reparameterized model file
    # on the search path for dependent source files.
    if model_info.filename and model_info.filename is not model_info.basefile:
        search_path.insert(0, dirname(model_info.filename))
    return [_search(search_path, f) for f in model_info.source]


def dll_timestamp(model_info):
    # type: (ModelInfo) -> int
    """
    Return a timestamp for the model corresponding to the most recently
    changed file or dependency.
    """
    # TODO: what about on-the-fly derived models---tag name with crc32
    # TODO: no longer used, and appears to be identical to ocl_timestamp
    # TODO: fails DRY; templates appear two places.
    model_templates = [joinpath(DATA_PATH, filename)
                       for filename in ('kernel_header.c', 'kernel_iq.c')]
    source_files = (model_sources(model_info)
                    + model_templates
                    + [model_info.basefile])
    # Include filename if different from basefile.
    if model_info.filename and model_info.filename is not model_info.basefile:
        source_files.append(model_info.filename)
    # Note: file may not exist when it is a standard model from library.zip
    # or when the model is generated dynamically.
    times = [getmtime(f) for f in source_files if f and exists(f)]
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
                       for filename in ('kernel_header.c', 'kernel_iq.c')]
    source_files = (model_sources(model_info)
                    + model_templates
                    + [model_info.basefile])
    # Include filename if different from basefile.
    if model_info.filename and model_info.filename is not model_info.basefile:
        source_files.append(model_info.filename)
    # Note: file may not exist when it is a standard model from library.zip
    # or when the model is generated dynamically.
    times = [getmtime(f) for f in source_files if exists(f)]
    newest = max(times) if times else 0
    return newest

def tag_source(source):
    # type: (str) -> str
    """
    Return a unique tag for the source code.
    """
    # Note: need 0xffffffff&val to force an unsigned 32-bit number
    try:
        source = source.encode('utf8')
    except AttributeError: # bytes has no encode attribute in python 3
        pass
    return "%08X"%(0xffffffff&crc32(source))

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
    # Assume complex numbers are represented as cdouble.
    # TODO: use lookahead so that "typedef complex double cdouble;" converts.
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
    """Check that floating point constants are identified and tagged with 'f'"""

    cases = """
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

    output = """
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
    interline_separator = "\n" + spaces
    return spaces + interline_separator.join(s.split("\n"))


_template_cache = {}  # type: Dict[str, Tuple[int, str, str]]
def load_template(filename):
    # type: (str) -> str
    """
    Load template file from sasmodels resource directory.
    """
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
def _gen_fn(model_info, name, pars):
    # type: (ModelInfo, str, List[Parameter]) -> str
    """
    Generate a function given pars and body.

    Returns the following string::

         double fn(double a, double b, ...);
         double fn(double a, double b, ...) {
             ....
         }
    """
    par_decl = ', '.join(p.as_function_argument() for p in pars) if pars else 'void'
    body = getattr(model_info, name)
    filename = model_info.basefile
    # Note: if symbol is defined strangely in the module then default it to 1
    lineno = model_info.lineno.get(name, 1)
    return _FN_TEMPLATE % {
        'name': name, 'pars': par_decl, 'body': body,
        'filename': filename.replace('\\', '/'), 'line': lineno,
    }


def _call_pars(pars, subs):
    # type: (str, List[Parameter]) -> List[str]
    """
    Return a list of *prefix+parameter* from parameter items.

    *pars* is the list of parameters from the base model.

    *subs* contains the translation equations with references to parameters
    from the new parameter table.  If there is no translation, then *subs*
    is just a list of references into the base table.
    """
    return [subs[p.id] for p in pars]

def _split_translation(translation):
    r"""
    Process the *translation* string, which is a sequence of assignments.

    Blanks and comments (c-style and python-style) are stripped.

    Conditional expressions should use C syntax (! || && ? :) not python.
    """
    # TODO: support translation for vector parameters
    # TODO: security concerns when code stored in model file
    assignments = []
    for line in translation.split('\n'):
        # strip comments and white space from line
        code = line.split('#', 1)[0].split('//', 1)[0].strip()
        if not code:
            continue
        # split assignment on first equals
        parts = code.split('=', 1)
        if len(parts) == 1:
            raise ValueError("translation expected 'var=expr' but got %r"%line)
        # store parts without assessing whether they are correct
        var, expr = parts
        assignments.append((var.strip(), expr.strip()))
    return assignments

# The phrase "(?<![.0-9])" in the following re indicates negative look-behind,
# designed to reject the like 314159e-5 and 1.0f
_IDENT_RE = re.compile(r"(?<![.0-9])([A-Za-z_][A-Za-z0-9_]*)")

def _build_translation(model_info, table_id="_v", var_prefix="_var_"):
    """
    Interpret parameter translation block, if any.

    *model_info* contains the parameter table and any translation definition
    for converting between model parameters and the base calculation model.

    *table_id* is the internal label used for the parameter call table.  It
    must be of the form "_table" matching whatever table variable name that
    appears in the macros such as CALL_VOLUME() and CALL_IQ().

    *var_prefix* is a tag to attach to intermediate variables to avoid
    collision with variables used inside kernel_iq.

    Returns:

    * *subs* = {name: expr, ...} parameter substitution table for calling
      into the kernel function
    * *translation* = "#define TRANSLATION_VARS(_v) _var_name = expr ...."
    * *validity* = "#define VALID(_v) ..."

    The returned *subs* is used to generate the substitions for CALL_VOLUME
    etc., parameters, via :func:`_call_pars`.

    The returned *translation* and *validity* macros need to be included inside
    the generated model files.  They are the same for all variants (1D, 2D,
    magnetic) so can be defined once along side the parameter table.  Even
    though they are expanded independently in each variant.
    """
    # Usual case of no translation.
    call_table = model_info.parameters.kernel_parameters
    if not model_info.translation:
        # Simple substitution {name: table_id.name}, an empty
        # TRANSLATION_VARS() macro and VALID() macro with simple subs.
        subs = {p.id: table_id+'.'+p.id for p in call_table}
        translation_vars = _build_translation_vars(table_id, [])
        validity = _build_validity_check(model_info.valid, table_id, subs)
        return subs, translation_vars, validity

    # Translation case.
    # Grab the list of assignment statements.
    assigns = _split_translation(model_info.translation)

    # Identify parameters and intermediates.  Parameters are from the calling
    # table, not the base table.  Only these have values defined by the caller
    # that can be used in expressions.  Need the list of base parameters to
    # identify substitutions.  And LHS variables that are not base parameters
    # are therefore intermediate parameters.  Note: leaving the main parameter
    # table as the caller parameter table because the rest of the code uses it
    # to build up and interpret the calling sequence.  The kernel functions
    # themselves of course use the base table.
    base_table = model_info.base.kernel_parameters
    call_pars = set(p.id for p in call_table)
    base_pars = set(p.id for p in base_table)
    variables = set(name for name, expr in assigns
                    if name not in call_pars and name not in base_pars)

    # Regular expression substition to tag the variables on each RHS.
    # Don't want to translate all symbols since they include C constants and
    # math functions.
    def id_sub(match_obj):
        name = match_obj.group(0)
        if name in variables:
            return var_prefix + name
        if name in call_pars:
            return table_id + '.' + name
        return name
    assigns = [(name, _IDENT_RE.sub(id_sub, eq)) for name, eq in assigns]

    # Create a substituion table assuming all parameters are directly assigned,
    # then update it with
    subs = {p.id: table_id+'.'+p.id for p in base_table}
    subs.update((name, eq) for name, eq in assigns if name in base_pars)

    # Build TRANSLATION_VARS macros assignments to variables.
    var_values = [(var_prefix+name, eq) for name, eq in assigns
                  if name in variables]
    translation_vars = _build_translation_vars(table_id, var_values)

    # Build VALID macro, with expressions using caller parameters.
    validity = _build_validity_check(model_info.valid, table_id, subs)

    return subs, translation_vars, validity

def _build_translation_vars(table_id, variables):
    r"""
    Build TRANSLATION_VARS macro for C which builds intermediate values.

    E. g.,

    ::

        #define TRANSLATION_VARS(_v) \
         const double _temporary_Re = cbrt(_v.volume/_v.eccentricity/M_4PI_3)
    """
    if variables:
        # Leave comma off def since last def doesn't have comma
        defs = ["const double %s = %s" % (name, eq)
                for name, eq in variables]
        return ("#define TRANSLATION_VARS(%s) \\\n  "%table_id
                + "; \\\n  ".join(defs))
    else:
        return '#define TRANSLATION_VARS(%s) do {} while(0)'%table_id

def _build_validity_check(eq, table_id, subs):
    """
    Substitute parameter expressions into validity test, returning the
    VALID(_table) macro.
    """
    if not eq:  # None or empty string
        return "#define VALID(%s) 1" % (table_id, )
    else:
        # TODO: maybe replace not/and/or with ! && ||
        # TODO: (hard) replace (a if b else c) with (b ? a : c)
        # TODO: reject a**n, or (hard) replace with pow(a, n)
        def id_sub(match_obj):
            name = match_obj.group(0)
            return "(%s)"%subs[name] if name in subs else name
        eq = _IDENT_RE.sub(id_sub, eq)
        return "#define VALID(%s) (%s)" % (table_id, eq)

_IQXY_PATTERN = re.compile(r"(^|\s)double\s+I(?P<mode>q(ac|abc|xy))\s*[(]",
                           flags=re.MULTILINE)
def find_xy_mode(source):
    # type: (List[str]) -> bool
    """
    Return the xy mode as qa, qac, qabc or qxy.

    Note this is not a C parser, and so can be easily confused by
    non-standard syntax.  Also, it will incorrectly identify the following
    as having 2D models::

        /*
        double Iqac(qab, qc, ...) { ... fill this in later ... }
        */

    If you want to comment out the function, use // on the front of the
    line::

        /*
        // double Iqac(qab, qc, ...) { ... fill this in later ... }
        */

    """
    for code in source:
        m = _IQXY_PATTERN.search(code)
        if m is not None:
            return m.group('mode')
    return 'qa'

# Note: not presently used.  Need to know whether Fq is available before
# trying to compile the source.  Leave the code here in case we decide that
# define have_Fq for each form factor is too tedious and error prone.
_FQ_PATTERN = re.compile(r"(^|\s)void\s+Fq[(]", flags=re.MULTILINE)
def contains_Fq(source):
    # type: (List[str]) -> bool
    """
    Return True if C source defines "void Fq(".
    """
    for code in source:
        if _FQ_PATTERN.search(code) is not None:
            return True
    return False

_SHELL_VOLUME_PATTERN = re.compile(r"(^|\s)double\s+shell_volume[(]",
                                   flags=re.MULTILINE)
def contains_shell_volume(source):
    # type: (List[str]) -> bool
    """
    Return True if C source defines "double shell_volume(".
    """
    for code in source:
        if _SHELL_VOLUME_PATTERN.search(code) is not None:
            return True
    return False

def _add_source(source, code, path, lineno=1):
    """
    Add a file to the list of source code chunks, tagged with path and line.
    """
    path = path.replace('\\', '/')
    source.append('#line %d "%s"' % (lineno, path))
    source.append(code)

def read_text(f):
    with open(f) as fid:
        return fid.read()

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

    # reference to table of calling parameters and table of called
    # parameters, which may be different if the model is reparameterized
    # by the caller.
    call_table = model_info.parameters
    base_table = model_info.base

    # Load templates and user code
    kernel_header = load_template('kernel_header.c')
    kernel_code = load_template('kernel_iq.c')
    user_code = [(f, read_text(f)) for f in model_sources(model_info)]

    # Build initial sources
    source = []
    _add_source(source, *kernel_header)
    for path, code in user_code:
        _add_source(source, code, path)
    if model_info.c_code:
        _add_source(source, model_info.c_code, model_info.basefile,
                    lineno=model_info.lineno.get('c_code', 1))

    # Make parameters for q, qx, qy so that we can use them in declarations
    q, qx, qy, qab, qa, qb, qc \
        = [Parameter(name=v) for v in 'q qx qy qab qa qb qc'.split()]

    # Generate form_volume, etc., function definitions from strings
    # given in the model.py file as the C function body only, such as
    # one might for the volume of a sphere:
    #     form_volume="return M_4PI_3*cube(radius)"
    if isinstance(model_info.form_volume, str):
        pars = call_table.form_volume_parameters
        source.append(_gen_fn(model_info, 'form_volume', pars))
    if isinstance(model_info.shell_volume, str):
        pars = call_table.form_volume_parameters
        source.append(_gen_fn(model_info, 'shell_volume', pars))
    if isinstance(model_info.Iq, str):
        pars = [q] + call_table.iq_parameters
        source.append(_gen_fn(model_info, 'Iq', pars))
    if isinstance(model_info.Iqxy, str):
        pars = [qx, qy] + call_table.iq_parameters + call_table.orientation_parameters
        source.append(_gen_fn(model_info, 'Iqxy', pars))
    if isinstance(model_info.Iqac, str):
        pars = [qab, qc] + call_table.iq_parameters
        source.append(_gen_fn(model_info, 'Iqac', pars))
    if isinstance(model_info.Iqabc, str):
        pars = [qa, qb, qc] + call_table.iq_parameters
        source.append(_gen_fn(model_info, 'Iqabc', pars))

    # Check for shell_volume function in source code
    is_hollow = contains_shell_volume(source)

    # What kind of 2D model do we need?  Is it consistent with the parameters?
    xy_mode = find_xy_mode(source)
    if xy_mode == 'qabc' and not base_table.is_asymmetric:
        raise ValueError("asymmetric oriented models need to define Iqabc")
    elif xy_mode == 'qac' and base_table.is_asymmetric:
        raise ValueError("symmetric oriented models need to define Iqac")
    elif not base_table.orientation_parameters and xy_mode in ('qac', 'qabc'):
        raise ValueError("Unexpected function I%s for unoriented shape"%xy_mode)
    elif base_table.orientation_parameters and xy_mode not in ('qac', 'qabc'):
        if xy_mode == 'qxy':
            logger.warning("oriented shapes should define Iqac or Iqabc")
        else:
            raise ValueError("Expected Iqac or Iqabc for oriented shape")

    # Process parameter substitutions
    subs, translation_vars, valid = _build_translation(model_info, '_v')

    # Validity check macro.  Putting it here so that line number
    # jiggering is short-lived.
    _add_source(source, valid, model_info.basefile,
                lineno=model_info.lineno.get('valid', 1))
    source.append(valid)

    # Define the parameter table
    lineno = getframeinfo(currentframe()).lineno + 2
    source.append('#line %d "sasmodels/generate.py"'%lineno)
    #source.append('introduce breakage in generate to test lineno reporting')
    source.append("#define PARAMETER_TABLE \\")
    source.append("\\\n".join(p.as_definition()
                              for p in call_table.kernel_parameters))

    # Define the function calls
    call_radius_effective = "#define CALL_RADIUS_EFFECTIVE(_mode, _v) 0.0"
    if base_table.form_volume_parameters:
        refs = _call_pars(base_table.form_volume_parameters, subs)
        if is_hollow:
            call_volume = (
                "#define CALL_VOLUME(_form, _shell, _v) "
                "do { _form = form_volume(%s); _shell = shell_volume(%s); } "
                "while (0)") % ((",".join(refs),)*2)
        else:
            call_volume = (
                "#define CALL_VOLUME(_form, _shell, _v) "
                "do { _form = _shell = form_volume(%s); } "
                "while (0)") % (",".join(refs))
        if model_info.radius_effective_modes:
            call_radius_effective = (
                "#define CALL_RADIUS_EFFECTIVE(_mode, _v) "
                "radius_effective(_mode, %s)") % (",".join(refs))
    else:
        # Model doesn't have volume.  We could make the kernel run a little
        # faster by not using/transferring the volume normalizations, but
        # the ifdef's reduce readability more than is worthwhile.
        call_volume = (
            "#define CALL_VOLUME(_form, _shell, _v) "
            "do { _form = _shell = 1.0; } while (0)")
    source.append(translation_vars)
    source.append(call_volume)
    source.append(call_radius_effective)
    model_refs = _call_pars(base_table.iq_parameters, subs)

    if model_info.have_Fq:
        pars = ",".join(["_q", "&_F1", "&_F2",] + model_refs)
        call_iq = "#define CALL_FQ(_q, _F1, _F2, _v) Fq(%s)" % pars
        clear_iq = "#undef CALL_FQ"
    else:
        pars = ",".join(["_q"] + model_refs)
        call_iq = "#define CALL_IQ(_q, _v) Iq(%s)" % pars
        clear_iq = "#undef CALL_IQ"
    if xy_mode == 'qabc':
        pars = ",".join(["_qa", "_qb", "_qc"] + model_refs)
        call_iqxy = "#define CALL_IQ_ABC(_qa,_qb,_qc,_v) Iqabc(%s)" % pars
        clear_iqxy = "#undef CALL_IQ_ABC"
    elif xy_mode == 'qac':
        pars = ",".join(["_qa", "_qc"] + model_refs)
        call_iqxy = "#define CALL_IQ_AC(_qa,_qc,_v) Iqac(%s)" % pars
        clear_iqxy = "#undef CALL_IQ_AC"
    elif xy_mode == 'qa' and not model_info.have_Fq:
        pars = ",".join(["_qa"] + model_refs)
        call_iqxy = "#define CALL_IQ_A(_qa,_v) Iq(%s)" % pars
        clear_iqxy = "#undef CALL_IQ_A"
    elif xy_mode == 'qa' and model_info.have_Fq:
        pars = ",".join(["_qa", "&_F1", "&_F2",] + model_refs)
        # Note: uses rare C construction (expr1, expr2) which computes
        # expr1 then expr2 and evaluates to expr2.  This allows us to
        # leave it looking like a function even though it is returning
        # its values by reference.
        call_iqxy = "#define CALL_FQ_A(_qa,_F1,_F2,_v) (Fq(%s),_F2)" % pars
        clear_iqxy = "#undef CALL_FQ_A"
    elif xy_mode == 'qxy':
        qxy_refs = _call_pars(base_table.orientation_parameters, subs)
        pars = ",".join(["_qx", "_qy"] + model_refs + qxy_refs)
        call_iqxy = "#define CALL_IQ_XY(_qx,_qy,_v) Iqxy(%s)" % pars
        clear_iqxy = "#undef CALL_IQ_XY"
        if base_table.orientation_parameters:
            call_iqxy += "\n#define HAVE_THETA"
            clear_iqxy += "\n#undef HAVE_THETA"
        if base_table.is_asymmetric:
            call_iqxy += "\n#define HAVE_PSI"
            clear_iqxy += "\n#undef HAVE_PSI"


    # TODO: test magnetism with translated sld parameters
    # By inspection it should work but intermediate variables aren't supported.
    magpars = [k-2 for k, p in enumerate(call_table.call_parameters)
               if p.type == 'sld']
    # Fill in definitions for numbers of parameters
    source.append("#define MAX_PD %s"%call_table.max_pd)
    source.append("#define NUM_PARS %d"%call_table.npars)
    source.append("#define NUM_VALUES %d" % call_table.nvalues)
    source.append("#define NUM_MAGNETIC %d" % call_table.nmagnetic)
    source.append("#define MAGNETIC_PARS %s"%",".join(str(k) for k in magpars))
    source.append("#define PROJECTION %d"%PROJECTION)
    wrappers = _kernels(kernel_code, call_iq, clear_iq,
                        call_iqxy, clear_iqxy, model_info.name)
    code = '\n'.join(source + wrappers[0] + wrappers[1] + wrappers[2])

    # Note: Identical code for dll and opencl.  This may not always be the case
    # so leave them as separate entries in the returned value.
    result = {'dll': code, 'opencl': code}
    return result


def _kernels(kernel, call_iq, clear_iq, call_iqxy, clear_iqxy, name):
    # type: (Dict[str, str], str, str, str, str, str) -> List[str]
    code = kernel[0]
    path = kernel[1].replace('\\', '/')
    iq = [
        # define the Iq kernel
        "#define KERNEL_NAME %s_Iq" % name,
        call_iq,
        '#line 1 "%s Iq"' % path,
        code,
        clear_iq,
        "#undef KERNEL_NAME",
        ]

    iqxy = [
        # define the Iqxy kernel from the same source with different #defines
        "#define KERNEL_NAME %s_Iqxy" % name,
        call_iqxy,
        '#line 1 "%s Iqxy"' % path,
        code,
        clear_iqxy,
        "#undef KERNEL_NAME",
    ]

    imagnetic = [
        # define the Imagnetic kernel
        "#define KERNEL_NAME %s_Imagnetic" % name,
        "#define MAGNETIC 1",
        call_iqxy,
        '#line 1 "%s Imagnetic"' % path,
        code,
        clear_iqxy,
        "#undef MAGNETIC",
        "#undef KERNEL_NAME",
    ]

    return iq, iqxy, imagnetic


def load_kernel_module(model_name):
    # type: (str) -> ModuleType
    """
    Return the kernel module named in *model_name*.

    If the name ends in *.py* then load it as a custom model using
    :func:`.custom.__init__.load_custom_kernel_module`, otherwise
    load it as a builtin from *sasmodels.models*.
    """
    if model_name.endswith('.py'):
        kernel_module = load_custom_kernel_module(model_name)
    else:
        try:
            from sasmodels import models
            __import__('sasmodels.models.'+model_name)
            kernel_module = getattr(models, model_name, None)
        except ImportError:
            # If the model isn't a built in model, try the plugin directory
            plugin_path = environ.get('SAS_MODELPATH', None)
            if plugin_path is not None:
                file_name = model_name.split(sep)[-1]
                model_name = plugin_path + sep + file_name + ".py"
                kernel_module = load_custom_kernel_module(model_name)
            else:
                raise
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
    if model_info.structure_factor:
        pars = model_info.parameters.kernel_parameters
    else:
        pars = (model_info.parameters.common_parameters
                + model_info.parameters.kernel_parameters)
    partable = make_partable(pars)
    subst = dict(id=model_info.id.replace('_', '-'),
                 name=model_info.name,
                 title=model_info.title,
                 parameters=partable,
                 returns=Sq_units if model_info.structure_factor else Iq_units,
                 docs=docs)
    return DOC_HEADER % subst


# TODO: need a single source for rst_prolog; it is also in doc/rst_prolog
RST_PROLOG = r"""\
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
    # type: (ModelInfo) -> str
    """
    Convert model docs directly to html.
    """
    from . import rst2html

    rst = make_doc(model_info)
    return rst2html.rst2html("".join((RST_ROLES, RST_PROLOG, rst)))

def view_html(model_name):
    # type: (str) -> None
    """
    Load the model definition and view its help.
    """
    from . import modelinfo
    kernel_module = load_kernel_module(model_name)
    info = modelinfo.make_model_info(kernel_module)
    view_html_from_info(info)

def view_html_from_info(info):
    # type: (ModelInfo) -> None
    """
    View the help for a loaded model definition.
    """
    from . import rst2html
    url = "file://"+dirname(info.filename)+"/"
    rst2html.view_html(make_html(info), url=url)

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
