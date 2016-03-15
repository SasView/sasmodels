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

    [*lb*, *ub*] are the hard limits on the parameter value, used to limit
    the polydispersity density function.  In the fit, the parameter limits
    given to the fit are the limits  on the central value of the parameter.
    If there is polydispersity, it will evaluate parameter values outside
    the fit limits, but not outside the hard limits specified in the model.
    If there are no limits, use +/-inf imported from numpy.

    *type* indicates how the parameter will be used.  "volume" parameters
    will be used in all functions.  "orientation" parameters will be used
    in *Iqxy* and *Imagnetic*.  "magnetic* parameters will be used in
    *Imagnetic* only.  If *type* is the empty string, the parameter will
    be used in all of *Iq*, *Iqxy* and *Imagnetic*.

    *description* is a short description of the parameter.  This will
    be displayed in the parameter table and used as a tool tip for the
    parameter value in the user interface.

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

# TODO: identify model files which have changed since loading and reload them.

import sys
from os.path import abspath, dirname, join as joinpath, exists, basename, \
    splitext
import re
import string

import numpy as np

#TODO: determine which functions are useful outside of generate
#__all__ = ["model_info", "make_doc", "make_source", "convert_type"]

C_KERNEL_TEMPLATE_PATH = joinpath(dirname(__file__), 'kernel_template.c')

F16 = np.dtype('float16')
F32 = np.dtype('float32')
F64 = np.dtype('float64')
try:  # CRUFT: older numpy does not support float128
    F128 = np.dtype('float128')
except TypeError:
    F128 = None

# Scale and background, which are parameters common to every form factor
COMMON_PARAMETERS = [
    ["scale", "", 1, [0, np.inf], "", "Source intensity"],
    ["background", "1/cm", 0, [0, np.inf], "", "Source background"],
    ]

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
        max(len(p[0]) for p in pars),
        max(len(p[-1]) for p in pars),
        max(len(format_units(p[1])) for p in pars),
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
            "%-*s" % (column_widths[0], p[0]),
            "%-*s" % (column_widths[1], p[-1]),
            "%-*s" % (column_widths[2], format_units(p[1])),
            "%*g" % (column_widths[3], p[2]),
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

# Pragmas for enable OpenCL features.  Be sure to protect them so that they
# still compile even if OpenCL is not present.
_F16_PRAGMA = """\
#if defined(__OPENCL_VERSION__) && !defined(cl_khr_fp16)
#  pragma OPENCL EXTENSION cl_khr_fp16: enable
#endif
"""

_F64_PRAGMA = """\
#if defined(__OPENCL_VERSION__) && !defined(cl_khr_fp64)
#  pragma OPENCL EXTENSION cl_khr_fp64: enable
#endif
"""

def convert_type(source, dtype):
    """
    Convert code from double precision to the desired type.

    Floating point constants are tagged with 'f' for single precision or 'L'
    for long double precision.
    """
    if dtype == F16:
        fbytes = 2
        source = _F16_PRAGMA + _convert_type(source, "half", "f")
    elif dtype == F32:
        fbytes = 4
        source = _convert_type(source, "float", "f")
    elif dtype == F64:
        fbytes = 8
        source = _F64_PRAGMA + source  # Source is already double
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


LOOP_OPEN = """\
for (int %(name)s_i=0; %(name)s_i < N%(name)s; %(name)s_i++) {
  const double %(name)s = loops[2*(%(name)s_i%(offset)s)];
  const double %(name)s_w = loops[2*(%(name)s_i%(offset)s)+1];\
"""
def build_polydispersity_loops(pd_pars):
    """
    Build polydispersity loops

    Returns loop opening and loop closing
    """
    depth = 4
    offset = ""
    loop_head = []
    loop_end = []
    for name in pd_pars:
        subst = {'name': name, 'offset': offset}
        loop_head.append(indent(LOOP_OPEN % subst, depth))
        loop_end.insert(0, (" "*depth) + "}")
        offset += '+N' + name
        depth += 2
    return "\n".join(loop_head), "\n".join(loop_end)

C_KERNEL_TEMPLATE = None
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

    # Load template
    global C_KERNEL_TEMPLATE
    if C_KERNEL_TEMPLATE is None:
        with open(C_KERNEL_TEMPLATE_PATH) as fid:
            C_KERNEL_TEMPLATE = fid.read()

    # Load additional sources
    source = [open(f).read() for f in model_sources(model_info)]

    # Prepare defines
    defines = []
    partype = model_info['partype']
    pd_1d = partype['pd-1d']
    pd_2d = partype['pd-2d']
    fixed_1d = partype['fixed-1d']
    fixed_2d = partype['fixed-1d']

    iq_parameters = [p[0]
                     for p in model_info['parameters'][2:]  # skip scale, background
                     if p[0] in set(fixed_1d + pd_1d)]
    iqxy_parameters = [p[0]
                       for p in model_info['parameters'][2:]  # skip scale, background
                       if p[0] in set(fixed_2d + pd_2d)]
    volume_parameters = [p[0]
                         for p in model_info['parameters']
                         if p[4] == 'volume']

    # Fill in defintions for volume parameters
    if volume_parameters:
        defines.append(('VOLUME_PARAMETERS',
                        ','.join(volume_parameters)))
        defines.append(('VOLUME_WEIGHT_PRODUCT',
                        '*'.join(p + '_w' for p in volume_parameters)))

    # Generate form_volume function from body only
    if model_info['form_volume'] is not None:
        if volume_parameters:
            vol_par_decl = ', '.join('double ' + p for p in volume_parameters)
        else:
            vol_par_decl = 'void'
        defines.append(('VOLUME_PARAMETER_DECLARATIONS',
                        vol_par_decl))
        fn = """\
double form_volume(VOLUME_PARAMETER_DECLARATIONS);
double form_volume(VOLUME_PARAMETER_DECLARATIONS) {
    %(body)s
}
""" % {'body':model_info['form_volume']}
        source.append(fn)

    # Fill in definitions for Iq parameters
    defines.append(('IQ_KERNEL_NAME', model_info['name'] + '_Iq'))
    defines.append(('IQ_PARAMETERS', ', '.join(iq_parameters)))
    if fixed_1d:
        defines.append(('IQ_FIXED_PARAMETER_DECLARATIONS',
                        ', \\\n    '.join('const double %s' % p for p in fixed_1d)))
    if pd_1d:
        defines.append(('IQ_WEIGHT_PRODUCT',
                        '*'.join(p + '_w' for p in pd_1d)))
        defines.append(('IQ_DISPERSION_LENGTH_DECLARATIONS',
                        ', \\\n    '.join('const int N%s' % p for p in pd_1d)))
        defines.append(('IQ_DISPERSION_LENGTH_SUM',
                        '+'.join('N' + p for p in pd_1d)))
        open_loops, close_loops = build_polydispersity_loops(pd_1d)
        defines.append(('IQ_OPEN_LOOPS',
                        open_loops.replace('\n', ' \\\n')))
        defines.append(('IQ_CLOSE_LOOPS',
                        close_loops.replace('\n', ' \\\n')))
    if model_info['Iq'] is not None:
        defines.append(('IQ_PARAMETER_DECLARATIONS',
                        ', '.join('double ' + p for p in iq_parameters)))
        fn = """\
double Iq(double q, IQ_PARAMETER_DECLARATIONS);
double Iq(double q, IQ_PARAMETER_DECLARATIONS) {
    %(body)s
}
""" % {'body':model_info['Iq']}
        source.append(fn)

    # Fill in definitions for Iqxy parameters
    defines.append(('IQXY_KERNEL_NAME', model_info['name'] + '_Iqxy'))
    defines.append(('IQXY_PARAMETERS', ', '.join(iqxy_parameters)))
    if fixed_2d:
        defines.append(('IQXY_FIXED_PARAMETER_DECLARATIONS',
                        ', \\\n    '.join('const double %s' % p for p in fixed_2d)))
    if pd_2d:
        defines.append(('IQXY_WEIGHT_PRODUCT',
                        '*'.join(p + '_w' for p in pd_2d)))
        defines.append(('IQXY_DISPERSION_LENGTH_DECLARATIONS',
                        ', \\\n    '.join('const int N%s' % p for p in pd_2d)))
        defines.append(('IQXY_DISPERSION_LENGTH_SUM',
                        '+'.join('N' + p for p in pd_2d)))
        open_loops, close_loops = build_polydispersity_loops(pd_2d)
        defines.append(('IQXY_OPEN_LOOPS',
                        open_loops.replace('\n', ' \\\n')))
        defines.append(('IQXY_CLOSE_LOOPS',
                        close_loops.replace('\n', ' \\\n')))
    if model_info['Iqxy'] is not None:
        defines.append(('IQXY_PARAMETER_DECLARATIONS',
                        ', '.join('double ' + p for p in iqxy_parameters)))
        fn = """\
double Iqxy(double qx, double qy, IQXY_PARAMETER_DECLARATIONS);
double Iqxy(double qx, double qy, IQXY_PARAMETER_DECLARATIONS) {
    %(body)s
}
""" % {'body':model_info['Iqxy']}
        source.append(fn)

    # Need to know if we have a theta parameter for Iqxy; it is not there
    # for the magnetic sphere model, for example, which has a magnetic
    # orientation but no shape orientation.
    if 'theta' in pd_2d:
        defines.append(('IQXY_HAS_THETA', '1'))

    #for d in defines: print(d)
    defines = '\n'.join('#define %s %s' % (k, v) for k, v in defines)
    sources = '\n\n'.join(source)
    return C_KERNEL_TEMPLATE % {
        'DEFINES': defines,
        'SOURCES': sources,
        }

def categorize_parameters(pars):
    """
    Build parameter categories out of the the parameter definitions.

    Returns a dictionary of categories.

    Note: these categories are subject to change, depending on the needs of
    the UI and the needs of the kernel calling function.

    The categories are as follows:

    * *volume* list of volume parameter names
    * *orientation* list of orientation parameters
    * *magnetic* list of magnetic parameters
    * *<empty string>* list of parameters that have no type info

    Each parameter is in one and only one category.

    The following derived categories are created:

    * *fixed-1d* list of non-polydisperse parameters for 1D models
    * *pd-1d* list of polydisperse parameters for 1D models
    * *fixed-2d* list of non-polydisperse parameters for 2D models
    * *pd-d2* list of polydisperse parameters for 2D models
    """
    partype = {
        'volume': [], 'orientation': [], 'magnetic': [], '': [],
        'fixed-1d': [], 'fixed-2d': [], 'pd-1d': [], 'pd-2d': [],
        'pd-rel': set(),
    }

    for p in pars:
        name, ptype = p[0], p[4]
        if ptype == 'volume':
            partype['pd-1d'].append(name)
            partype['pd-2d'].append(name)
            partype['pd-rel'].add(name)
        elif ptype == 'magnetic':
            partype['fixed-2d'].append(name)
        elif ptype == 'orientation':
            partype['pd-2d'].append(name)
        elif ptype == '':
            partype['fixed-1d'].append(name)
            partype['fixed-2d'].append(name)
        else:
            raise ValueError("unknown parameter type %r" % ptype)
        partype[ptype].append(name)

    return partype

def process_parameters(model_info):
    """
    Process parameter block, precalculating parameter details.
    """
    # Fill in the derived attributes
    model_info['limits'] = dict((p[0], p[3]) for p in model_info['parameters'])
    model_info['partype'] = categorize_parameters(model_info['parameters'])
    model_info['defaults'] = dict((p[0], p[2]) for p in model_info['parameters'])
    if model_info.get('demo', None) is None:
        model_info['demo'] = model_info['defaults']

def make_model_info(kernel_module):
    """
    Interpret the model definition file, categorizing the parameters.

    The module can be loaded with a normal python import statement if you
    know which module you need, or with __import__('sasmodels.model.'+name)
    if the name is in a string.

    The *model_info* structure contains the following fields:

    * *id* is the id of the kernel
    * *name* is the display name of the kernel
    * *title* is a short description of the kernel
    * *description* is a long description of the kernel (this doesn't seem
      very useful since the Help button on the model page brings you directly
      to the documentation page)
    * *docs* is the docstring from the module.  Use :func:`make_doc` to
    * *category* specifies the model location in the docs
    * *parameters* is the model parameter table
    * *single* is True if the model allows single precision
    * *defaults* is the *{parameter: value}* table built from the parameter
      description table.
    * *limits* is the *{parameter: [min, max]}* table built from the
      parameter description table.
    * *partypes* categorizes the model parameters. See
      :func:`categorize_parameters` for details.
    * *demo* contains the *{parameter: value}* map used in compare (and maybe
      for the demo plot, if plots aren't set up to use the default values).
      If *demo* is not given in the file, then the default values will be used.
    * *tests* is a set of tests that must pass
    * *source* is the list of library files to include in the C model build
    * *Iq*, *Iqxy*, *form_volume*, *ER*, and *VR* are python functions
      implementing the kernel for the module, or None if they are not
      defined in python
    * *oldname* is the model name in pre-4.0 Sasview
    * *oldpars* is the *{new: old}* parameter translation table
      from pre-4.0 Sasview
    * *composition* is None if the model is independent, otherwise it is a
      tuple with composition type ('product' or 'mixture') and a list of
      *model_info* blocks for the composition objects.  This allows us to
      build complete product and mixture models from just the info.
    """
    # TODO: maybe turn model_info into a class ModelDefinition
    #print(kernelfile)
    category = getattr(kernel_module, 'category', None)
    parameters = COMMON_PARAMETERS + kernel_module.parameters
    # Default the demo parameters to the starting values for the individual
    # parameters if an explicit demo parameter set has not been specified.
    demo_parameters = getattr(kernel_module, 'demo', None)
    filename = abspath(kernel_module.__file__)
    kernel_id = splitext(basename(filename))[0]
    name = getattr(kernel_module, 'name', None)
    single = getattr(kernel_module, 'single', True)
    if name is None:
        name = " ".join(w.capitalize() for w in kernel_id.split('_'))
    model_info = dict(
        id=kernel_id,  # string used to load the kernel
        filename=abspath(kernel_module.__file__),
        name=name,
        title=kernel_module.title,
        description=kernel_module.description,
        docs=kernel_module.__doc__,
        category=category,
        parameters=parameters,
        composition=None,
        single=single,
        demo=demo_parameters,
        source=getattr(kernel_module, 'source', []),
        oldname=getattr(kernel_module, 'oldname', None),
        oldpars=getattr(kernel_module, 'oldpars', {}),
        tests=getattr(kernel_module, 'tests', []),
        )
    process_parameters(model_info)
    # Fill in attributes which default to None
    model_info.update((k, getattr(kernel_module, k, None))
                      for k in ('ER', 'VR', 'form_volume', 'Iq', 'Iqxy'))
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
    is_Sq = ("structure-factor" in model_info['category'])
    #docs = kernel_module.__doc__
    docs = convert_section_titles_to_boldface(model_info['docs'])
    subst = dict(id=model_info['id'].replace('_', '-'),
                 name=model_info['name'],
                 title=model_info['title'],
                 parameters=make_partable(model_info['parameters']),
                 returns=Sq_units if is_Sq else Iq_units,
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
