"""
SAS model constructor.

Small angle scattering models are defined by a set of kernel functions:

    *Iq(q, p1, p2, ...)* returns the scattering at q for a form with
    particular dimensions averaged over all orientations.

    *Iqxy(qx, qy, p1, p2, ...)* returns the scattering at qx,qy for a form
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
functions written for OpenCL.  Floating point values should be
declared as *real*.  Depending on how the function is called, a macro
will replace *real* with *float* or *double*.  Unfortunately, MacOSX
is picky about floating point constants, which should be defined with
value + 'f' if they are of type *float* or just a bare value if they
are of type *double*.  The solution is a macro *REAL(value)* which
adds the 'f' if compiling for single precision floating point.  This
does make the code ugly, and may someday be replaced by a clever
regular expression which does the same job.  OpenCL has a *sincos*
function which can improve performance when both the *sin* and *cos*
values are needed for a particular argument.  Since this function
does not exist in C-99, all use of *sincos* should be replaced by the
macro *SINCOS(value,sn,cn)* where *sn* and *cn* are previously declared
*real* values.  *value* may be an expression.  When compiled for systems
without OpenCL, *SINCOS* will be replaced by *sin* and *cos* calls.  All
functions need prototype declarations even if the are defined before they
are used -- another present from MacOSX.  OpenCL does not support
*#include* preprocessor directives; instead the includes must be listed
in the kernel metadata, with functions defined before they are used.
The included files should be listed using relative path to the kernel
source file, or if using one of the standard models, relative to the
sasmodels source files.

*ER* and *VR* are python functions which operate on parameter vectors.
The constructor code will generate the necessary vectors for computing
them with the desired polydispersity.

The kernel module must set variables defining the kernel meta data:

    *name* is the model name

    *title* is a short description of the model, suitable for a tool tip,
    or a one line model summary in a table of models.

    *description* is an extended description of the model to be displayed
    while the model parameters are being edited.

    *parameters* is a list of parameters.  Each parameter is itself a
    list containing *name*, *units*, *default value*,
    [*lower bound*, *upper bound*], *type* and *description*.

    *source* is the list of C-99 source files that must be joined to
    create the OpenCL kernel functions.  The files defining the functions
    need to be listed before the files which use the functions.

    *ER* is a python function defining the effective radius.  If it is
    not present, the effective radius is 0.

    *VR* is a python function defining the volume ratio.  If it is not
    present, the volume ratio is 1.

The doc string at the start of the kernel module will be used to
construct the model documentation web pages.  Embedded figures should
appear in the subdirectory "img" beside the model definition, and tagged
with the kernel module name to avoid collision with other models.  Some
file systems are case-sensitive, so only use lower case characters for
file names and extensions.

Parameters are defined as follows:

    *name* is the name that will be used in the call to the kernel
    function and the name that will be displayed to the user.  Names
    should be lower case, with words separated by underscore.  If
    acronyms are used, the whole acronym should be upper case.

    *units* should be one of *degrees* for angles, *Ang* for lengths,
    *1e-6/Ang^2* for SLDs.

    *default value* will be the initial value for  the model when it
    is selected, or when an initial value is not otherwise specified.

    *limits* are the valid bounds of the parameter, used to limit the
    polydispersity density function.   In the fit, the parameter limits
    given to the fit are the limits  on the central value of the parameter.
    If there is polydispersity, it will evaluate parameter values outside
    the fit limits, but not outside the hard limits specified in the model.
    If there are no limits, use +/-inf imported from numpy.

    *type* indicates how the parameter will be used.  "volume" parameters
    will be used in all functions.  "orientation" parameters will be used
    in *Iqxy* and *Imagnetic*.  "magnetic* parameters will be used in
    *Imagnetic* only.  If *type* is none, the parameter will be used in
    all of *Iq*, *Iqxy* and *Imagnetic*.  This will probably be a
    is the empty string if the parameter is used in all model calculations,
    it is "volu

    *description* is a short description of the parameter.  This will
    be displayed in the parameter table and used as a tool tip for the
    parameter value in the user interface.

The function :func:`make` loads the metadata from the module and returns
the kernel source.  The function :func:`doc` extracts
"""

# TODO: identify model files which have changed since loading and reload them.

__all__ = ["make, doc"]

import os.path

import numpy as np

F64 = np.dtype('float64')
F32 = np.dtype('float32')

# Scale and background, which are parameters common to every form factor
COMMON_PARAMETERS = [
    [ "scale", "", 1, [0, np.inf], "", "Source intensity" ],
    [ "background", "1/cm", 0, [0, np.inf], "", "Source background" ],
    ]


# Conversion from units defined in the parameter table for each model
# to units displayed in the sphinx documentation.
RST_UNITS = {
    "Ang": "|Ang|",
    "1/Ang^2": "|Ang^-2|",
    "1e-6/Ang^2": "|1e-6Ang^-2|",
    "degrees": "degree",
    "1/cm": "|cm^-1|",
    "": "None",
    }

# Headers for the parameters tables in th sphinx documentation
PARTABLE_HEADERS = [
    "Parameter name",
    "Units",
    "Default value",
    ]

PARTABLE_VALUE_WIDTH = 10

# Header included before every kernel.
# This makes sure that the appropriate math constants are defined, and
KERNEL_HEADER = """\
// GENERATED CODE --- DO NOT EDIT ---
// Code is produced by sasmodels.gen from sasmodels/models/MODEL.c

#ifdef __OPENCL_VERSION__
# define USE_OPENCL
#endif

// If opencl is not available, then we are compiling a C function
// Note: if using a C++ compiler, then define kernel as extern "C"
#ifndef USE_OPENCL
#  include <math.h>
#  define REAL(x) (x)
#  ifndef real
#      define real double
#  endif
#  define global
#  define local
#  define constant const
#  define kernel
#  define SINCOS(angle,svar,cvar) do {svar=sin(angle);cvar=cos(angle);} while (0)
#  define powr(a,b) pow(a,b)
#else
#  ifdef USE_SINCOS
#    define SINCOS(angle,svar,cvar) svar=sincos(angle,&cvar)
#  else
#    define SINCOS(angle,svar,cvar) do {svar=sin(angle);cvar=cos(angle);} while (0)
#  endif
#endif

// Standard mathematical constants, prefixed with M_:
//   E, LOG2E, LOG10E, LN2, LN10, PI, PI_2, PI_4, 1_PI, 2_PI,
//   2_SQRTPI, SQRT2, SQRT1_2
// OpenCL defines M_constant_F for float constants, and nothing if double
// is not enabled on the card, which is why these constants may be missing
#ifndef M_PI
#  define M_PI REAL(3.141592653589793)
#endif
#ifndef M_PI_2
#  define M_PI_2 REAL(1.570796326794897)
#endif
#ifndef M_PI_4
#  define M_PI_4 REAL(0.7853981633974483)
#endif

// Non-standard pi/180, used for converting between degrees and radians
#ifndef M_PI_180
#  define M_PI_180 REAL(0.017453292519943295)
#endif
"""


# The I(q) kernel and the I(qx, qy) kernel have one and two q parameters
# respectively, so the template builder will need to do extra work to
# declare, initialize and pass the q parameters.
KERNEL_1D = {
    'fn': "Iq",
    'q_par_decl': "global const real *q,",
    'qinit': "const real qi = q[i];",
    'qcall': "qi",
    }

KERNEL_2D = {
    'fn': "Iqxy",
    'q_par_decl': "global const real *qx,\n    global const real *qy,",
    'qinit': "const real qxi = qx[i];\n    const real qyi = qy[i];",
    'qcall': "qxi, qyi",
    }

# Generic kernel template for opencl/openmp.
# This defines the opencl kernel that is available to the host.  The same
# structure is used for Iq and Iqxy kernels, so extra flexibility is needed
# for q parameters.  The polydispersity loop is built elsewhere and
# substituted into this template.
KERNEL_TEMPLATE = """\
kernel void %(name)s(
    %(q_par_decl)s
    global real *result,
#ifdef USE_OPENCL
    global real *loops_g,
#else
    const int Nq,
#endif
    local real *loops,
    const real cutoff,
    %(par_decl)s
    )
{
#ifdef USE_OPENCL
  // copy loops info to local memory
  event_t e = async_work_group_copy(loops, loops_g, (%(pd_length)s)*2, 0);
  wait_group_events(1, &e);

  int i = get_global_id(0);
  int Nq = get_global_size(0);
#endif

#ifdef USE_OPENCL
  if (i < Nq)
#else
  #pragma omp parallel for
  for (int i=0; i < Nq; i++)
#endif
  {
    %(qinit)s
    real ret=REAL(0.0), norm=REAL(0.0);
    real vol=REAL(0.0), norm_vol=REAL(0.0);
%(loops)s
    if (vol*norm_vol != REAL(0.0)) {
      ret *= norm_vol/vol;
    }
    result[i] = scale*ret/norm+background;
  }
}
"""

# Polydispersity loop level.
# This pulls the parameter value and weight from the looping vector in order
# in preperation for a nested loop.
LOOP_OPEN="""\
for (int %(name)s_i=0; %(name)s_i < N%(name)s; %(name)s_i++) {
  const real %(name)s = loops[2*(%(name)s_i%(offset)s)];
  const real %(name)s_w = loops[2*(%(name)s_i%(offset)s)+1];"""

# Polydispersity loop body.
# This computes the weight, and if it is sufficient, calls the scattering
# function and adds it to the total.  If there is a volume normalization,
# it will also be added here.
LOOP_BODY="""\
const real weight = %(weight_product)s;
if (weight > cutoff) {
  ret += weight*%(fn)s(%(qcall)s, %(pcall)s);
  norm += weight;
  %(volume_norm)s
}"""

# Volume normalization.
# If there are "volume" polydispersity parameters, then these will be used
# to call the form_volume function from the user supplied kernel, and accumulate
# a normalized weight.
VOLUME_NORM="""const real vol_weight = %(weight)s;
  vol += vol_weight*form_volume(%(pars)s);
  norm_vol += vol_weight;"""

# Documentation header for the module, giving the model name, its short
# description and its parameter table.  The remainder of the doc comes
# from the module docstring.
DOC_HEADER=""".. _%(name)s:

%(name)s
=======================================================

%(title)s

%(parameters)s

The returned value is scaled to units of |cm^-1|.

%(docs)s
"""

def indent(s, depth):
    """
    Indent a string of text with *depth* additional spaces on each line.
    """
    spaces = " "*depth
    sep = "\n"+spaces
    return spaces + sep.join(s.split("\n"))


def kernel_name(info, is_2D):
    return info['name'] + "_" + ("Iqxy" if is_2D else "Iq")


def make_kernel(info, is_2D):
    """
    Build a kernel call from metadata supplied by the user.

    *info* is the json object defined in the kernel file.

    *form* is either "Iq" or "Iqxy".

    This does not create a complete OpenCL kernel source, only the top
    level kernel call with polydispersity and a call to the appropriate
    Iq or Iqxy function.
    """

    # If we are building the Iqxy kernel, we need to propagate qx,qy
    # parameters, otherwise we can
    dim = "2d" if is_2D else "1d"
    fixed_pars = info['partype']['fixed-'+dim]
    pd_pars = info['partype']['pd-'+dim]
    vol_pars = info['partype']['volume']
    q_pars = KERNEL_2D if is_2D else KERNEL_1D

    # Build polydispersity loops
    depth = 4
    offset = ""
    loop_head = []
    loop_end = []
    for name in pd_pars:
        subst = { 'name': name, 'offset': offset }
        loop_head.append(indent(LOOP_OPEN%subst, depth))
        loop_end.insert(0, (" "*depth) + "}")
        offset += '+N'+name
        depth += 2

    # The volume parameters in the inner loop are used to call the volume()
    # function in the kernel, with the parameters defined in vol_pars and the
    # weight product defined in weight.  If there are no volume parameters,
    # then there will be no volume normalization.
    if vol_pars:
        subst = {
            'weight': "*".join(p+"_w" for p in vol_pars),
            'pars': ", ".join(vol_pars),
            }
        volume_norm = VOLUME_NORM%subst
    else:
        volume_norm = ""

    # Define the inner loop function call
    # The parameters to the f(q,p1,p2...) call should occur in the same
    # order as given in the parameter info structure.  This may be different
    # from the parameter order in the call to the kernel since the kernel
    # call places all fixed parameters before all polydisperse parameters.
    fq_pars = [p[0] for p in info['parameters'][len(COMMON_PARAMETERS):]
               if p[0] in set(fixed_pars+pd_pars)]
    subst = {
        'weight_product': "*".join(p+"_w" for p in pd_pars),
        'volume_norm': volume_norm,
        'fn': q_pars['fn'],
        'qcall': q_pars['qcall'],
        'pcall': ", ".join(fq_pars), # skip scale and background
        }
    loop_body = [indent(LOOP_BODY%subst, depth)]
    loops = "\n".join(loop_head+loop_body+loop_end)

    # declarations for non-pd followed by pd pars
    # e.g.,
    #     const real sld,
    #     const int Nradius
    fixed_par_decl = ",\n    ".join("const real %s"%p for p in fixed_pars)
    pd_par_decl = ",\n    ".join("const int N%s"%p for p in pd_pars)
    if fixed_par_decl and pd_par_decl:
        par_decl = ",\n    ".join((fixed_par_decl, pd_par_decl))
    elif fixed_par_decl:
        par_decl = fixed_par_decl
    else:
        par_decl = pd_par_decl

    # Finally, put the pieces together in the kernel.
    subst = {
        # kernel name is, e.g., cylinder_Iq
        'name': kernel_name(info, is_2D),
        # to declare, e.g., global real q[],
        'q_par_decl': q_pars['q_par_decl'],
        # to declare, e.g., real sld, int Nradius, int Nlength
        'par_decl': par_decl,
        # to copy global to local pd pars we need, e.g., Nradius+Nlength
        'pd_length': "+".join('N'+p for p in pd_pars),
        # the q initializers, e.g., real qi = q[i];
        'qinit': q_pars['qinit'],
        # the actual polydispersity loop
        'loops': loops,
        }
    kernel = KERNEL_TEMPLATE%subst
    return kernel

def make_partable(info):
    pars = info['parameters']
    column_widths = [
        max(len(p[0]) for p in pars),
        max(len(RST_UNITS[p[1]]) for p in pars),
        PARTABLE_VALUE_WIDTH,
        ]
    column_widths = [max(w, len(h))
                     for w,h in zip(column_widths, PARTABLE_HEADERS)]

    sep = " ".join("="*w for w in column_widths)
    lines = [
        sep,
        " ".join("%-*s"%(w,h) for w,h in zip(column_widths, PARTABLE_HEADERS)),
        sep,
        ]
    for p in pars:
        lines.append(" ".join([
            "%-*s"%(column_widths[0],p[0]),
            "%-*s"%(column_widths[1],RST_UNITS[p[1]]),
            "%*g"%(column_widths[2],p[2]),
            ]))
    lines.append(sep)
    return "\n".join(lines)

def _search(search_path, filename):
    """
    Find *filename* in *search_path*.

    Raises ValueError if file does not exist.
    """
    for path in search_path:
        target = os.path.join(path, filename)
        if os.path.exists(target):
            return target
    raise ValueError("%r not found in %s"%(filename, search_path))

def make_model(search_path, info):
    kernel_Iq = make_kernel(info, is_2D=False)
    kernel_Iqxy = make_kernel(info, is_2D=True)
    source = [open(_search(search_path, f)).read() for f in info['source']]
    kernel = "\n\n".join([KERNEL_HEADER]+source+[kernel_Iq, kernel_Iqxy])
    return kernel

def categorize_parameters(pars):
    """
    Build parameter categories out of the the parameter definitions.

    Returns a dictionary of categories.

    The function call sequence consists of q inputs and the return vector,
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
    """
    partype = {
        'volume': [], 'orientation': [], 'magnetic': [], '': [],
        'fixed-1d': [], 'fixed-2d': [], 'pd-1d': [], 'pd-2d': [],
        'pd-rel': set(),
    }

    for p in pars:
        name,ptype = p[0],p[4]
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
            raise ValueError("unknown parameter type %r"%ptype)
        partype[ptype].append(name)

    return partype

def make(kernel_module):
    """
    Build an OpenCL/ctypes function from the definition in *kernel_module*.

    The module can be loaded with a normal python import statement if you
    know which module you need, or with __import__('sasmodels.model.'+name)
    if the name is in a string.
    """
    # TODO: allow Iq and Iqxy to be defined in python
    from os.path import abspath, dirname, join as joinpath
    #print kernelfile
    info = dict(
        filename = abspath(kernel_module.__file__),
        name = kernel_module.name,
        title = kernel_module.title,
        source = kernel_module.source,
        description = kernel_module.description,
        parameters = COMMON_PARAMETERS + kernel_module.parameters,
        ER = getattr(kernel_module, 'ER', None),
        VR = getattr(kernel_module, 'VR', None),
        )
    info['limits'] = dict((p[0],p[3]) for p in info['parameters'])
    info['partype'] = categorize_parameters(info['parameters'])

    search_path = [ dirname(info['filename']),
                    abspath(joinpath(dirname(__file__),'models')) ]
    source = make_model(search_path, info)

    return source, info

def doc(kernel_module):
    """
    Return the documentation for the model.
    """
    subst = dict(name=kernel_module.name,
                 title=kernel_module.title,
                 parameters=make_partable(kernel_module.parameters),
                 doc=kernel_module.__doc__)
    return DOC_HEADER%subst


def demo_time():
    import datetime
    tic = datetime.datetime.now()
    toc = lambda: (datetime.datetime.now()-tic).total_seconds()
    path = os.path.dirname("__file__")
    doc, c = make_model(os.path.join(path, "models", "cylinder.c"))
    print "time:",toc()

def demo():
    from os.path import join as joinpath, dirname
    c, info, doc = make_model(joinpath(dirname(__file__), "models", "cylinder.c"))
    #print doc
    #print c

if __name__ == "__main__":
    demo()
