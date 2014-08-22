__all__ = ["make_opencl"]

import os.path

import numpy as np

from .jsonutil import relaxed_loads

F64 = np.dtype('float64')
F32 = np.dtype('float32')

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

# Scale and background, which are parameters common to every form factor
COMMON_PARAMETERS = [
    [ "scale", "", 0, [0, np.inf], "", "Source intensity" ],
    [ "background", "1/cm", 0, [0, np.inf], "", "Source background" ],
    ]


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
IQ_KERNEL = {
    'fn': "Iq",
    'q_par_decl': "global const real *q,",
    'qinit': "const real qi = q[i];",
    'qcall': "qi",
    }

IQXY_KERNEL = {
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
# to call the volume function from the user supplied kernel, and accumulate
# a normalized weight.
VOLUME_NORM="""const real vol_weight = %(weight)s;
  vol += vol_weight*volume(%(pars)s);
  norm_vol += vol_weight;"""

def indent(s, depth):
    """
    Indent a string of text with *depth* additional spaces on each line.
    """
    spaces = " "*depth
    sep = "\n"+spaces
    return spaces + sep.join(s.split("\n"))


def make_kernel(meta, form):
    """
    Build a kernel call from metadata supplied by the user.

    *meta* is the json object defined in the kernel file.

    *form* is either "Iq" or "Iqxy".

    This does not create a complete OpenCL kernel source, only the top
    level kernel call with polydispersity and a call to the appropriate
    Iq or Iqxy function.
    """

    # If we are building the Iqxy kernel, we need to propagate qx,qy
    # parameters, otherwise we can
    if form == "Iqxy":
        qpars = IQXY_KERNEL
    else:
        qpars = IQ_KERNEL

    depth = 4
    offset = ""
    loop_head = []
    loop_end = []
    vol_pars = []
    fixed_pars = []
    pd_pars = []
    fn_pars = []
    for i,p in enumerate(meta['parameters']):
        name = p[0]
        ptype = p[4]
        if ptype == "volume":
            vol_pars.append(name)
        elif ptype == "orientation":
            if form != "Iqxy": continue  # no orientation for 1D kernels
        elif ptype == "magnetic":
            raise NotImplementedError("no magnetic parameters yet")
        if name not in ['scale','background']: fn_pars.append(name)
        if ptype == "":
            fixed_pars.append(name)
            continue
        else:
            pd_pars.append(name)
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
    subst = {
        'weight_product': "*".join(p+"_w" for p in pd_pars),
        'volume_norm': volume_norm,
        'fn': qpars['fn'],
        'qcall': qpars['qcall'],
        'pcall': ", ".join(fn_pars),
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
        'name': "_".join((meta['name'], qpars['fn'])),
        # to declare, e.g., global real q[],
        'q_par_decl': qpars['q_par_decl'],
        # to declare, e.g., real sld, int Nradius, int Nlength
        'par_decl': par_decl,
        # to copy global to local pd pars we need, e.g., Nradius+Nlength
        'pd_length': "+".join('N'+p for p in pd_pars),
        # the q initializers, e.g., real qi = q[i];
        'qinit': qpars['qinit'],
        # the actual polydispersity loop
        'loops': loops,
        }
    kernel = KERNEL_TEMPLATE%subst
    return kernel

def make_partable(meta):
    pars = meta['parameters']
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

def make_doc(kernelfile, meta, doc):
    doc = doc%{'parameters': make_partable(meta)}
    return doc

def make_model(kernelfile, meta, source):
    kernel_Iq = make_kernel(meta, "Iq")
    kernel_Iqxy = make_kernel(meta, "Iqxy")
    path = os.path.dirname(kernelfile)
    extra = [open("%s/%s"%(path,f)).read() for f in meta['include']]
    kernel = "\n\n".join([KERNEL_HEADER]+extra+[source, kernel_Iq, kernel_Iqxy])
    return kernel

def parse_file(kernelfile):
    source = open(kernelfile).read()

    # select parameters out of the source file
    parts = source.split("PARAMETERS")
    if len(parts) != 3:
        raise ValueError("PARAMETERS block missing from %r"%kernelfile)
    meta_source = parts[1].strip()
    try:
        meta = relaxed_loads(meta_source)
    except:
        print "in json text:"
        print "\n".join("%2d: %s"%(i+1,s)
                        for i,s in enumerate(meta_source.split('\n')))
        raise
        #raise ValueError("PARAMETERS block could not be parsed from %r"%kernelfile)
    meta['parameters'] = COMMON_PARAMETERS + meta['parameters']
    meta['filename'] = kernelfile

    # select documentation out of the source file
    parts = source.split("DOCUMENTATION")
    if len(parts) == 3:
        doc = make_doc(kernelfile, meta, parts[1].strip())
    elif len(parts) == 1:
        raise ValueError("DOCUMENTATION block is missing from %r"%kernelfile)
    else:
        raise ValueError("DOCUMENTATION block incorrect from %r"%kernelfile)

    return source, meta, doc

def make(kernelfile):
    """
    Build an OpenCL function from the source in *kernelfile*.

    The kernel file needs to define metadata about the parameters.  This
    will be a JSON definition containing
    """
    #print kernelfile
    source, meta, doc = parse_file(kernelfile)
    doc = make_doc(kernelfile, meta, doc)
    model = make_model(kernelfile, meta, source)
    return model, meta, doc



# Convert from python float to C float or double, depending on dtype
FLOAT_CONVERTER = {
    F32: np.float32,
    F64: np.float64,
    }

def kernel_name(meta, is_2D):
    return meta['name'] + "_" + ("Iqxy" if is_2D else "Iq")


def kernel_pars(pars, par_info, is_2D, dtype=F32):
    """
    Convert parameter dictionary into arguments for the kernel.

    *pars* is a dictionary of *{ name: value }*, with *name_pd* for the
    polydispersity width, *name_pd_n* for the number of pd steps, and
    *name_pd_nsigma* for the polydispersity limits.

    *par_info* is the parameter info structure from the kernel metadata.

    *is_2D* is True if the dataset represents 2D data, with the corresponding
    orientation parameters.

    *dtype* is F32 or F64, the numpy single and double precision floating
    point dtypes.  These should not be the strings.
    """
    from .weights import GaussianDispersion
    real = np.float32 if dtype == F32 else np.float64
    fixed = []
    parts = []
    for p in par_info['parameters']:
        name, ptype = p[0],p[4]
        value = pars[name]
        if ptype == "":
            fixed.append(real(value))
        elif is_2D or ptype != "orientation":
            limits, width = p[3], pars[name+'_pd']
            n, nsigma = pars[name+'_pd_n'], pars[name+'_pd_nsigma']
            relative = (ptype != "orientation")
            dist = GaussianDispersion(int(n), width, nsigma)
            # Make sure that weights are normalized to peaks at 1 so that
            # the tolerance term can be used properly on truncated distributions
            v,w = dist.get_weights(value, limits[0], limits[1], relative)
            parts.append((v, w/w.max()))
    loops = np.hstack(parts)
    loops = np.ascontiguousarray(loops.T, dtype).flatten()
    loopsN = [np.uint32(len(p[0])) for p in parts]
    return fixed, loops, loopsN


def demo_time():
    import datetime
    tic = datetime.datetime.now()
    toc = lambda: (datetime.datetime.now()-tic).total_seconds()
    path = os.path.dirname("__file__")
    doc, c = make_model(os.path.join(path, "models", "cylinder.c"))
    print "time:",toc()

def demo():
    from os.path import join as joinpath, dirname
    c, meta, doc = make_model(joinpath(dirname(__file__), "models", "cylinder.c"))
    #print doc
    #print c

if __name__ == "__main__":
    demo()
