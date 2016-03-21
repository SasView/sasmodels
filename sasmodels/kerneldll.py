r"""
DLL driver for C kernels

The global attribute *ALLOW_SINGLE_PRECISION_DLLS* should be set to *True* if
you wish to allow single precision floating point evaluation for the compiled
models, otherwise it defaults to *False*.

The compiler command line is stored in the attribute *COMPILE*, with string
substitutions for %(source)s and %(output)s indicating what to compile and
where to store it.  The actual command is system dependent.

On windows systems, you have a choice of compilers.  *MinGW* is the GNU
compiler toolchain, available in packages such as anaconda and PythonXY,
or available stand alone. This toolchain has had difficulties on some
systems, and may or may not work for you.  In order to build DLLs, *gcc*
must be on your path.  If the environment variable *SAS_OPENMP* is given
then -fopenmp is added to the compiler flags.  This requires a version
of MinGW compiled with OpenMP support.

An alternative toolchain uses the Microsoft Visual C++ compiler, available
free from microsoft:

    `<http://www.microsoft.com/en-us/download/details.aspx?id=44266>`_

Again, this requires that the compiler is available on your path.  This is
done by running vcvarsall.bat in a windows terminal.  Install locations are
system dependent, such as:

    C:\Program Files (x86)\Common Files\Microsoft\Visual C++ for Python\9.0\vcvarsall.bat

or maybe

    C:\Users\yourname\AppData\Local\Programs\Common\Microsoft\Visual C++ for Python\9.0\vcvarsall.bat

And again, the environment variable *SAS_OPENMP* controls whether OpenMP is
used to compile the C code.  This requires the Microsoft vcomp90.dll library,
which doesn't seem to be included with the compiler, nor does there appear
to be a public download location.  There may be one on your machine already
in a location such as:

    C:\Windows\winsxs\x86_microsoft.vc90.openmp*\vcomp90.dll

If you copy this onto your path, such as the python directory or the install
directory for this application, then OpenMP should be supported.
"""
from __future__ import print_function

import sys
import os
import tempfile
import ctypes as ct
from ctypes import c_void_p, c_int32, c_longdouble, c_double, c_float

import numpy as np

from . import generate
from .kernelpy import PyInput, PyModel
from .exception import annotate_exception

# Compiler platform details
if sys.platform == 'darwin':
    #COMPILE = "gcc-mp-4.7 -shared -fPIC -std=c99 -fopenmp -O2 -Wall %s -o %s -lm -lgomp"
    COMPILE = "gcc -shared -fPIC -std=c99 -O2 -Wall %(source)s -o %(output)s -lm"
elif os.name == 'nt':
    # call vcvarsall.bat before compiling to set path, headers, libs, etc.
    if "VCINSTALLDIR" in os.environ:
        # MSVC compiler is available, so use it.  OpenMP requires a copy of
        # vcomp90.dll on the path.  One may be found here:
        #       C:/Windows/winsxs/x86_microsoft.vc90.openmp*/vcomp90.dll
        # Copy this to the python directory and uncomment the OpenMP COMPILE
        # TODO: remove intermediate OBJ file created in the directory
        # TODO: maybe don't use randomized name for the c file
        CC = "cl /nologo /Ox /MD /W3 /GS- /DNDEBUG /Tp%(source)s "
        LN = "/link /DLL /INCREMENTAL:NO /MANIFEST /OUT:%(output)s"
        if "SAS_OPENMP" in os.environ:
            COMPILE = " ".join((CC, "/openmp", LN))
        else:
            COMPILE = " ".join((CC, LN))
    else:
        COMPILE = "gcc -shared -fPIC -std=c99 -O2 -Wall %(source)s -o %(output)s -lm"
        if "SAS_OPENMP" in os.environ:
            COMPILE += " -fopenmp"
else:
    COMPILE = "cc -shared -fPIC -fopenmp -std=c99 -O2 -Wall %(source)s -o %(output)s -lm"

DLL_PATH = tempfile.gettempdir()

ALLOW_SINGLE_PRECISION_DLLS = True


def dll_path(model_info, dtype="double"):
    """
    Path to the compiled model defined by *model_info*.
    """
    from os.path import join as joinpath, split as splitpath, splitext
    basename = splitext(splitpath(model_info['filename'])[1])[0]
    if np.dtype(dtype) == generate.F32:
        basename += "32"
    elif np.dtype(dtype) == generate.F64:
        basename += "64"
    else:
        basename += "128"
    return joinpath(DLL_PATH, basename+'.so')


def make_dll(source, model_info, dtype="double"):
    """
    Load the compiled model defined by *kernel_module*.

    Recompile if any files are newer than the model file.

    *dtype* is a numpy floating point precision specifier indicating whether
    the model should be single or double precision.  The default is double
    precision.

    The DLL is not loaded until the kernel is called so models can
    be defined without using too many resources.

    Set *sasmodels.kerneldll.DLL_PATH* to the compiled dll output path.
    The default is the system temporary directory.

    Set *sasmodels.ALLOW_SINGLE_PRECISION_DLLS* to True if single precision
    models are allowed as DLLs.
    """
    if callable(model_info.get('Iq', None)):
        return PyModel(model_info)

    dtype = np.dtype(dtype)
    if dtype == generate.F16:
        raise ValueError("16 bit floats not supported")
    if dtype == generate.F32 and not ALLOW_SINGLE_PRECISION_DLLS:
        dtype = generate.F64  # Force 64-bit dll

    if dtype == generate.F32: # 32-bit dll
        tempfile_prefix = 'sas_' + model_info['name'] + '32_'
    elif dtype == generate.F64:
        tempfile_prefix = 'sas_' + model_info['name'] + '64_'
    else:
        tempfile_prefix = 'sas_' + model_info['name'] + '128_'

    source = generate.convert_type(source, dtype)
    source_files = (generate.model_sources(model_info)
                    + [model_info['filename']]
                    + generate.model_templates())
    dll = dll_path(model_info, dtype)
    newest = max(os.path.getmtime(f) for f in source_files)
    if not os.path.exists(dll) or os.path.getmtime(dll) < newest:
        # Replace with a proper temp file
        fid, filename = tempfile.mkstemp(suffix=".c", prefix=tempfile_prefix)
        os.fdopen(fid, "w").write(source)
        command = COMPILE%{"source":filename, "output":dll}
        print("Compile command: "+command)
        status = os.system(command)
        if status != 0 or not os.path.exists(dll):
            raise RuntimeError("compile failed.  File is in %r"%filename)
        else:
            ## comment the following to keep the generated c file
            os.unlink(filename)
            #print("saving compiled file in %r"%filename)
    return dll


def load_dll(source, model_info, dtype="double"):
    """
    Create and load a dll corresponding to the source, info pair returned
    from :func:`sasmodels.generate.make` compiled for the target precision.

    See :func:`make_dll` for details on controlling the dll path and the
    allowed floating point precision.
    """
    filename = make_dll(source, model_info, dtype=dtype)
    return DllModel(filename, model_info, dtype=dtype)

class DllModel(object):
    """
    ctypes wrapper for a single model.

    *source* and *model_info* are the model source and interface as returned
    from :func:`gen.make`.

    *dtype* is the desired model precision.  Any numpy dtype for single
    or double precision floats will do, such as 'f', 'float32' or 'single'
    for single and 'd', 'float64' or 'double' for double.  Double precision
    is an optional extension which may not be available on all devices.

    Call :meth:`release` when done with the kernel.
    """
    def __init__(self, dllpath, model_info, dtype=generate.F32):
        self.info = model_info
        self.dllpath = dllpath
        self.dll = None
        self.dtype = np.dtype(dtype)

    def _load_dll(self):
        #print("dll", self.dllpath)
        try:
            self.dll = ct.CDLL(self.dllpath)
        except Exception as exc:
            annotate_exception(exc, "while loading "+self.dllpath)
            raise

        fp = (c_float if self.dtype == generate.F32
              else c_double if self.dtype == generate.F64
              else c_longdouble)

        # int, int, int, int*, double*, double*, double*, double*, double*, double
        argtypes = [c_int32]*3 + [c_void_p]*5 + [fp]
        self.Iq = self.dll[generate.kernel_name(self.info, False)]
        self.Iqxy = self.dll[generate.kernel_name(self.info, True)]
        self.Iq.argtypes = argtypes
        self.Iqxy.argtypes = argtypes

    def __getstate__(self):
        return self.info, self.dllpath

    def __setstate__(self, state):
        self.info, self.dllpath = state
        self.dll = None

    def __call__(self, q_vectors):
        q_input = PyInput(q_vectors, self.dtype)
        if self.dll is None: self._load_dll()
        kernel = self.Iqxy if q_input.is_2d else self.Iq
        return DllKernel(kernel, self.info, q_input)

    def release(self):
        """
        Release any resources associated with the model.
        """
        pass # TODO: should release the dll


class DllKernel(object):
    """
    Callable SAS kernel.

    *kernel* is the c function to call.

    *model_info* is the module information

    *q_input* is the DllInput q vectors at which the kernel should be
    evaluated.

    The resulting call method takes the *pars*, a list of values for
    the fixed parameters to the kernel, and *pd_pars*, a list of (value, weight)
    vectors for the polydisperse parameters.  *cutoff* determines the
    integration limits: any points with combined weight less than *cutoff*
    will not be calculated.

    Call :meth:`release` when done with the kernel instance.
    """
    def __init__(self, kernel, model_info, q_input):
        self.info = model_info
        self.q_input = q_input
        self.dtype = q_input.dtype
        self.kernel = kernel
        self.res = np.empty(q_input.nq+3, q_input.dtype)
        self.dim = '2d' if q_input.is_2d else '1d'

        # In dll kernel, but not in opencl kernel
        self.p_res = self.res.ctypes.data

    def __call__(self, details, weights, values, cutoff):
        real = (np.float32 if self.q_input.dtype == generate.F32
                else np.float64 if self.q_input.dtype == generate.F64
                else np.float128)
        if details.dtype != np.int32 or weights.dtype != real or values.dtype != real:
            raise TypeError("numeric type does not match the kernel type")
        #details = np.asarray(details, dtype='int32')
        #weights = np.asarray(weights, dtype=real)
        #values = np.asarray(values, dtype=real)
        #TODO: How can I access max_pd and is this the way to do it?
        #max_pd = model_info['max_pd']
        max_pd = 1
        args = [
            self.q_input.nq, # nq
            #TODO: pd_start will need to be changed
            0, # pd_start
            details[4*max_pd-1], # pd_stop pd_stride[MAX_PD]
            details.ctypes.data, # problem
            weights.ctypes.data,  # weights
            values.ctypes.data,  #pars
            self.q_input.q_pointers[0], #q
            self.p_res,   # results
            real(cutoff), # cutoff
            ]
        self.kernel(*args)
        return self.res[:-3]

    def release(self):
        """
        Release any resources associated with the kernel.
        """
        pass
