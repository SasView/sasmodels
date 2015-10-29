"""
C types wrapper for sasview models.

The global attribute *ALLOW_SINGLE_PRECISION_DLLS* should be set to *True* if
you wish to allow single precision floating point evaluation for the compiled
models, otherwise it defaults to *False*.
"""

import sys
import os
import tempfile
import ctypes as ct
from ctypes import c_void_p, c_int, c_double, c_float

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
        # MSVC compiler is available, so use it.
        # TODO: remove intermediate OBJ file created in the directory
        # TODO: maybe don't use randomized name for the c file
        #COMPILE = "cl /nologo /Ox /MD /W3 /GS- /DNDEBUG /Tp%(source)s /openmp /link /DLL /INCREMENTAL:NO /MANIFEST /OUT:%(output)s"
        # Can't find VCOMP90.DLL (don't know why), so remove openmp support
        # from windows compiler build
        COMPILE = "cl /nologo /Ox /MD /W3 /GS- /DNDEBUG /Tp%(source)s /link /DLL /INCREMENTAL:NO /MANIFEST /OUT:%(output)s"
    else:
        #COMPILE = "gcc -shared -fPIC -std=c99 -fopenmp -O2 -Wall %(source)s -o %(output)s -lm"
        COMPILE = "gcc -shared -fPIC -std=c99 -O2 -Wall %(source)s -o %(output)s -lm"
else:
    COMPILE = "cc -shared -fPIC -std=c99 -fopenmp -O2 -Wall %(source)s -o %(output)s -lm"

DLL_PATH = tempfile.gettempdir()

ALLOW_SINGLE_PRECISION_DLLS = False


def dll_path(info, dtype="double"):
    """
    Path to the compiled model defined by *info*.
    """
    from os.path import join as joinpath, split as splitpath, splitext
    basename = splitext(splitpath(info['filename'])[1])[0]
    if np.dtype(dtype) == generate.F32:
        basename += "32"
    return joinpath(DLL_PATH, basename+'.so')


def make_dll(source, info, dtype="double"):
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
    if not ALLOW_SINGLE_PRECISION_DLLS: dtype = "double"   # Force 64-bit dll
    dtype = np.dtype(dtype)

    if callable(info.get('Iq',None)):
        return PyModel(info)

    if dtype == generate.F32: # 32-bit dll
        source = generate.use_single(source)
        tempfile_prefix = 'sas_'+info['name']+'32_'
    else:
        tempfile_prefix = 'sas_'+info['name']+'_'

    source_files = generate.sources(info) + [info['filename']]
    dll= dll_path(info, dtype)
    newest = max(os.path.getmtime(f) for f in source_files)
    if not os.path.exists(dll) or os.path.getmtime(dll)<newest:
        # Replace with a proper temp file
        fid, filename = tempfile.mkstemp(suffix=".c",prefix=tempfile_prefix)
        os.fdopen(fid,"w").write(source)
        command = COMPILE%{"source":filename, "output":dll}
        print "Compile command:",command
        status = os.system(command)
        if status != 0 or not os.path.exists(dll):
            raise RuntimeError("compile failed.  File is in %r"%filename)
        else:
            ## uncomment the following to keep the generated c file
            os.unlink(filename); print "saving compiled file in %r"%filename
    return dll


def load_dll(source, info, dtype="double"):
    """
    Create and load a dll corresponding to the source,info pair returned
    from :func:`sasmodels.generate.make` compiled for the target precision.

    See :func:`make_dll` for details on controlling the dll path and the
    allowed floating point precision.
    """
    filename = make_dll(source, info, dtype=dtype)
    return DllModel(filename, info, dtype=dtype)


IQ_ARGS = [c_void_p, c_void_p, c_int]
IQXY_ARGS = [c_void_p, c_void_p, c_void_p, c_int]

class DllModel(object):
    """
    ctypes wrapper for a single model.

    *source* and *info* are the model source and interface as returned
    from :func:`gen.make`.

    *dtype* is the desired model precision.  Any numpy dtype for single
    or double precision floats will do, such as 'f', 'float32' or 'single'
    for single and 'd', 'float64' or 'double' for double.  Double precision
    is an optional extension which may not be available on all devices.

    Call :meth:`release` when done with the kernel.
    """
    def __init__(self, dllpath, info, dtype=generate.F32):
        self.info = info
        self.dllpath = dllpath
        self.dll = None
        self.dtype = np.dtype(dtype)

    def _load_dll(self):
        Nfixed1d = len(self.info['partype']['fixed-1d'])
        Nfixed2d = len(self.info['partype']['fixed-2d'])
        Npd1d = len(self.info['partype']['pd-1d'])
        Npd2d = len(self.info['partype']['pd-2d'])

        #print "dll",self.dllpath
        try:
            self.dll = ct.CDLL(self.dllpath)
        except Exception, exc:
            annotate_exception(exc, "while loading "+self.dllpath)
            raise

        fp = c_float if self.dtype == generate.F32 else c_double
        pd_args_1d = [c_void_p, fp] + [c_int]*Npd1d if Npd1d else []
        pd_args_2d= [c_void_p, fp] + [c_int]*Npd2d if Npd2d else []
        self.Iq = self.dll[generate.kernel_name(self.info, False)]
        self.Iq.argtypes = IQ_ARGS + pd_args_1d + [fp]*Nfixed1d

        self.Iqxy = self.dll[generate.kernel_name(self.info, True)]
        self.Iqxy.argtypes = IQXY_ARGS + pd_args_2d + [fp]*Nfixed2d

    def __getstate__(self):
        return {'info': self.info, 'dllpath': self.dllpath, 'dll': None}

    def __setstate__(self, state):
        self.__dict__ = state

    def __call__(self, q_input):
        if self.dtype != q_input.dtype:
            raise TypeError("data is %s kernel is %s" % (q_input.dtype, self.dtype))
        if self.dll is None: self._load_dll()
        kernel = self.Iqxy if q_input.is_2D else self.Iq
        return DllKernel(kernel, self.info, q_input)

    # pylint: disable=no-self-use
    def make_input(self, q_vectors):
        """
        Make q input vectors available to the model.

        Note that each model needs its own q vector even if the case of
        mixture models because some models may be OpenCL, some may be
        ctypes and some may be pure python.
        """
        return PyInput(q_vectors, dtype=self.dtype)

    def release(self):
        pass # TODO: should release the dll


class DllKernel(object):
    """
    Callable SAS kernel.

    *kernel* is the c function to call.

    *info* is the module information

    *q_input* is the DllInput q vectors at which the kernel should be
    evaluated.

    The resulting call method takes the *pars*, a list of values for
    the fixed parameters to the kernel, and *pd_pars*, a list of (value,weight)
    vectors for the polydisperse parameters.  *cutoff* determines the
    integration limits: any points with combined weight less than *cutoff*
    will not be calculated.

    Call :meth:`release` when done with the kernel instance.
    """
    def __init__(self, kernel, info, q_input):
        self.info = info
        self.q_input = q_input
        self.kernel = kernel
        self.res = np.empty(q_input.nq, q_input.dtype)
        dim = '2d' if q_input.is_2D else '1d'
        self.fixed_pars = info['partype']['fixed-'+dim]
        self.pd_pars = info['partype']['pd-'+dim]

        # In dll kernel, but not in opencl kernel
        self.p_res = self.res.ctypes.data

    def __call__(self, fixed_pars, pd_pars, cutoff):
        real = np.float32 if self.q_input.dtype == generate.F32 else np.float64

        nq = c_int(self.q_input.nq)
        if pd_pars:
            cutoff = real(cutoff)
            loops_N = [np.uint32(len(p[0])) for p in pd_pars]
            loops = np.hstack(pd_pars)
            loops = np.ascontiguousarray(loops.T, self.q_input.dtype).flatten()
            p_loops = loops.ctypes.data
            dispersed = [p_loops, cutoff] + loops_N
        else:
            dispersed = []
        fixed = [real(p) for p in fixed_pars]
        args = self.q_input.q_pointers + [self.p_res, nq] + dispersed + fixed
        #print pars
        self.kernel(*args)

        return self.res

    def release(self):
        pass
