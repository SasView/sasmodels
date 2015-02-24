"""
C types wrapper for sasview models.
"""
import sys
import os
import tempfile
import ctypes as ct
from ctypes import c_void_p, c_int, c_double

import numpy as np

from . import generate
from .kernelpy import PyInput, PyModel

from .generate import F32, F64
# Compiler platform details
if sys.platform == 'darwin':
    #COMPILE = "gcc-mp-4.7 -shared -fPIC -std=c99 -fopenmp -O2 -Wall %s -o %s -lm -lgomp"
    COMPILE = "gcc -shared -fPIC -std=c99 -O2 -Wall %(source)s -o %(output)s -lm"
elif os.name == 'nt':
    # make sure vcvarsall.bat is called first in order to set compiler, headers, lib paths, etc.
    if "VCINSTALLDIR" in os.environ:
        # MSVC compiler is available, so use it.
        # TODO: remove intermediate OBJ file created in the directory
        # TODO: maybe don't use randomized name for the c file
        COMPILE = "cl /nologo /Ox /MD /W3 /GS- /DNDEBUG /Tp%(source)s /openmp /link /DLL /INCREMENTAL:NO /MANIFEST /OUT:%(output)s"
        # Can't find VCOMP90.DLL (don't know why), so remove openmp support from windows compiler build
        #COMPILE = "cl /nologo /Ox /MD /W3 /GS- /DNDEBUG /Tp%(source)s /link /DLL /INCREMENTAL:NO /MANIFEST /OUT:%(output)s"
    else:
        #COMPILE = "gcc -shared -fPIC -std=c99 -fopenmp -O2 -Wall %(source)s -o %(output)s -lm"
        COMPILE = "gcc -shared -fPIC -std=c99 -O2 -Wall %(source)s -o %(output)s -lm"
else:
    COMPILE = "cc -shared -fPIC -std=c99 -fopenmp -O2 -Wall %(source)s -o %(output)s -lm"

DLL_PATH = tempfile.gettempdir()


def dll_path(info):
    """
    Path to the compiled model defined by *info*.
    """
    from os.path import join as joinpath, split as splitpath, splitext
    basename = splitext(splitpath(info['filename'])[1])[0]
    return joinpath(DLL_PATH, basename+'.so')


def load_model(kernel_module, dtype=None):
    """
    Load the compiled model defined by *kernel_module*.

    Recompile if any files are newer than the model file.

    *dtype* is ignored.  Compiled files are always double.

    The DLL is not loaded until the kernel is called so models an
    be defined without using too many resources.
    """
    source, info = generate.make(kernel_module)
    if callable(info.get('Iq',None)):
        return PyModel(info)
    source_files = generate.sources(info) + [info['filename']]
    newest = max(os.path.getmtime(f) for f in source_files)
    dllpath = dll_path(info)
    if not os.path.exists(dllpath) or os.path.getmtime(dllpath)<newest:
        # Replace with a proper temp file
        fid, filename = tempfile.mkstemp(suffix=".c",prefix="sas_"+info['name'])
        os.fdopen(fid,"w").write(source)
        command = COMPILE%{"source":filename, "output":dllpath}
        print "Compile command:",command
        status = os.system(command)
        if status != 0:
            raise RuntimeError("compile failed.  File is in %r"%filename)
        else:
            ## uncomment the following to keep the generated c file
            #os.unlink(filename); print "saving compiled file in %r"%filename
            pass
    return DllModel(dllpath, info)


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
    def __init__(self, dllpath, info):
        self.info = info
        self.dllpath = dllpath
        self.dll = None

    def _load_dll(self):
        Nfixed1d = len(self.info['partype']['fixed-1d'])
        Nfixed2d = len(self.info['partype']['fixed-2d'])
        Npd1d = len(self.info['partype']['pd-1d'])
        Npd2d = len(self.info['partype']['pd-2d'])

        #print "dll",self.dllpath
        self.dll = ct.CDLL(self.dllpath)

        pd_args_1d = [c_void_p, c_double] + [c_int]*Npd1d if Npd1d else []
        pd_args_2d= [c_void_p, c_double] + [c_int]*Npd2d if Npd2d else []
        self.Iq = self.dll[generate.kernel_name(self.info, False)]
        self.Iq.argtypes = IQ_ARGS + pd_args_1d + [c_double]*Nfixed1d

        self.Iqxy = self.dll[generate.kernel_name(self.info, True)]
        self.Iqxy.argtypes = IQXY_ARGS + pd_args_2d + [c_double]*Nfixed2d

    def __getstate__(self):
        return {'info': self.info, 'dllpath': self.dllpath, 'dll': None}

    def __setstate__(self, state):
        self.__dict__ = state

    def __call__(self, input):
        if self.dll is None: self._load_dll()
        kernel = self.Iqxy if input.is_2D else self.Iq
        return DllKernel(kernel, self.info, input)

    def make_input(self, q_vectors):
        """
        Make q input vectors available to the model.

        Note that each model needs its own q vector even if the case of
        mixture models because some models may be OpenCL, some may be
        ctypes and some may be pure python.
        """
        return PyInput(q_vectors, dtype=F64)

    def release(self):
        pass # TODO: should release the dll


class DllKernel(object):
    """
    Callable SAS kernel.

    *kernel* is the c function to call.

    *info* is the module information

    *input* is the DllInput q vectors at which the kernel should be
    evaluated.

    The resulting call method takes the *pars*, a list of values for
    the fixed parameters to the kernel, and *pd_pars*, a list of (value,weight)
    vectors for the polydisperse parameters.  *cutoff* determines the
    integration limits: any points with combined weight less than *cutoff*
    will not be calculated.

    Call :meth:`release` when done with the kernel instance.
    """
    def __init__(self, kernel, info, input):
        self.info = info
        self.input = input
        self.kernel = kernel
        self.res = np.empty(input.nq, input.dtype)
        dim = '2d' if input.is_2D else '1d'
        self.fixed_pars = info['partype']['fixed-'+dim]
        self.pd_pars = info['partype']['pd-'+dim]

        # In dll kernel, but not in opencl kernel
        self.p_res = self.res.ctypes.data

    def __call__(self, fixed_pars, pd_pars, cutoff):
        real = np.float32 if self.input.dtype == F32 else np.float64

        nq = c_int(self.input.nq)
        if pd_pars:
            cutoff = real(cutoff)
            loops_N = [np.uint32(len(p[0])) for p in pd_pars]
            loops = np.hstack(pd_pars)
            loops = np.ascontiguousarray(loops.T, self.input.dtype).flatten()
            p_loops = loops.ctypes.data
            dispersed = [p_loops, cutoff] + loops_N
        else:
            dispersed = []
        fixed = [real(p) for p in fixed_pars]
        args = self.input.q_pointers + [self.p_res, nq] + dispersed + fixed
        #print pars
        self.kernel(*args)

        return self.res

    def release(self):
        pass
