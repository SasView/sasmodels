"""
C types wrapper for sasview models.
"""
import sys
import os
import ctypes as ct
from ctypes import c_void_p, c_int, c_double

import numpy as np

from . import gen

from .gen import F32, F64
# Compiler platform details
if sys.platform == 'darwin':
    COMPILE = "gcc-mp-4.7 -shared -fPIC -std=c99 -fopenmp -O2 -Wall %s -o %s -lm -lgomp"
elif os.name == 'nt':
    COMPILE = "gcc -shared -fPIC -std=c99 -fopenmp -O2 -Wall %s -o %s -lm"
else:
    COMPILE = "cc -shared -fPIC -std=c99 -fopenmp -O2 -Wall %s -o %s -lm"
DLL_PATH = "/tmp"


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
    import tempfile

    source, info = gen.make(kernel_module)
    source_files = gen.sources(info) + [info['filename']]
    newest = max(os.path.getmtime(f) for f in source_files)
    dllpath = dll_path(info)
    if not os.path.exists(dllpath) or os.path.getmtime(dllpath)<newest:
        # Replace with a proper temp file
        fid, filename = tempfile.mkstemp(suffix=".c",prefix="sas_"+info['name'])
        os.fdopen(fid,"w").write(source)
        status = os.system(COMPILE%(filename, dllpath))
        if status != 0:
            print "compile failed.  File is in %r"%filename
        else:
            ## uncomment the following to keep the generated c file
            #os.unlink(filename); print "saving compiled file in %r"%filename
            pass
    return DllModel(dllpath, info)


IQ_ARGS = [c_void_p, c_void_p, c_int, c_void_p, c_double]
IQXY_ARGS = [c_void_p, c_void_p, c_void_p, c_int, c_void_p, c_double]

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

        self.dll = ct.CDLL(self.dllpath)

        self.Iq = self.dll[gen.kernel_name(self.info, False)]
        self.Iq.argtypes = IQ_ARGS + [c_double]*Nfixed1d + [c_int]*Npd1d

        self.Iqxy = self.dll[gen.kernel_name(self.info, True)]
        self.Iqxy.argtypes = IQXY_ARGS + [c_double]*Nfixed2d + [c_int]*Npd2d

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

        This only needs to be done once for all models that operate on the
        same input.  So for example, if you are adding two different models
        together to compare to a data set, then only one model needs to
        needs to call make_input, so long as the models have the same dtype.
        """
        return DllInput(q_vectors)

    def release(self):
        pass # TODO: should release the dll


class DllInput(object):
    """
    Make q data available to the gpu.

    *q_vectors* is a list of q vectors, which will be *[q]* for 1-D data,
    and *[qx, qy]* for 2-D data.  Internally, the vectors will be reallocated
    to get the best performance on OpenCL, which may involve shifting and
    stretching the array to better match the memory architecture.  Additional
    points will be evaluated with *q=1e-3*.

    *dtype* is the data type for the q vectors. The data type should be
    set to match that of the kernel, which is an attribute of
    :class:`GpuProgram`.  Note that not all kernels support double
    precision, so even if the program was created for double precision,
    the *GpuProgram.dtype* may be single precision.

    Call :meth:`release` when complete.  Even if not called directly, the
    buffer will be released when the data object is freed.
    """
    def __init__(self, q_vectors):
        self.nq = q_vectors[0].size
        self.dtype = np.dtype('double')
        self.is_2D = (len(q_vectors) == 2)
        self.q_vectors = [np.ascontiguousarray(q,self.dtype) for q in q_vectors]
        self.q_pointers = [q.ctypes.data for q in q_vectors]

    def release(self):
        self.q_vectors = []

class DllKernel(object):
    """
    Callable SAS kernel.

    *kernel* is the DllKernel object to call.

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

    def __call__(self, pars, pd_pars, cutoff):
        real = np.float32 if self.input.dtype == F32 else np.float64
        fixed = [real(p) for p in pars]
        cutoff = real(cutoff)
        loops = np.hstack(pd_pars)
        loops = np.ascontiguousarray(loops.T, self.input.dtype).flatten()
        loops_N = [np.uint32(len(p[0])) for p in pd_pars]

        nq = c_int(self.input.nq)
        p_loops = loops.ctypes.data
        args = self.input.q_pointers + [self.p_res, nq, p_loops, cutoff] + fixed + loops_N
        #print pars
        self.kernel(*args)

        return self.res

    def release(self):
        pass
