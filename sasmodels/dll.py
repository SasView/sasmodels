"""
C types wrapper for sasview models.
"""

import ctypes as ct
from ctypes import c_void_p, c_int, c_double

import numpy as np

from . import gen

from .gen import F32, F64

IQ_ARGS = [c_void_p, c_void_p, c_int, c_void_p, c_double]
IQXY_ARGS = [c_void_p, c_void_p, c_void_p, c_int, c_void_p, c_double]

class DllModel(object):
    """
    ctypes wrapper for a single model.

    *source* and *meta* are the model source and interface as returned
    from :func:`gen.make`.

    *dtype* is the desired model precision.  Any numpy dtype for single
    or double precision floats will do, such as 'f', 'float32' or 'single'
    for single and 'd', 'float64' or 'double' for double.  Double precision
    is an optional extension which may not be available on all devices.
    """
    def __init__(self, dllpath, meta):
        self.meta = meta
        self.dll = ct.CDLL(dllpath)
        self.Iq = self.dll[gen.kernel_name(self.meta, False)]
        self.Iqxy = self.dll[gen.kernel_name(self.meta, True)]


        self.PARS = dict((p[0],p[2]) for p in meta['parameters'])
        self.PD_PARS = [p[0] for p in meta['parameters'] if p[4] != ""]

        # Determine the set of fixed and polydisperse parameters
        Nfixed = len([p[0] for p in meta['parameters'] if p[4] == ""])
        N1D = len([p for p in meta['parameters'] if p[4]=="volume"])
        N2D = len([p for p in meta['parameters'] if p[4]!=""])
        self.Iq.argtypes = IQ_ARGS + [c_double]*Nfixed + [c_int]*N1D
        self.Iqxy.argtypes = IQXY_ARGS + [c_double]*Nfixed + [c_int]*N2D

    def __call__(self, input, cutoff=1e-5):
        kernel = self.Iqxy if input.is_2D else self.Iq
        return DllKernel(kernel, self.meta, input, cutoff)

    def make_input(self, q_vectors):
        """
        Make q input vectors available to the model.

        This only needs to be done once for all models that operate on the
        same input.  So for example, if you are adding two different models
        together to compare to a data set, then only one model needs to
        needs to call make_input, so long as the models have the same dtype.
        """
        return DllInput(q_vectors)


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
    def __init__(self, kernel, meta, input, cutoff):
        self.cutoff = cutoff
        self.input = input
        self.kernel = kernel
        self.meta = meta

        self.res = np.empty(input.nq, input.dtype)
        self.p_res = self.res.ctypes.data

        # Determine the set of fixed and polydisperse parameters
        self.fixed_pars = [p[0] for p in meta['parameters'] if p[4] == ""]
        self.pd_pars = [p for p in meta['parameters']
               if p[4]=="volume" or (p[4]=="orientation" and input.is_2D)]

    def eval(self, pars):
        fixed, loops, loop_n = \
            gen.kernel_pars(pars, self.meta, self.input.is_2D, dtype=self.input.dtype)
        real = np.float32 if self.input.dtype == F32 else np.float64
        nq = c_int(self.input.nq)
        cutoff = real(self.cutoff)

        p_loops = loops.ctypes.data
        pars = self.input.q_pointers + [self.p_res, nq, p_loops, cutoff] + fixed + loop_n
        #print pars
        self.kernel(*pars)

        return self.res

    def release(self):
        pass

    def __del__(self):
        self.release()
