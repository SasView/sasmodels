"""
GPU support through OpenCL

There should be a single GPU environment running on the system.  This
environment is constructed on the first call to :func:`env`, and the
same environment is returned on each call.

After retrieving the environment, the next step is to create the kernel.
This is done with a call to :meth:`GpuEnvironment.make_kernel`, which
returns the type of data used by the kernel.

Next a :class:`GpuData` object should be created with the correct kind
of data.  This data object can be used by multiple kernels, for example,
if the target model is a weighted sum of multiple kernels.  The data
should include any extra evaluation points required to compute the proper
data smearing.  This need not match the square grid for 2D data if there
is an index saying which q points are active.

Together the GpuData, the program, and a device form a :class:`GpuKernel`.
This kernel is used during fitting, receiving new sets of parameters and
evaluating them.  The output value is stored in an output buffer on the
devices, where it can be combined with other structure factors and form
factors and have instrumental resolution effects applied.


"""
import warnings

import numpy as np
import pyopencl as cl
from pyopencl import mem_flags as mf

from . import gen

from .gen import F32, F64

F32_DEFS = """\
#define REAL(x) (x##f)
#define real float
"""

F64_DEFS = """\
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define REAL(x) (x)
#define real double
"""

ENV = None
def environment():
    """
    Returns a singleton :class:`GpuEnvironment`.

    This provides an OpenCL context and one queue per device.
    """
    global ENV
    if ENV is None:
        ENV = GpuEnvironment()
    return ENV

def has_double(device):
    """
    Return true if device supports double precision.
    """
    return "cl_khr_fp64" in device.extensions


def _stretch_input(vector, dtype, extra=1e-3, boundary=128):
    """
    Stretch an input vector to the correct boundary.

    Performance on the kernels can drop by a factor of two or more if the
    number of values to compute does not fall on a nice power of two
    boundary.  A good choice for the boundary value is the
    min_data_type_align_size property of the OpenCL device.  The usual
    value of 128 gives a working size as a multiple of 32.  The trailing
    additional vector elements are given a value of *extra*, and so
    f(*extra*) will be computed for each of them.  The returned array
    will thus be a subset of the computed array.
    """
    boundary // dtype.itemsize
    remainder = vector.size%boundary
    size = vector.size + (boundary - remainder if remainder != 0 else 0)
    if size != vector.size:
        vector = np.hstack((vector, [extra]*(size-vector.size)))
    return np.ascontiguousarray(vector, dtype=dtype)


def compile_model(context, source, dtype):
    """
    Build a model to run on the gpu.

    Returns the compiled program and its type.  The returned type will
    be float32 even if the desired type is float64 if any of the
    devices in the context do not support the cl_khr_fp64 extension.
    """
    dtype = np.dtype(dtype)
    if dtype==F64 and not all(has_double(d) for d in context.devices):
        warnings.warn(RuntimeWarning("Double precision not support; using single precision instead"))
        dtype = F32

    header = F64_DEFS if dtype == F64 else F32_DEFS
    # Note: USE_SINCOS makes the intel cpu slower under opencl
    if context.devices[0].type == cl.device_type.GPU:
        header += "#define USE_SINCOS\n"
    program  = cl.Program(context, header+source).build()
    return program, dtype


def make_result(self, size):
    self.res = np.empty(size, dtype=self.dtype)
    self.res_b = cl.Buffer(self.program.context, mf.READ_WRITE, self.res.nbytes)
    return self.res, self.res_b


# for now, this returns one device in the context
# TODO: create a context that contains all devices on all platforms
class GpuEnvironment(object):
    """
    GPU context, with possibly many devices, and one queue per device.
    """
    def __init__(self):
        self.context = cl.create_some_context()
        self.queues = [cl.CommandQueue(self.context, d)
                       for d in self.context.devices]
        self.boundary = max(d.min_data_type_align_size
                            for d in self.context.devices)

class GpuModel(object):
    """
    GPU wrapper for a single model.

    *source* and *meta* are the model source and interface as returned
    from :func:`gen.make`.

    *dtype* is the desired model precision.  Any numpy dtype for single
    or double precision floats will do, such as 'f', 'float32' or 'single'
    for single and 'd', 'float64' or 'double' for double.  Double precision
    is an optional extension which may not be available on all devices.
    """
    def __init__(self, source, meta, dtype=F32):
        context = environment().context
        self.meta = meta
        self.program, self.dtype = compile_model(context, source, dtype)

        #TODO: need better interface
        self.PARS = dict((p[0],p[2]) for p in meta['parameters'])
        self.PD_PARS = [p[0] for p in meta['parameters'] if p[4] != ""]

    def __call__(self, input, cutoff=1e-5):
        if self.dtype != input.dtype:
            raise TypeError("data and kernel have different types")
        kernel_name = gen.kernel_name(self.meta, input.is_2D)
        kernel = getattr(self.program, kernel_name)
        return GpuKernel(kernel, self.meta, input, cutoff)

    def make_input(self, q_vectors):
        """
        Make q input vectors available to the model.

        This only needs to be done once for all models that operate on the
        same input.  So for example, if you are adding two different models
        together to compare to a data set, then only one model needs to
        needs to call make_input, so long as the models have the same dtype.
        """
        # Note: the weird interface, where make_input doesn't care which
        # model calls it, allows us to ask the model to define the data
        # and the caller does not need to know if it is opencl or ctypes.
        # The returned data object is opaque.
        return GpuInput(q_vectors, dtype=self.dtype)

# TODO: check that we don't need a destructor for buffers which go out of scope
class GpuInput(object):
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
    def __init__(self, q_vectors, dtype=F32):
        env = environment()
        self.nq = q_vectors[0].size
        self.dtype = np.dtype(dtype)
        self.is_2D = (len(q_vectors) == 2)
        self.q_vectors = [
            _stretch_input(q, self.dtype, boundary=env.boundary)
            for q in q_vectors
        ]
        self.q_buffers = [
            cl.Buffer(env.context,  mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=q)
            for q in self.q_vectors
        ]
        self.global_size = [self.q_vectors[0].size]

    def release(self):
        for b in self.q_buffers:
            b.release()
        self.q_buffers = []

class GpuKernel(object):
    def __init__(self, kernel, meta, input, cutoff):
        env = environment()

        self.cutoff = cutoff
        self.input = input
        self.kernel = kernel
        self.meta = meta

        # Inputs and outputs for each kernel call
        self.loops_b = [cl.Buffer(env.context, mf.READ_WRITE,
                                  1024*input.dtype.itemsize)
                        for _ in env.queues]
        self.res_b = [cl.Buffer(env.context, mf.READ_WRITE,
                                input.global_size[0]*input.dtype.itemsize)
                      for _ in env.queues]

        # Note: may be shorter than res_b if global_size != nq
        self.res = np.empty(input.nq, input.dtype)

        # Determine the set of fixed and polydisperse parameters
        self.fixed_pars = [p[0] for p in meta['parameters'] if p[4] == ""]
        self.pd_pars = [p for p in meta['parameters']
               if p[4]=="volume" or (p[4]=="orientation" and input.is_2D)]

    def eval(self, pars):
        device_num = 0
        res_bi = self.res_b[device_num]
        queuei = environment().queues[device_num]
        fixed, loops, loop_n = \
            gen.kernel_pars(pars, self.meta, self.input.is_2D, dtype=self.input.dtype)
        loops_l = cl.LocalMemory(len(loops.data))
        real = np.float32 if self.input.dtype == F32 else np.float64
        cutoff = real(self.cutoff)

        loops_bi = self.loops_b[device_num]
        cl.enqueue_copy(queuei, loops_bi, loops)
        #ctx = environment().context
        #loops_bi = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=loops)
        pars = self.input.q_buffers + [res_bi,loops_bi,loops_l,cutoff] + fixed + loop_n
        self.kernel(queuei, self.input.global_size, None, *pars)
        cl.enqueue_copy(queuei, self.res, res_bi)

        return self.res

    def release(self):
        for b in self.loops_b:
            b.release()
        self.loops_b = []
        for b in self.res_b:
            b.release()
        self.res_b = []

    def __del__(self):
        self.release()
