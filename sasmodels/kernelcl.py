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
import os
import warnings

import numpy as np

try:
    import pyopencl as cl
    # Ask OpenCL for the default context so that we know that one exists
    cl.create_some_context(interactive=False)
except Exception, exc:
    warnings.warn(str(exc))
    raise RuntimeError("OpenCL not available")

from pyopencl import mem_flags as mf

from . import generate

F64_DEFS = """\
#ifdef cl_khr_fp64
#  pragma OPENCL EXTENSION cl_khr_fp64: enable
#endif
"""

# The max loops number is limited by the amount of local memory available
# on the device.  You don't want to make this value too big because it will
# waste resources, nor too small because it may interfere with users trying
# to do their polydispersity calculations.  A value of 1024 should be much
# larger than necessary given that cost grows as npts^k where k is the number
# of polydisperse parameters.
MAX_LOOPS = 2048


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

def get_warp(kernel, queue):
    """
    Return the size of an execution batch for *kernel* running on *queue*.
    """
    return kernel.get_work_group_info(
        cl.kernel_work_group_info.PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
        queue.device)

def _stretch_input(vector, dtype, extra=1e-3, boundary=32):
    """
    Stretch an input vector to the correct boundary.

    Performance on the kernels can drop by a factor of two or more if the
    number of values to compute does not fall on a nice power of two
    boundary.   The trailing additional vector elements are given a
    value of *extra*, and so f(*extra*) will be computed for each of
    them.  The returned array will thus be a subset of the computed array.

    *boundary* should be a power of 2 which is at least 32 for good
    performance on current platforms (as of Jan 2015).  It should
    probably be the max of get_warp(kernel,queue) and
    device.min_data_type_align_size//4.
    """
    remainder = vector.size % boundary
    if remainder != 0:
        size = vector.size + (boundary - remainder)
        vector = np.hstack((vector, [extra] * (size - vector.size)))
    return np.ascontiguousarray(vector, dtype=dtype)


def compile_model(context, source, dtype):
    """
    Build a model to run on the gpu.

    Returns the compiled program and its type.  The returned type will
    be float32 even if the desired type is float64 if any of the
    devices in the context do not support the cl_khr_fp64 extension.
    """
    dtype = np.dtype(dtype)
    if dtype == generate.F64 and not all(has_double(d) for d in context.devices):
        raise RuntimeError("Double precision not supported for devices")

    header = F64_DEFS if dtype == generate.F64 else ""
    if dtype == generate.F32:
        source = generate.use_single(source)
    # Note: USE_SINCOS makes the intel cpu slower under opencl
    if context.devices[0].type == cl.device_type.GPU:
        header += "#define USE_SINCOS\n"
    program = cl.Program(context, header + source).build()
    return program


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
        # find gpu context
        #self.context = cl.create_some_context()

        self.context = None
        if 'PYOPENCL_CTX' in os.environ:
            self._create_some_context()

        if not self.context:
            self.context = _get_default_context()

        # Byte boundary for data alignment
        #self.data_boundary = max(d.min_data_type_align_size
        #                         for d in self.context.devices)
        self.queues = [cl.CommandQueue(self.context, d)
                       for d in self.context.devices]
        self.has_double = all(has_double(d) for d in self.context.devices)
        self.compiled = {}

    def _create_some_context(self):
        try:
            self.context = cl.create_some_context(interactive=False)
        except Exception, exc:
            warnings.warn(str(exc))
            warnings.warn("pyopencl.create_some_context() failed")
            warnings.warn("the environment variable 'PYOPENCL_CTX' might not be set correctly")

    def compile_program(self, name, source, dtype):
        if name not in self.compiled:
            #print "compiling",name
            self.compiled[name] = compile_model(self.context, source, dtype)
        return self.compiled[name]

    def release_program(self, name):
        if name in self.compiled:
            self.compiled[name].release()
            del self.compiled[name]

def _get_default_context():
    default = None
    for platform in cl.get_platforms():
        for device in platform.get_devices():
            if device.type == cl.device_type.GPU:
                return cl.Context([device])
            if default is None:
                default = device

    if not default:
        raise RuntimeError("OpenCL device not found")

    return cl.Context([default])


class GpuModel(object):
    """
    GPU wrapper for a single model.

    *source* and *info* are the model source and interface as returned
    from :func:`gen.make`.

    *dtype* is the desired model precision.  Any numpy dtype for single
    or double precision floats will do, such as 'f', 'float32' or 'single'
    for single and 'd', 'float64' or 'double' for double.  Double precision
    is an optional extension which may not be available on all devices.
    """
    def __init__(self, source, info, dtype=generate.F32):
        self.info = info
        self.source = source
        self.dtype = np.dtype(dtype)
        self.program = None # delay program creation

    def __getstate__(self):
        state = self.__dict__.copy()
        state['program'] = None
        return state

    def __setstate__(self, state):
        self.__dict__ = state.copy()

    def __call__(self, q_input):
        if self.dtype != q_input.dtype:
            raise TypeError("data is %s kernel is %s" % (q_input.dtype, self.dtype))
        if self.program is None:
            compiler = environment().compile_program
            self.program = compiler(self.info['name'], self.source, self.dtype)
        kernel_name = generate.kernel_name(self.info, q_input.is_2D)
        kernel = getattr(self.program, kernel_name)
        return GpuKernel(kernel, self.info, q_input)

    def release(self):
        if self.program is not None:
            environment().release_program(self.info['name'])
            self.program = None

    def make_input(self, q_vectors):
        """
        Make q input vectors available to the model.

        Note that each model needs its own q vector even if the case of
        mixture models because some models may be OpenCL, some may be
        ctypes and some may be pure python.
        """
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
    def __init__(self, q_vectors, dtype=generate.F32):
        env = environment()
        self.nq = q_vectors[0].size
        self.dtype = np.dtype(dtype)
        self.is_2D = (len(q_vectors) == 2)
        # TODO: stretch input based on get_warp()
        # not doing it now since warp depends on kernel, which is not known
        # at this point, so instead using 32, which is good on the set of
        # architectures tested so far.
        self.q_vectors = [_stretch_input(q, self.dtype, 32) for q in q_vectors]
        self.q_buffers = [
            cl.Buffer(env.context, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=q)
            for q in self.q_vectors
        ]
        self.global_size = [self.q_vectors[0].size]

    def release(self):
        for b in self.q_buffers:
            b.release()
        self.q_buffers = []

class GpuKernel(object):
    """
    Callable SAS kernel.

    *kernel* is the GpuKernel object to call.

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
        self.q_input = q_input
        self.kernel = kernel
        self.info = info
        self.res = np.empty(q_input.nq, q_input.dtype)
        dim = '2d' if q_input.is_2D else '1d'
        self.fixed_pars = info['partype']['fixed-' + dim]
        self.pd_pars = info['partype']['pd-' + dim]

        # Inputs and outputs for each kernel call
        # Note: res may be shorter than res_b if global_size != nq
        env = environment()
        self.loops_b = [cl.Buffer(env.context, mf.READ_WRITE,
                                  2 * MAX_LOOPS * q_input.dtype.itemsize)
                        for _ in env.queues]
        self.res_b = [cl.Buffer(env.context, mf.READ_WRITE,
                                q_input.global_size[0] * q_input.dtype.itemsize)
                      for _ in env.queues]


    def __call__(self, fixed_pars, pd_pars, cutoff=1e-5):
        real = np.float32 if self.q_input.dtype == generate.F32 else np.float64

        device_num = 0
        queuei = environment().queues[device_num]
        res_bi = self.res_b[device_num]
        nq = np.uint32(self.q_input.nq)
        if pd_pars:
            cutoff = real(cutoff)
            loops_N = [np.uint32(len(p[0])) for p in pd_pars]
            loops = np.hstack(pd_pars) \
                if pd_pars else np.empty(0, dtype=self.q_input.dtype)
            loops = np.ascontiguousarray(loops.T, self.q_input.dtype).flatten()
            #print "loops",Nloops, loops

            #import sys; print >>sys.stderr,"opencl eval",pars
            #print "opencl eval",pars
            if len(loops) > 2 * MAX_LOOPS:
                raise ValueError("too many polydispersity points")

            loops_bi = self.loops_b[device_num]
            cl.enqueue_copy(queuei, loops_bi, loops)
            loops_l = cl.LocalMemory(len(loops.data))
            #ctx = environment().context
            #loops_bi = cl.Buffer(ctx, mf.READ_ONLY|mf.COPY_HOST_PTR, hostbuf=loops)
            dispersed = [loops_bi, loops_l, cutoff] + loops_N
        else:
            dispersed = []
        fixed = [real(p) for p in fixed_pars]
        args = self.q_input.q_buffers + [res_bi, nq] + dispersed + fixed
        self.kernel(queuei, self.q_input.global_size, None, *args)
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
