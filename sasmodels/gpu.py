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
import numpy as np
import pyopencl as cl
from pyopencl import mem_flags as mf

from . import gen

F32_DEFS = """\
#define REAL(x) (x##f)
#define real float
"""

F64_DEFS = """\
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define REAL(x) (x)
#define real double
"""

# The max loops number is limited by the amount of local memory available
# on the device.  You don't want to make this value too big because it will
# waste resources, nor too small because it may interfere with users trying
# to do their polydispersity calculations.  A value of 1024 should be much
# larger than necessary given that cost grows as npts^k where k is the number
# of polydisperse parameters.
MAX_LOOPS = 2048

def load_model(kernel_module, dtype="single"):
    """
    Load the OpenCL model defined by *kernel_module*.

    Access to the OpenCL device is delayed until the kernel is called
    so models can be defined without using too many resources.
    """
    source, info = gen.make(kernel_module)
    ## for debugging, save source to a .cl file, edit it, and reload as model
    #open(info['name']+'.cl','w').write(source)
    #source = open(info['name']+'.cl','r').read()
    return GpuModel(source, info, dtype)

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
    if dtype==gen.F64 and not all(has_double(d) for d in context.devices):
        raise RuntimeError("Double precision not supported for devices")

    header = F64_DEFS if dtype == gen.F64 else F32_DEFS
    # Note: USE_SINCOS makes the intel cpu slower under opencl
    if context.devices[0].type == cl.device_type.GPU:
        header += "#define USE_SINCOS\n"
    program  = cl.Program(context, header+source).build()
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
        self.context = cl.create_some_context()
        self.queues = [cl.CommandQueue(self.context, d)
                       for d in self.context.devices]
        self.boundary = max(d.min_data_type_align_size
                            for d in self.context.devices)
        self.has_double = all(has_double(d) for d in self.context.devices)
        self.compiled = {}

    def compile_program(self, name, source, dtype):
        if name not in self.compiled:
            #print "compiling",name
            self.compiled[name] = compile_model(self.context, source, dtype)
        return self.compiled[name]

    def release_program(self, name):
        if name in self.compiled:
            self.compiled[name].release()
            del self.compiled[name]

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
    def __init__(self, source, info, dtype=gen.F32):
        self.info = info
        self.source = source
        self.dtype = dtype
        self.program = None # delay program creation

    def __getstate__(self):
        state = self.__dict__.copy()
        state['program'] = None
        return state

    def __setstate__(self, state):
        self.__dict__ = state.copy()

    def __call__(self, input):
        if self.dtype != input.dtype:
            raise TypeError("data and kernel have different types")
        if self.program is None:
            self.program = environment().compile_program(self.info['name'],self.source, self.dtype)
        kernel_name = gen.kernel_name(self.info, input.is_2D)
        kernel = getattr(self.program, kernel_name)
        return GpuKernel(kernel, self.info, input)

    def release(self):
        if self.program is not None:
            environment().release_program(self.info['name'])
            self.program = None

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
    def __init__(self, q_vectors, dtype=gen.F32):
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
    """
    Callable SAS kernel.

    *kernel* is the GpuKernel object to call.

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
        self.input = input
        self.kernel = kernel
        self.info = info
        self.res = np.empty(input.nq, input.dtype)
        dim = '2d' if input.is_2D else '1d'
        self.fixed_pars = info['partype']['fixed-'+dim]
        self.pd_pars = info['partype']['pd-'+dim]

        # Inputs and outputs for each kernel call
        # Note: res may be shorter than res_b if global_size != nq
        env = environment()
        self.loops_b = [cl.Buffer(env.context, mf.READ_WRITE,
                                  2*MAX_LOOPS*input.dtype.itemsize)
                        for _ in env.queues]
        self.res_b = [cl.Buffer(env.context, mf.READ_WRITE,
                                input.global_size[0]*input.dtype.itemsize)
                      for _ in env.queues]


    def __call__(self, pars, pd_pars, cutoff=1e-5):
        real = np.float32 if self.input.dtype == gen.F32 else np.float64
        fixed = [real(p) for p in pars]
        cutoff = real(cutoff)
        loops = np.hstack(pd_pars)
        loops = np.ascontiguousarray(loops.T, self.input.dtype).flatten()
        Nloops = [np.uint32(len(p[0])) for p in pd_pars]
        #print "loops",Nloops, loops

        #import sys; print >>sys.stderr,"opencl eval",pars
        #print "opencl eval",pars
        if len(loops) > 2*MAX_LOOPS:
            raise ValueError("too many polydispersity points")
        device_num = 0
        res_bi = self.res_b[device_num]
        queuei = environment().queues[device_num]
        loops_bi = self.loops_b[device_num]
        loops_l = cl.LocalMemory(len(loops.data))
        cl.enqueue_copy(queuei, loops_bi, loops)
        #ctx = environment().context
        #loops_bi = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=loops)
        args = self.input.q_buffers + [res_bi,loops_bi,loops_l,cutoff] + fixed + Nloops
        self.kernel(queuei, self.input.global_size, None, *args)
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
