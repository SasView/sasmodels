"""
GPU driver for C kernels

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

In order to use OpenCL for your models, you will need OpenCL drivers for
your machine.  These should be available from your graphics card vendor.
Intel provides OpenCL drivers for CPUs as well as their integrated HD
graphics chipsets.  AMD also provides drivers for Intel CPUs, but as of
this writing the performance is lacking compared to the Intel drivers.
NVidia combines drivers for CUDA and OpenCL in one package.  The result
is a bit messy if you have multiple drivers installed.  You can see which
drivers are available by starting python and running:

    import pyopencl as cl
    cl.create_some_context(interactive=True)

Once you have done that, it will show the available drivers which you
can select.  It will then tell you that you can use these drivers
automatically by setting the PYOPENCL_CTX environment variable.

Some graphics cards have multiple devices on the same card.  You cannot
yet use both of them concurrently to evaluate models, but you can run
the program twice using a different device for each session.

OpenCL kernels are compiled when needed by the device driver.  Some
drivers produce compiler output even when there is no error.  You
can see the output by setting PYOPENCL_COMPILER_OUTPUT=1.  It should be
harmless, albeit annoying.
"""
import os
import warnings

import numpy as np

try:
    import pyopencl as cl
    # Ask OpenCL for the default context so that we know that one exists
    cl.create_some_context(interactive=False)
except Exception as exc:
    warnings.warn(str(exc))
    raise RuntimeError("OpenCL not available")

from pyopencl import mem_flags as mf
from pyopencl.characterize import get_fast_inaccurate_build_options

from . import generate

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

def has_type(device, dtype):
    """
    Return true if device supports the requested precision.
    """
    if dtype == generate.F32:
        return True
    elif dtype == generate.F64:
        return "cl_khr_fp64" in device.extensions
    elif dtype == generate.F16:
        return "cl_khr_fp16" in device.extensions
    else:
        return False

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


def compile_model(context, source, dtype, fast=False):
    """
    Build a model to run on the gpu.

    Returns the compiled program and its type.  The returned type will
    be float32 even if the desired type is float64 if any of the
    devices in the context do not support the cl_khr_fp64 extension.
    """
    dtype = np.dtype(dtype)
    if not all(has_type(d, dtype) for d in context.devices):
        raise RuntimeError("%s not supported for devices"%dtype)

    source = generate.convert_type(source, dtype)
    # Note: USE_SINCOS makes the intel cpu slower under opencl
    if context.devices[0].type == cl.device_type.GPU:
        source = "#define USE_SINCOS\n" + source
    options = (get_fast_inaccurate_build_options(context.devices[0])
               if fast else [])
    program = cl.Program(context, source).build(options=options)
    return program


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
        self.queues = [cl.CommandQueue(context, context.devices[0])
                       for context in self.context]
        self.compiled = {}

    def has_type(self, dtype):
        """
        Return True if all devices support a given type.
        """
        dtype = generate.F32 if dtype == 'fast' else np.dtype(dtype)
        return any(has_type(d, dtype)
                   for context in self.context
                   for d in context.devices)

    def get_queue(self, dtype):
        """
        Return a command queue for the kernels of type dtype.
        """
        for context, queue in zip(self.context, self.queues):
            if all(has_type(d, dtype) for d in context.devices):
                return queue

    def get_context(self, dtype):
        """
        Return a OpenCL context for the kernels of type dtype.
        """
        for context, queue in zip(self.context, self.queues):
            if all(has_type(d, dtype) for d in context.devices):
                return context

    def _create_some_context(self):
        """
        Protected call to cl.create_some_context without interactivity.  Use
        this if PYOPENCL_CTX is set in the environment.  Sets the *context*
        attribute.
        """
        try:
            self.context = [cl.create_some_context(interactive=False)]
        except Exception as exc:
            warnings.warn(str(exc))
            warnings.warn("pyopencl.create_some_context() failed")
            warnings.warn("the environment variable 'PYOPENCL_CTX' might not be set correctly")

    def compile_program(self, name, source, dtype, fast=False):
        """
        Compile the program for the device in the given context.
        """
        key = "%s-%s-%s"%(name, dtype, fast)
        if key not in self.compiled:
            #print("compiling",name)
            dtype = np.dtype(dtype)
            program = compile_model(self.get_context(dtype), source, dtype, fast)
            self.compiled[key] = program
        return self.compiled[key]

    def release_program(self, name):
        """
        Free memory associated with the program on the device.
        """
        if name in self.compiled:
            self.compiled[name].release()
            del self.compiled[name]

def _get_default_context():
    """
    Get an OpenCL context, preferring GPU over CPU, and preferring Intel
    drivers over AMD drivers.
    """
    # Note: on mobile devices there is automatic clock scaling if either the
    # CPU or the GPU is underutilized; probably doesn't affect us, but we if
    # it did, it would mean that putting a busy loop on the CPU while the GPU
    # is running may increase throughput.
    #
    # Macbook pro, base install:
    #     {'Apple': [Intel CPU, NVIDIA GPU]}
    # Macbook pro, base install:
    #     {'Apple': [Intel CPU, Intel GPU]}
    # 2 x nvidia 295 with Intel and NVIDIA opencl drivers installed
    #     {'Intel': [CPU], 'NVIDIA': [GPU, GPU, GPU, GPU]}
    gpu, cpu = None, None
    for platform in cl.get_platforms():
        # AMD provides a much weaker CPU driver than Intel/Apple, so avoid it.
        # If someone has bothered to install the AMD/NVIDIA drivers, prefer them over the integrated
        # graphics driver that may have been supplied with the CPU chipset.
        preferred_cpu = platform.vendor.startswith('Intel') or platform.vendor.startswith('Apple')
        preferred_gpu = platform.vendor.startswith('Advanced') or platform.vendor.startswith('NVIDIA')
        for device in platform.get_devices():
            if device.type == cl.device_type.GPU:
                # If the existing type is not GPU then it will be CUSTOM or ACCELERATOR,
                # so don't override it.
                if gpu is None or (preferred_gpu and gpu.type == cl.device_type.GPU):
                    gpu = device
            elif device.type == cl.device_type.CPU:
                if cpu is None or preferred_cpu:
                    cpu = device
            else:
                # System has cl.device_type.ACCELERATOR or cl.device_type.CUSTOM
                # Intel Phi for example registers as an accelerator
                # Since the user installed a custom device on their system and went through the
                # pain of sorting out OpenCL drivers for it, lets assume they really do want to
                # use it as their primary compute device.
                gpu = device

    # order the devices by gpu then by cpu; when searching for an available device by data type they
    # will be checked in this order, which means that if the gpu supports double then the cpu will never
    # be used (though we may make it possible to explicitly request the cpu at some point).
    devices = []
    if gpu is not None:
        devices.append(gpu)
    if cpu is not None:
        devices.append(cpu)
    return [cl.Context([d]) for d in devices]


class GpuModel(object):
    """
    GPU wrapper for a single model.

    *source* and *model_info* are the model source and interface as returned
    from :func:`generate.make_source` and :func:`generate.make_model_info`.

    *dtype* is the desired model precision.  Any numpy dtype for single
    or double precision floats will do, such as 'f', 'float32' or 'single'
    for single and 'd', 'float64' or 'double' for double.  Double precision
    is an optional extension which may not be available on all devices.
    Half precision ('float16','half') may be available on some devices.
    Fast precision ('fast') is a loose version of single precision, indicating
    that the compiler is allowed to take shortcuts.
    """
    def __init__(self, source, model_info, dtype=generate.F32):
        self.info = model_info
        self.source = source
        self.dtype = generate.F32 if dtype == 'fast' else np.dtype(dtype)
        self.fast = (dtype == 'fast')
        self.program = None # delay program creation

    def __getstate__(self):
        return self.info, self.source, self.dtype, self.fast

    def __setstate__(self, state):
        self.info, self.source, self.dtype, self.fast = state
        self.program = None

    def __call__(self, q_vectors):
        if self.program is None:
            compiler = environment().compile_program
            self.program = compiler(self.info['name'], self.source, self.dtype,
                                    self.fast)
        is_2d = len(q_vectors) == 2
        kernel_name = generate.kernel_name(self.info, is_2d)
        kernel = getattr(self.program, kernel_name)
        return GpuKernel(kernel, self.info, q_vectors, self.dtype)

    def release(self):
        """
        Free the resources associated with the model.
        """
        if self.program is not None:
            environment().release_program(self.info['name'])
            self.program = None

    def __del__(self):
        self.release()

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
        # TODO: do we ever need double precision q?
        env = environment()
        self.nq = q_vectors[0].size
        self.dtype = np.dtype(dtype)
        self.is_2d = (len(q_vectors) == 2)
        # TODO: stretch input based on get_warp()
        # not doing it now since warp depends on kernel, which is not known
        # at this point, so instead using 32, which is good on the set of
        # architectures tested so far.
        self.q_vectors = [_stretch_input(q, self.dtype, 32) for q in q_vectors]
        context = env.get_context(self.dtype)
        self.global_size = [self.q_vectors[0].size]
        #print("creating inputs of size", self.global_size)
        self.q_buffers = [
            cl.Buffer(context, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=q)
            for q in self.q_vectors
        ]

    def release(self):
        """
        Free the memory.
        """
        for b in self.q_buffers:
            b.release()
        self.q_buffers = []

    def __del__(self):
        self.release()

class GpuKernel(object):
    """
    Callable SAS kernel.

    *kernel* is the GpuKernel object to call

    *model_info* is the module information

    *q_vectors* is the q vectors at which the kernel should be evaluated

    *dtype* is the kernel precision

    The resulting call method takes the *pars*, a list of values for
    the fixed parameters to the kernel, and *pd_pars*, a list of (value,weight)
    vectors for the polydisperse parameters.  *cutoff* determines the
    integration limits: any points with combined weight less than *cutoff*
    will not be calculated.

    Call :meth:`release` when done with the kernel instance.
    """
    def __init__(self, kernel, model_info, q_vectors, dtype):
        q_input = GpuInput(q_vectors, dtype)
        self.kernel = kernel
        self.info = model_info
        self.res = np.empty(q_input.nq, q_input.dtype)
        dim = '2d' if q_input.is_2d else '1d'
        self.fixed_pars = model_info['partype']['fixed-' + dim]
        self.pd_pars = model_info['partype']['pd-' + dim]

        # Inputs and outputs for each kernel call
        # Note: res may be shorter than res_b if global_size != nq
        env = environment()
        self.queue = env.get_queue(dtype)
        self.loops_b = cl.Buffer(self.queue.context, mf.READ_WRITE,
                                 2 * MAX_LOOPS * q_input.dtype.itemsize)
        self.res_b = cl.Buffer(self.queue.context, mf.READ_WRITE,
                               q_input.global_size[0] * q_input.dtype.itemsize)
        self.q_input = q_input

        self._need_release = [self.loops_b, self.res_b, self.q_input]

    def __call__(self, fixed_pars, pd_pars, cutoff):
        real = (np.float32 if self.q_input.dtype == generate.F32
                else np.float64 if self.q_input.dtype == generate.F64
                else np.float16 if self.q_input.dtype == generate.F16
                else np.float32)  # will never get here, so use np.float32

        #print "pars", fixed_pars, pd_pars
        res_bi = self.res_b
        nq = np.uint32(self.q_input.nq)
        if pd_pars:
            cutoff = real(cutoff)
            loops_N = [np.uint32(len(p[0])) for p in pd_pars]
            loops = np.hstack(pd_pars) \
                if pd_pars else np.empty(0, dtype=self.q_input.dtype)
            loops = np.ascontiguousarray(loops.T, self.q_input.dtype).flatten()
            #print("loops",Nloops, loops)

            #import sys; print("opencl eval",pars)
            #print("opencl eval",pars)
            if len(loops) > 2 * MAX_LOOPS:
                raise ValueError("too many polydispersity points")

            loops_bi = self.loops_b
            cl.enqueue_copy(self.queue, loops_bi, loops)
            loops_l = cl.LocalMemory(len(loops.data))
            #ctx = environment().context
            #loops_bi = cl.Buffer(ctx, mf.READ_ONLY|mf.COPY_HOST_PTR, hostbuf=loops)
            dispersed = [loops_bi, loops_l, cutoff] + loops_N
        else:
            dispersed = []
        fixed = [real(p) for p in fixed_pars]
        args = self.q_input.q_buffers + [res_bi, nq] + dispersed + fixed
        self.kernel(self.queue, self.q_input.global_size, None, *args)
        cl.enqueue_copy(self.queue, self.res, res_bi)

        return self.res

    def release(self):
        """
        Release resources associated with the kernel.
        """
        for v in self._need_release:
            v.release()
        self._need_release = []

    def __del__(self):
        self.release()
