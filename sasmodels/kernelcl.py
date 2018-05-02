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
automatically by setting the SAS_OPENCL environment variable, which is
PYOPENCL_CTX equivalent but not conflicting with other pyopnecl programs.

Some graphics cards have multiple devices on the same card.  You cannot
yet use both of them concurrently to evaluate models, but you can run
the program twice using a different device for each session.

OpenCL kernels are compiled when needed by the device driver.  Some
drivers produce compiler output even when there is no error.  You
can see the output by setting PYOPENCL_COMPILER_OUTPUT=1.  It should be
harmless, albeit annoying.
"""
from __future__ import print_function

import os
import warnings
import logging
import time

import numpy as np  # type: ignore


# Attempt to setup opencl. This may fail if the opencl package is not
# installed or if it is installed but there are no devices available.
try:
    import pyopencl as cl  # type: ignore
    from pyopencl import mem_flags as mf
    from pyopencl.characterize import get_fast_inaccurate_build_options
    # Ask OpenCL for the default context so that we know that one exists
    cl.create_some_context(interactive=False)
    HAVE_OPENCL = True
    OPENCL_ERROR = ""
except Exception as exc:
    HAVE_OPENCL = False
    OPENCL_ERROR = str(exc)

from . import generate
from .kernel import KernelModel, Kernel

# pylint: disable=unused-import
try:
    from typing import Tuple, Callable, Any
    from .modelinfo import ModelInfo
    from .details import CallDetails
except ImportError:
    pass
# pylint: enable=unused-import

# CRUFT: pyopencl < 2017.1  (as of June 2016 needs quotes around include path)
def quote_path(v):
    """
    Quote the path if it is not already quoted.

    If v starts with '-', then assume that it is a -I option or similar
    and do not quote it.  This is fragile:  -Ipath with space needs to
    be quoted.
    """
    return '"'+v+'"' if v and ' ' in v and not v[0] in "\"'-" else v

def fix_pyopencl_include():
    """
    Monkey patch pyopencl to allow spaces in include file path.
    """
    import pyopencl as cl
    if hasattr(cl, '_DEFAULT_INCLUDE_OPTIONS'):
        cl._DEFAULT_INCLUDE_OPTIONS = [quote_path(v) for v in cl._DEFAULT_INCLUDE_OPTIONS]

if HAVE_OPENCL:
    fix_pyopencl_include()

# The max loops number is limited by the amount of local memory available
# on the device.  You don't want to make this value too big because it will
# waste resources, nor too small because it may interfere with users trying
# to do their polydispersity calculations.  A value of 1024 should be much
# larger than necessary given that cost grows as npts^k where k is the number
# of polydisperse parameters.
MAX_LOOPS = 2048


# Pragmas for enable OpenCL features.  Be sure to protect them so that they
# still compile even if OpenCL is not present.
_F16_PRAGMA = """\
#if defined(__OPENCL_VERSION__) // && !defined(cl_khr_fp16)
#  pragma OPENCL EXTENSION cl_khr_fp16: enable
#endif
"""

_F64_PRAGMA = """\
#if defined(__OPENCL_VERSION__) // && !defined(cl_khr_fp64)
#  pragma OPENCL EXTENSION cl_khr_fp64: enable
#endif
"""

def use_opencl():
    return HAVE_OPENCL and os.environ.get("SAS_OPENCL", "").lower() != "none"

ENV = None
def reset_environment():
    """
    Call to create a new OpenCL context, such as after a change to SAS_OPENCL.
    """
    global ENV
    ENV = GpuEnvironment() if use_opencl() else None

def environment():
    # type: () -> "GpuEnvironment"
    """
    Returns a singleton :class:`GpuEnvironment`.

    This provides an OpenCL context and one queue per device.
    """
    if ENV is None:
        if not HAVE_OPENCL:
            raise RuntimeError("OpenCL startup failed with ***"
                               + OPENCL_ERROR + "***; using C compiler instead")
        reset_environment()
        if ENV is None:
            raise RuntimeError("SAS_OPENCL=None in environment")
    return ENV

def has_type(device, dtype):
    # type: (cl.Device, np.dtype) -> bool
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
    # type: (cl.Kernel, cl.CommandQueue) -> int
    """
    Return the size of an execution batch for *kernel* running on *queue*.
    """
    return kernel.get_work_group_info(
        cl.kernel_work_group_info.PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
        queue.device)

def _stretch_input(vector, dtype, extra=1e-3, boundary=32):
    # type: (np.ndarray, np.dtype, float, int) -> np.ndarray
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
    # type: (cl.Context, str, np.dtype, bool) -> cl.Program
    """
    Build a model to run on the gpu.

    Returns the compiled program and its type.

    Raises an error if the desired precision is not available.
    """
    dtype = np.dtype(dtype)
    if not all(has_type(d, dtype) for d in context.devices):
        raise RuntimeError("%s not supported for devices"%dtype)

    source_list = [generate.convert_type(source, dtype)]

    if dtype == generate.F16:
        source_list.insert(0, _F16_PRAGMA)
    elif dtype == generate.F64:
        source_list.insert(0, _F64_PRAGMA)

    # Note: USE_SINCOS makes the intel cpu slower under opencl
    if context.devices[0].type == cl.device_type.GPU:
        source_list.insert(0, "#define USE_SINCOS\n")
    options = (get_fast_inaccurate_build_options(context.devices[0])
               if fast else [])
    source = "\n".join(source_list)
    program = cl.Program(context, source).build(options=options)
    #print("done with "+program)
    return program


# for now, this returns one device in the context
# TODO: create a context that contains all devices on all platforms
class GpuEnvironment(object):
    """
    GPU context, with possibly many devices, and one queue per device.
    """
    def __init__(self):
        # type: () -> None
        # find gpu context
        #self.context = cl.create_some_context()

        self.context = None
        if 'SAS_OPENCL' in os.environ:
            #Setting PYOPENCL_CTX as a SAS_OPENCL to create cl context
            os.environ["PYOPENCL_CTX"] = os.environ["SAS_OPENCL"]
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
        # type: (np.dtype) -> bool
        """
        Return True if all devices support a given type.
        """
        return any(has_type(d, dtype)
                   for context in self.context
                   for d in context.devices)

    def get_queue(self, dtype):
        # type: (np.dtype) -> cl.CommandQueue
        """
        Return a command queue for the kernels of type dtype.
        """
        for context, queue in zip(self.context, self.queues):
            if all(has_type(d, dtype) for d in context.devices):
                return queue

    def get_context(self, dtype):
        # type: (np.dtype) -> cl.Context
        """
        Return a OpenCL context for the kernels of type dtype.
        """
        for context in self.context:
            if all(has_type(d, dtype) for d in context.devices):
                return context

    def _create_some_context(self):
        # type: () -> cl.Context
        """
        Protected call to cl.create_some_context without interactivity.  Use
        this if SAS_OPENCL is set in the environment.  Sets the *context*
        attribute.
        """
        try:
            self.context = [cl.create_some_context(interactive=False)]
        except Exception as exc:
            warnings.warn(str(exc))
            warnings.warn("pyopencl.create_some_context() failed")
            warnings.warn("the environment variable 'SAS_OPENCL' might not be set correctly")

    def compile_program(self, name, source, dtype, fast, timestamp):
        # type: (str, str, np.dtype, bool, float) -> cl.Program
        """
        Compile the program for the device in the given context.
        """
        # Note: PyOpenCL caches based on md5 hash of source, options and device
        # so we don't really need to cache things for ourselves.  I'll do so
        # anyway just to save some data munging time.
        tag = generate.tag_source(source)
        key = "%s-%s-%s%s"%(name, dtype, tag, ("-fast" if fast else ""))
        # Check timestamp on program
        program, program_timestamp = self.compiled.get(key, (None, np.inf))
        if program_timestamp < timestamp:
            del self.compiled[key]
        if key not in self.compiled:
            context = self.get_context(dtype)
            logging.info("building %s for OpenCL %s", key,
                         context.devices[0].name.strip())
            program = compile_model(self.get_context(dtype),
                                    str(source), dtype, fast)
            self.compiled[key] = (program, timestamp)
        return program

def _get_default_context():
    # type: () -> List[cl.Context]
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
        # If someone has bothered to install the AMD/NVIDIA drivers, prefer
        # them over the integrated graphics driver that may have been supplied
        # with the CPU chipset.
        preferred_cpu = (platform.vendor.startswith('Intel')
                         or platform.vendor.startswith('Apple'))
        preferred_gpu = (platform.vendor.startswith('Advanced')
                         or platform.vendor.startswith('NVIDIA'))
        for device in platform.get_devices():
            if device.type == cl.device_type.GPU:
                # If the existing type is not GPU then it will be CUSTOM
                # or ACCELERATOR so don't override it.
                if gpu is None or (preferred_gpu and gpu.type == cl.device_type.GPU):
                    gpu = device
            elif device.type == cl.device_type.CPU:
                if cpu is None or preferred_cpu:
                    cpu = device
            else:
                # System has cl.device_type.ACCELERATOR or cl.device_type.CUSTOM
                # Intel Phi for example registers as an accelerator
                # Since the user installed a custom device on their system
                # and went through the pain of sorting out OpenCL drivers for
                # it, lets assume they really do want to use it as their
                # primary compute device.
                gpu = device

    # order the devices by gpu then by cpu; when searching for an available
    # device by data type they will be checked in this order, which means
    # that if the gpu supports double then the cpu will never be used (though
    # we may make it possible to explicitly request the cpu at some point).
    devices = []
    if gpu is not None:
        devices.append(gpu)
    if cpu is not None:
        devices.append(cpu)
    return [cl.Context([d]) for d in devices]


class GpuModel(KernelModel):
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
    def __init__(self, source, model_info, dtype=generate.F32, fast=False):
        # type: (Dict[str,str], ModelInfo, np.dtype, bool) -> None
        self.info = model_info
        self.source = source
        self.dtype = dtype
        self.fast = fast
        self.program = None # delay program creation
        self._kernels = None

    def __getstate__(self):
        # type: () -> Tuple[ModelInfo, str, np.dtype, bool]
        return self.info, self.source, self.dtype, self.fast

    def __setstate__(self, state):
        # type: (Tuple[ModelInfo, str, np.dtype, bool]) -> None
        self.info, self.source, self.dtype, self.fast = state
        self.program = None

    def make_kernel(self, q_vectors):
        # type: (List[np.ndarray]) -> "GpuKernel"
        if self.program is None:
            compile_program = environment().compile_program
            timestamp = generate.ocl_timestamp(self.info)
            self.program = compile_program(
                self.info.name,
                self.source['opencl'],
                self.dtype,
                self.fast,
                timestamp)
            variants = ['Iq', 'Iqxy', 'Imagnetic']
            names = [generate.kernel_name(self.info, k) for k in variants]
            kernels = [getattr(self.program, k) for k in names]
            self._kernels = dict((k, v) for k, v in zip(variants, kernels))
        is_2d = len(q_vectors) == 2
        if is_2d:
            kernel = [self._kernels['Iqxy'], self._kernels['Imagnetic']]
        else:
            kernel = [self._kernels['Iq']]*2
        return GpuKernel(kernel, self.dtype, self.info, q_vectors)

    def release(self):
        # type: () -> None
        """
        Free the resources associated with the model.
        """
        if self.program is not None:
            self.program = None

    def __del__(self):
        # type: () -> None
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
        # type: (List[np.ndarray], np.dtype) -> None
        # TODO: do we ever need double precision q?
        env = environment()
        self.nq = q_vectors[0].size
        self.dtype = np.dtype(dtype)
        self.is_2d = (len(q_vectors) == 2)
        # TODO: stretch input based on get_warp()
        # not doing it now since warp depends on kernel, which is not known
        # at this point, so instead using 32, which is good on the set of
        # architectures tested so far.
        if self.is_2d:
            # Note: 16 rather than 15 because result is 1 longer than input.
            width = ((self.nq+16)//16)*16
            self.q = np.empty((width, 2), dtype=dtype)
            self.q[:self.nq, 0] = q_vectors[0]
            self.q[:self.nq, 1] = q_vectors[1]
        else:
            # Note: 32 rather than 31 because result is 1 longer than input.
            width = ((self.nq+32)//32)*32
            self.q = np.empty(width, dtype=dtype)
            self.q[:self.nq] = q_vectors[0]
        self.global_size = [self.q.shape[0]]
        context = env.get_context(self.dtype)
        #print("creating inputs of size", self.global_size)
        self.q_b = cl.Buffer(context, mf.READ_ONLY | mf.COPY_HOST_PTR,
                             hostbuf=self.q)

    def release(self):
        # type: () -> None
        """
        Free the memory.
        """
        if self.q_b is not None:
            self.q_b.release()
            self.q_b = None

    def __del__(self):
        # type: () -> None
        self.release()

class GpuKernel(Kernel):
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
    def __init__(self, kernel, dtype, model_info, q_vectors):
        # type: (cl.Kernel, np.dtype, ModelInfo, List[np.ndarray]) -> None
        q_input = GpuInput(q_vectors, dtype)
        self.kernel = kernel
        self.info = model_info
        self.dtype = dtype
        self.dim = '2d' if q_input.is_2d else '1d'
        # plus three for the normalization values
        self.result = np.empty(q_input.nq+1, dtype)

        # Inputs and outputs for each kernel call
        # Note: res may be shorter than res_b if global_size != nq
        env = environment()
        self.queue = env.get_queue(dtype)

        self.result_b = cl.Buffer(self.queue.context, mf.READ_WRITE,
                                  q_input.global_size[0] * dtype.itemsize)
        self.q_input = q_input # allocated by GpuInput above

        self._need_release = [self.result_b, self.q_input]
        self.real = (np.float32 if dtype == generate.F32
                     else np.float64 if dtype == generate.F64
                     else np.float16 if dtype == generate.F16
                     else np.float32)  # will never get here, so use np.float32

    def __call__(self, call_details, values, cutoff, magnetic):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool) -> np.ndarray
        context = self.queue.context
        # Arrange data transfer to card
        details_b = cl.Buffer(context, mf.READ_ONLY | mf.COPY_HOST_PTR,
                              hostbuf=call_details.buffer)
        values_b = cl.Buffer(context, mf.READ_ONLY | mf.COPY_HOST_PTR,
                             hostbuf=values)

        kernel = self.kernel[1 if magnetic else 0]
        args = [
            np.uint32(self.q_input.nq), None, None,
            details_b, values_b, self.q_input.q_b, self.result_b,
            self.real(cutoff),
        ]
        #print("Calling OpenCL")
        #call_details.show(values)
        # Call kernel and retrieve results
        wait_for = None
        last_nap = time.clock()
        step = 1000000//self.q_input.nq + 1
        for start in range(0, call_details.num_eval, step):
            stop = min(start + step, call_details.num_eval)
            #print("queuing",start,stop)
            args[1:3] = [np.int32(start), np.int32(stop)]
            wait_for = [kernel(self.queue, self.q_input.global_size, None,
                               *args, wait_for=wait_for)]
            if stop < call_details.num_eval:
                # Allow other processes to run
                wait_for[0].wait()
                current_time = time.clock()
                if current_time - last_nap > 0.5:
                    time.sleep(0.05)
                    last_nap = current_time
        cl.enqueue_copy(self.queue, self.result, self.result_b)
        #print("result", self.result)

        # Free buffers
        for v in (details_b, values_b):
            if v is not None:
                v.release()

        pd_norm = self.result[self.q_input.nq]
        scale = values[0]/(pd_norm if pd_norm != 0.0 else 1.0)
        background = values[1]
        #print("scale",scale,values[0],self.result[self.q_input.nq],background)
        return scale*self.result[:self.q_input.nq] + background
        # return self.result[:self.q_input.nq]

    def release(self):
        # type: () -> None
        """
        Release resources associated with the kernel.
        """
        for v in self._need_release:
            v.release()
        self._need_release = []

    def __del__(self):
        # type: () -> None
        self.release()
