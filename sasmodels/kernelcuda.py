"""
GPU driver for C kernels (with CUDA)

To select cuda, use SAS_OPENCL=cuda, or SAS_OPENCL=cuda:n for a particular
device number.  If no device number is specified, then look for CUDA_DEVICE=n
or a file ~/.cuda-device containing n for the device number.  Otherwise, try
all available device numbers.

TODO: docs are out of date

There should be a single GPU environment running on the system.  This
environment is constructed on the first call to :func:`environment`, and the
same environment is returned on each call.

After retrieving the environment, the next step is to create the kernel.
This is done with a call to :meth:`GpuEnvironment.compile_program`, which
returns the type of data used by the kernel.

Next a :class:`GpuInput` object should be created with the correct kind
of data.  This data object can be used by multiple kernels, for example,
if the target model is a weighted sum of multiple kernels.  The data
should include any extra evaluation points required to compute the proper
data smearing.  This need not match the square grid for 2D data if there
is an index saying which q points are active.

Together the GpuInput, the program, and a device form a :class:`GpuKernel`.
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
import logging
import time
import re
import atexit

import numpy as np  # type: ignore


# Attempt to setup CUDA. This may fail if the pycuda package is not
# installed or if it is installed but there are no devices available.
try:
    import pycuda.driver as cuda  # type: ignore
    from pycuda.compiler import SourceModule
    from pycuda.tools import make_default_context, clear_context_caches
    # Ask CUDA for the default context (so that we know that one exists)
    # then immediately throw it away in case the user doesn't want it.
    # Note: cribbed from pycuda.autoinit
    cuda.init()
    context = make_default_context()
    context.pop()
    clear_context_caches()
    del context
    HAVE_CUDA = True
    CUDA_ERROR = ""
except Exception as exc:
    HAVE_CUDA = False
    CUDA_ERROR = str(exc)

from . import generate
from .kernel import KernelModel, Kernel

# pylint: disable=unused-import
try:
    from typing import Tuple, Callable, Any, Dict, List
    from .modelinfo import ModelInfo
    from .details import CallDetails
except ImportError:
    pass
# pylint: enable=unused-import

# The max loops number is limited by the amount of local memory available
# on the device.  You don't want to make this value too big because it will
# waste resources, nor too small because it may interfere with users trying
# to do their polydispersity calculations.  A value of 1024 should be much
# larger than necessary given that cost grows as npts^k where k is the number
# of polydisperse parameters.
MAX_LOOPS = 2048


def use_cuda():
    # type: () -> bool
    """Returns True if CUDA is the default compute engine."""
    sas_opencl = os.environ.get("SAS_OPENCL", "CUDA").lower()
    return HAVE_CUDA and sas_opencl.startswith("cuda")


ENV = None
def reset_environment():
    # type: () -> None
    """
    Call to create a new OpenCL context, such as after a change to SAS_OPENCL.
    """
    global ENV
    # Free any previous allocated context.
    if ENV is not None and ENV.context is not None:
        ENV.release()
    ENV = GpuEnvironment() if use_cuda() else None


def environment():
    # type: () -> "GpuEnvironment"
    """
    Returns a singleton :class:`GpuEnvironment`.

    This provides an OpenCL context and one queue per device.
    """
    if ENV is None:
        if not HAVE_CUDA:
            raise RuntimeError("CUDA startup failed with ***"
                               + CUDA_ERROR + "***; using C compiler instead")
        reset_environment()
        if ENV is None:
            raise RuntimeError("SAS_OPENCL=None in environment")
    return ENV


# PyTest is not freeing ENV, so make sure it gets freed.
atexit.register(lambda: ENV.release() if ENV is not None else None)


def has_type(dtype):
    # type: (np.dtype) -> bool
    """
    Return true if device supports the requested precision.
    """
    # Assume the NVIDIA card supports 32-bit and 64-bit floats.
    # TODO: Check if pycuda support F16.
    return dtype in (generate.F32, generate.F64)


FUNCTION_PATTERN = re.compile(r"""^
  (?P<space>\s*)                       # Initial space.
  (?P<qualifiers>^(?:\s*\b\w+\b\s*)+)  # One or more qualifiers before function.
  (?P<function>\s*\b\w+\b\s*[(])       # Function name plus open parens.
  """, re.VERBOSE|re.MULTILINE)

MARKED_PATTERN = re.compile(r"""
  \b(return|else|kernel|device|__device__)\b
  """, re.VERBOSE|re.MULTILINE)


def _add_device_tag(match):
    # type: (None) -> str
    # Note: Should be re.Match, but that isn't a simple type.
    """
    replace qualifiers with __device__ qualifiers if needed
    """
    qualifiers = match.group("qualifiers")
    if MARKED_PATTERN.search(qualifiers):
        start, end = match.span()
        return match.string[start:end]
    else:
        function = match.group("function")
        space = match.group("space")
        return "".join((space, "__device__ ", qualifiers, function))


def mark_device_functions(source):
    # type: (str) -> str
    """
    Mark all function declarations as __device__ functions (except kernel).
    """
    return FUNCTION_PATTERN.sub(_add_device_tag, source)


def show_device_functions(source):
    # type: (str) -> str
    """
    Show all discovered function declarations, but don't change any.
    """
    for match in FUNCTION_PATTERN.finditer(source):
        print(match.group('qualifiers').replace('\n', r'\n'),
              match.group('function'), '(')
    return source


def compile_model(source, dtype, fast=False):
    # type: (str, np.dtype, bool) -> SourceModule
    """
    Build a model to run on the gpu.

    Returns the compiled program and its type.  The returned type will
    be float32 even if the desired type is float64 if any of the
    devices in the context do not support the cl_khr_fp64 extension.
    """
    dtype = np.dtype(dtype)
    if not has_type(dtype):
        raise RuntimeError("%s not supported for devices"%dtype)

    source_list = [generate.convert_type(source, dtype)]

    source_list.insert(0, "#define USE_SINCOS\n")
    source = "\n".join(source_list)
    #source = show_device_functions(source)
    source = mark_device_functions(source)
    #with open('/tmp/kernel.cu', 'w') as fd: fd.write(source)
    #print(source)
    #options = ['--verbose', '-E']
    options = ['--use_fast_math'] if fast else None
    program = SourceModule(source, no_extern_c=True, options=options) #, include_dirs=[...])

    #print("done with "+program)
    return program


# For now, this returns one device in the context.
# TODO: Create a context that contains all devices on all platforms.
class GpuEnvironment(object):
    """
    GPU context for CUDA.
    """
    context = None # type: cuda.Context
    def __init__(self, devnum=None):
        # type: (int) -> None
        env = os.environ.get("SAS_OPENCL", "").lower()
        if devnum is None and env.startswith("cuda:"):
            devnum = int(env[5:])

        # Set the global context to the particular device number if one is
        # given, otherwise use the default context.  Perhaps this will be set
        # by an environment variable within autoinit.
        if devnum is not None:
            self.context = cuda.Device(devnum).make_context()
        else:
            self.context = make_default_context()

        ## Byte boundary for data alignment.
        #self.data_boundary = max(d.min_data_type_align_size
        #                         for d in self.context.devices)

        # Cache for compiled programs, and for items in context.
        self.compiled = {}

    def release(self):
        """Free the CUDA device associated with this context."""
        if self.context is not None:
            self.context.pop()
            self.context = None

    def __del__(self):
        self.release()

    def has_type(self, dtype):
        # type: (np.dtype) -> bool
        """
        Return True if all devices support a given type.
        """
        return has_type(dtype)

    def compile_program(self, name, source, dtype, fast, timestamp):
        # type: (str, str, np.dtype, bool, float) -> SourceModule
        """
        Compile the program for the device in the given context.
        """
        # Note: PyCuda (probably) caches but I'll do so as well just to
        # save some data munging time.
        tag = generate.tag_source(source)
        key = "%s-%s-%s%s"%(name, dtype, tag, ("-fast" if fast else ""))
        # Check timestamp on program.
        program, program_timestamp = self.compiled.get(key, (None, np.inf))
        if program_timestamp < timestamp:
            del self.compiled[key]
        if key not in self.compiled:
            logging.info("building %s for CUDA", key)
            program = compile_model(str(source), dtype, fast)
            self.compiled[key] = (program, timestamp)
        return program


class GpuModel(KernelModel):
    """
    GPU wrapper for a single model.

    *source* and *model_info* are the model source and interface as returned
    from :func:`.generate.make_source` and :func:`.modelinfo.make_model_info`.

    *dtype* is the desired model precision.  Any numpy dtype for single
    or double precision floats will do, such as 'f', 'float32' or 'single'
    for single and 'd', 'float64' or 'double' for double.  Double precision
    is an optional extension which may not be available on all devices.
    Half precision ('float16','half') may be available on some devices.
    Fast precision ('fast') is a loose version of single precision, indicating
    that the compiler is allowed to take shortcuts.
    """
    info = None  # type: ModelInfo
    source = ""  # type: str
    dtype = None  # type: np.dtype
    fast = False  # type: bool
    _program = None  # type: SourceModule
    _kernels = None  # type: Dict[str, cuda.Function]

    def __init__(self, source, model_info, dtype=generate.F32, fast=False):
        # type: (Dict[str,str], ModelInfo, np.dtype, bool) -> None
        self.info = model_info
        self.source = source
        self.dtype = dtype
        self.fast = fast

    def __getstate__(self):
        # type: () -> Tuple[ModelInfo, str, np.dtype, bool]
        return self.info, self.source, self.dtype, self.fast

    def __setstate__(self, state):
        # type: (Tuple[ModelInfo, str, np.dtype, bool]) -> None
        self.info, self.source, self.dtype, self.fast = state
        self._program = self._kernels = None

    def make_kernel(self, q_vectors):
        # type: (List[np.ndarray]) -> "GpuKernel"
        return GpuKernel(self, q_vectors)

    def get_function(self, name):
        # type: (str) -> cuda.Function
        """
        Fetch the kernel from the environment by name, compiling it if it
        does not already exist.
        """
        if self._program is None:
            self._prepare_program()
        return self._kernels[name]

    def _prepare_program(self):
        # type: (str) -> None
        env = environment()
        timestamp = generate.ocl_timestamp(self.info)
        program = env.compile_program(
            self.info.name,
            self.source['opencl'],
            self.dtype,
            self.fast,
            timestamp)
        variants = ['Iq', 'Iqxy', 'Imagnetic']
        names = [generate.kernel_name(self.info, k) for k in variants]
        functions = [program.get_function(k) for k in names]
        self._kernels = {k: v for k, v in zip(variants, functions)}
        # Keep a handle to program so GC doesn't collect.
        self._program = program


# TODO: Check that we don't need a destructor for buffers which go out of scope.
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
    :class:`GpuModel`.  Note that not all kernels support double
    precision, so even if the program was created for double precision,
    the *GpuModel.dtype* may be single precision.

    Call :meth:`release` when complete.  Even if not called directly, the
    buffer will be released when the data object is freed.
    """
    def __init__(self, q_vectors, dtype=generate.F32):
        # type: (List[np.ndarray], np.dtype) -> None
        # TODO: Do we ever need double precision q?
        self.nq = q_vectors[0].size
        self.dtype = np.dtype(dtype)
        self.is_2d = (len(q_vectors) == 2)
        # TODO: Stretch input based on get_warp().
        # Not doing it now since warp depends on kernel, which is not known
        # at this point, so instead using 32, which is good on the set of
        # architectures tested so far.
        if self.is_2d:
            width = ((self.nq+15)//16)*16
            self.q = np.empty((width, 2), dtype=dtype)
            self.q[:self.nq, 0] = q_vectors[0]
            self.q[:self.nq, 1] = q_vectors[1]
        else:
            width = ((self.nq+31)//32)*32
            self.q = np.empty(width, dtype=dtype)
            self.q[:self.nq] = q_vectors[0]
        self.global_size = [self.q.shape[0]]
        #print("creating inputs of size", self.global_size)

        # Transfer input value to GPU.
        self.q_b = cuda.to_device(self.q)

    def release(self):
        # type: () -> None
        """
        Free the buffer associated with the q value.
        """
        if self.q_b is not None:
            self.q_b.free()
            self.q_b = None

    def __del__(self):
        # type: () -> None
        self.release()


class GpuKernel(Kernel):
    """
    Callable SAS kernel.

    *model* is the GpuModel object to call

    The kernel is derived from :class:`.kernel.Kernel`, providing the
    *_call_kernel()* method to evaluate the kernel for a given set of
    parameters. Because of the need to move the q values to the GPU before
    evaluation, the kernel is instantiated for a particular set of q vectors,
    and can be called many times without transfering q each time.

    Call :meth:`release` when done with the kernel instance.
    """
    #: SAS model information structure.
    info = None  # type: ModelInfo
    #: Kernel precision.
    dtype = None  # type: np.dtype
    #: Kernel dimensions (1d or 2d).
    dim = ""  # type: str
    #: Calculation results, updated after each call to *_call_kernel()*.
    result = None  # type: np.ndarray

    def __init__(self, model, q_vectors):
        # type: (GpuModel, List[np.ndarray]) -> None
        dtype = model.dtype
        self.q_input = GpuInput(q_vectors, dtype)
        self._model = model

        # Attributes accessed from the outside.
        self.dim = '2d' if self.q_input.is_2d else '1d'
        self.info = model.info
        self.dtype = dtype

        # Converter to translate input to target type.
        self._as_dtype = np.float64 if dtype == generate.F64 else np.float32

        # Holding place for the returned value.
        nout = 2 if self.info.have_Fq and self.dim == '1d' else 1
        extra_q = 4  # Total weight, form volume, shell volume and R_eff.
        self.result = np.empty(self.q_input.nq*nout + extra_q, dtype)

        # Allocate result value on GPU.
        width = ((self.result.size+31)//32)*32 * self.dtype.itemsize
        self._result_b = cuda.mem_alloc(width)

    def _call_kernel(self, call_details, values, cutoff, magnetic,
                     radius_effective_mode):
        # type: (CallDetails, np.ndarray, float, bool, int) -> None

        # Arrange data transfer to card.
        details_b = cuda.to_device(call_details.buffer)
        values_b = cuda.to_device(values)

        # Setup kernel function and arguments.
        name = 'Iq' if self.dim == '1d' else 'Imagnetic' if magnetic else 'Iqxy'
        kernel = self._model.get_function(name)
        kernel_args = [
            np.uint32(self.q_input.nq),  # Number of inputs.
            None,  # Placeholder for pd_start.
            None,  # Placeholder for pd_stop.
            details_b,  # Problem definition.
            values_b,  # Parameter values.
            self.q_input.q_b,  # Q values.
            self._result_b,   # Result storage.
            self._as_dtype(cutoff),  # Probability cutoff.
            np.uint32(radius_effective_mode),  # R_eff mode.
        ]
        grid = partition(self.q_input.nq)

        # Call kernel and retrieve results.
        #print("Calling CUDA")
        #call_details.show(values)
        last_nap = time.clock()
        step = 100000000//self.q_input.nq + 1
        #step = 1000000000
        for start in range(0, call_details.num_eval, step):
            stop = min(start + step, call_details.num_eval)
            #print("queuing",start,stop)
            kernel_args[1:3] = [np.int32(start), np.int32(stop)]
            kernel(*kernel_args, **grid)
            if stop < call_details.num_eval:
                sync()
                # Allow other processes to run.
                current_time = time.clock()
                if current_time - last_nap > 0.5:
                    time.sleep(0.001)
                    last_nap = current_time
        sync()
        cuda.memcpy_dtoh(self.result, self._result_b)
        #print("result", self.result)

        details_b.free()
        values_b.free()

    def release(self):
        # type: () -> None
        """
        Release resources associated with the kernel.
        """
        self.q_input.release()
        if self._result_b is not None:
            self._result_b.free()
            self._result_b = None

    def __del__(self):
        # type: () -> None
        self.release()


def sync():
    """
    Overview:
        Waits for operation in the current context to complete.

    Note: Maybe context.synchronize() is sufficient.
    """
    # Create an event with which to synchronize.
    done = cuda.Event()

    # Schedule an event trigger on the GPU.
    done.record()

    # Make sure we don't hog resource while waiting to sync.
    while not done.query():
        time.sleep(0.01)

    # Block until the GPU executes the kernel.
    done.synchronize()

    # Clean up the event; I don't think they can be reused.
    del done


def partition(n):
    """
    Constructs block and grid arguments for *n* elements.
    """
    max_gx, max_gy = 65535, 65535
    blocksize = 32
    #max_gx, max_gy = 5, 65536
    #blocksize = 3
    block = (blocksize, 1, 1)
    num_blocks = int((n+blocksize-1)/blocksize)
    if num_blocks < max_gx:
        grid = (num_blocks, 1)
    else:
        gx = max_gx
        gy = (num_blocks + max_gx - 1) / max_gx
        if gy >= max_gy:
            raise ValueError("vector is too large")
        grid = (gx, gy)
    #print("block", block, "grid", grid)
    #print("waste", block[0]*block[1]*block[2]*grid[0]*grid[1] - n)
    return dict(block=block, grid=grid)
