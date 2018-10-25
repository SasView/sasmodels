"""
GPU driver for C kernels (with CUDA)

To select cuda, use SAS_OPENCL=cuda, or SAS_OPENCL=cuda:n for a particular
device number.  If no device number is specified, then look for CUDA_DEVICE=n
or a file ~/.cuda-device containing n for the device number.  Otherwise, try
all available device numbers.

TODO: docs are out of date

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
import re

import numpy as np  # type: ignore


# Attempt to setup cuda. This may fail if the pycuda package is not
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
    from typing import Tuple, Callable, Any
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
    env = os.environ.get("SAS_OPENCL", "").lower()
    return HAVE_CUDA and (env == "" or env.startswith("cuda"))

ENV = None
def reset_environment():
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

def has_type(dtype):
    # type: (np.dtype) -> bool
    """
    Return true if device supports the requested precision.
    """
    # Assume the nvidia card supports 32-bit and 64-bit floats.
    # TODO: check if pycuda support F16
    return dtype in (generate.F32, generate.F64)


FUNCTION_PATTERN = re.compile(r"""^
  (?P<space>\s*)                   # initial space
  (?P<qualifiers>^(?:\s*\b\w+\b\s*)+) # one or more qualifiers before function
  (?P<function>\s*\b\w+\b\s*[(])      # function name plus open parens
  """, re.VERBOSE|re.MULTILINE)

MARKED_PATTERN = re.compile(r"""
  \b(return|else|kernel|device|__device__)\b
  """, re.VERBOSE|re.MULTILINE)

def _add_device_tag(match):
    # type: (None) -> str
    # Note: should be re.Match, but that isn't a simple type
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
        print(match.group('qualifiers').replace('\n',r'\n'), match.group('function'), '(')
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
    program = SourceModule(source, no_extern_c=True, options=options) # include_dirs=[...]

    #print("done with "+program)
    return program


# for now, this returns one device in the context
# TODO: create a context that contains all devices on all platforms
class GpuEnvironment(object):
    """
    GPU context, with possibly many devices, and one queue per device.
    """
    context = None # type: cuda.Context
    def __init__(self, devnum=None):
        # type: (int) -> None
        # Byte boundary for data alignment
        #self.data_boundary = max(d.min_data_type_align_size
        #                         for d in self.context.devices)
        self.compiled = {}
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

    def release(self):
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
            logging.info("building %s for CUDA", key)
            program = compile_model(str(source), dtype, fast)
            self.compiled[key] = (program, timestamp)
        return program

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
    info = None # type: ModelInfo
    source = "" # type: str
    dtype = None # type: np.dtype
    fast = False # type: bool
    program = None # type: SourceModule
    _kernels = None # type: List[cuda.Function]

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
            kernels = [self.program.get_function(k) for k in names]
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
        #print("creating inputs of size", self.global_size)
        self.q_b = cuda.to_device(self.q)

    def release(self):
        # type: () -> None
        """
        Free the memory.
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
        self.result_b = cuda.mem_alloc(q_input.global_size[0] * dtype.itemsize)
        self.q_input = q_input # allocated by GpuInput above

        self._need_release = [self.result_b]
        self.real = (np.float32 if dtype == generate.F32
                     else np.float64 if dtype == generate.F64
                     else np.float16 if dtype == generate.F16
                     else np.float32)  # will never get here, so use np.float32

    def __call__(self, call_details, values, cutoff, magnetic):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool) -> np.ndarray
        # Arrange data transfer to card
        details_b = cuda.to_device(call_details.buffer)
        values_b = cuda.to_device(values)

        kernel = self.kernel[1 if magnetic else 0]
        args = [
            np.uint32(self.q_input.nq), None, None,
            details_b, values_b, self.q_input.q_b, self.result_b,
            self.real(cutoff),
        ]
        grid = partition(self.q_input.nq)
        #print("Calling OpenCL")
        #call_details.show(values)
        # Call kernel and retrieve results
        last_nap = time.clock()
        step = 100000000//self.q_input.nq + 1
        #step = 1000000000
        for start in range(0, call_details.num_eval, step):
            stop = min(start + step, call_details.num_eval)
            #print("queuing",start,stop)
            args[1:3] = [np.int32(start), np.int32(stop)]
            kernel(*args, **grid)
            if stop < call_details.num_eval:
                sync()
                # Allow other processes to run
                current_time = time.clock()
                if current_time - last_nap > 0.5:
                    time.sleep(0.001)
                    last_nap = current_time
        sync()
        cuda.memcpy_dtoh(self.result, self.result_b)
        #print("result", self.result)

        details_b.free()
        values_b.free()

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
        for p in self._need_release:
            p.free()
        self._need_release = []

    def __del__(self):
        # type: () -> None
        self.release()


def sync():
    """
    Overview:
        Waits for operation in the current context to complete.

    Note: Maybe context.synchronize() is sufficient.
    """
    #return # The following works in C++; don't know what pycuda is doing
    # Create an event with which to synchronize
    done = cuda.Event()

    # Schedule an event trigger on the GPU.
    done.record()

    #line added to not hog resources
    while not done.query():
        time.sleep(0.01)

    # Block until the GPU executes the kernel.
    done.synchronize()
    # Clean up the event; I don't think they can be reused.
    del done


def partition(n):
    '''
    Overview:
        Auto grids the thread blocks to achieve some level of calculation
    efficiency.
    '''
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
