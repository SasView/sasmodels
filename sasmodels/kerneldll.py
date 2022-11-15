r"""
DLL driver for C kernels

If the environment variable *SAS_OPENMP* is set, then sasmodels
will attempt to compile with OpenMP flags so that the model can use all
available kernels.  This may or may not be available on your compiler
toolchain.  Depending on operating system and environment.

Windows does not have provide a compiler with the operating system.
Instead, we assume that TinyCC is installed and available.  This can
be done with a simple pip command if it is not already available::

    pip install tinycc

If Microsoft Visual C++ is available (because VCINSTALLDIR is
defined in the environment), then that will be used instead.
Microsoft Visual C++ for Python is available from Microsoft:

    `<http://www.microsoft.com/en-us/download/details.aspx?id=44266>`_

If neither compiler is available, sasmodels will check for *MinGW*,
the GNU compiler toolchain. This available in packages such as Anaconda
and PythonXY, or available stand alone. This toolchain has had
difficulties on some systems, and may or may not work for you.

You can control which compiler to use by setting SAS_COMPILER in the
environment:

  - tinycc (Windows): use the TinyCC compiler shipped with SasView
  - msvc (Windows): use the Microsoft Visual C++ compiler
  - mingw (Windows): use the MinGW GNU cc compiler
  - unix (Linux): use the system cc compiler.
  - unix (Mac): use the clang compiler. You will need XCode installed, and
    the XCode command line tools. Mac comes with OpenCL drivers, so generally
    this will not be needed.

Both *msvc* and *mingw* require that the compiler is available on your path.
For *msvc*, this can done by running vcvarsall.bat in a windows terminal.
Install locations are system dependent, such as:

    C:\Program Files (x86)\Common Files\Microsoft\Visual
    C++ for Python\9.0\vcvarsall.bat

or maybe

    C:\Users\yourname\AppData\Local\Programs\Common\Microsoft\Visual
    C++ for Python\9.0\vcvarsall.bat

OpenMP for *msvc* requires the Microsoft vcomp90.dll library, which doesn't
seem to be included with the compiler, nor does there appear to be a public
download location.  There may be one on your machine already in a location
such as:

    C:\Windows\winsxs\x86_microsoft.vc90.openmp*\vcomp90.dll

If you copy this to somewhere on your path, such as the python directory or
the install directory for this application, then OpenMP should be supported.

For full control of the compiler, define a function
*compile_command(source,output)* which takes the name of the source file
and the name of the output file and returns a compile command that can be
evaluated in the shell.  For even more control, replace the entire
*compile(source,output)* function.

The global attribute *ALLOW_SINGLE_PRECISION_DLLS* should be set to *False* if
you wish to prevent single precision floating point evaluation for the compiled
models, otherwise set it defaults to *True*.
"""
from __future__ import print_function

import sys
import os
from os.path import join as joinpath, splitext
import subprocess
import shlex
import tempfile
import ctypes as ct  # type: ignore
import _ctypes as _ct
import logging

import numpy as np  # type: ignore

try:
    import tinycc
except ImportError:
    tinycc = None

from . import generate
from .kernel import KernelModel, Kernel
from .kernelpy import PyInput
from .exception import annotate_exception
from .generate import F16, F32, F64

# pylint: disable=unused-import
try:
    from typing import Tuple, Callable, Any, List
    from .modelinfo import ModelInfo
    from .details import CallDetails
except ImportError:
    pass
# pylint: enable=unused-import

# Compiler output is a byte stream that needs to be decode in python 3.
decode = (lambda s: s) if sys.version_info[0] < 3 else (lambda s: s.decode('utf8'))

if "SAS_DLL_PATH" in os.environ:
    SAS_DLL_PATH = os.environ["SAS_DLL_PATH"]
else:
    # Assume the default location of module DLLs is in .sasmodels/compiled_models.
    SAS_DLL_PATH = os.path.join(os.path.expanduser("~"), ".sasmodels", "compiled_models")

if "SAS_COMPILER" in os.environ:
    COMPILER = os.environ["SAS_COMPILER"]
elif os.name == 'nt':
    if tinycc is not None:
        COMPILER = "tinycc"
    elif "VCINSTALLDIR" in os.environ:
        # If vcvarsall.bat has been called, then VCINSTALLDIR is in the
        # environment and we can use the MSVC compiler.  Otherwise, if
        # tinycc is available then use it.  Otherwise, hope that mingw
        # is available.
        COMPILER = "msvc"
    else:
        COMPILER = "mingw"
else:
    COMPILER = "unix"

ARCH = "" if ct.sizeof(ct.c_void_p) > 4 else "x86"  # 4 byte pointers on x86.
if COMPILER == "unix":
    # Generic unix compile.
    # On Mac users will need the X code command line tools installed.
    #COMPILE = "gcc-mp-4.7 -shared -fPIC -std=c99 -fopenmp -O2 -Wall %s -o %s -lm -lgomp"
    CC = os.environ.get("CC", "cc")
    CPPFLAGS = os.environ.get("CPPFLAGS", "")
    CFLAGS = os.environ.get("CFLAGS", "-std=c99 -O2 -Wall")
    LDFLAGS = os.environ.get("LDFLAGS", "")
    SOFLAGS = "-fPIC -shared"
    compiler_vars = (CC, CPPFLAGS, CFLAGS, LDFLAGS, SOFLAGS)
    compiler = [val for var in compiler_vars for val in shlex.split(var)]
    LIBS = ["-lm"] + shlex.split(os.environ.get("LIBS", ""))

    # Add OpenMP support if not running on a Mac.
    if sys.platform != "darwin":
        # OpenMP seems to be broken on gcc 5.4.0 (ubuntu 16.04.9).
        # Shut it off for all unix until we can investigate.
        #CC.append("-fopenmp")
        pass
    def compile_command(source, output):
        """unix compiler command"""
        return compiler + [source, "-o", output] + LIBS
elif COMPILER == "msvc":
    # Call vcvarsall.bat before compiling to set path, headers, libs, etc.
    # MSVC compiler is available, so use it.  OpenMP requires a copy of
    # vcomp90.dll on the path.  One may be found here:
    #       C:/Windows/winsxs/x86_microsoft.vc90.openmp*/vcomp90.dll
    # Copy this to the python directory and uncomment the OpenMP COMPILE.
    # TODO: Remove intermediate OBJ file created in the directory.
    # TODO: Maybe don't use randomized name for the c file.
    # TODO: Maybe ask distutils to find MSVC.
    CC = "cl /nologo /Ox /MD /W3 /GS- /DNDEBUG".split()
    if "SAS_OPENMP" in os.environ:
        CC.append("/openmp")
    LN = "/link /DLL /INCREMENTAL:NO /MANIFEST".split()
    def compile_command(source, output):
        """MSVC compiler command"""
        return CC + ["/Tp%s"%source] + LN + ["/OUT:%s"%output]
elif COMPILER == "tinycc":
    # TinyCC compiler.
    CC = [tinycc.TCC] + "-shared -rdynamic -Wall".split()
    def compile_command(source, output):
        """tinycc compiler command"""
        return CC + [source, "-o", output]
elif COMPILER == "mingw":
    # MinGW compiler.
    CC = "gcc -shared -std=c99 -O2 -Wall".split()
    if "SAS_OPENMP" in os.environ:
        CC.append("-fopenmp")
    def compile_command(source, output):
        """mingw compiler command"""
        return CC + [source, "-o", output, "-lm"]

ALLOW_SINGLE_PRECISION_DLLS = True


def compile_model(source, output):
    # type: (str, str) -> None
    """
    Compile *source* producing *output*.

    Raises RuntimeError if the compile failed or the output wasn't produced.
    """
    command = compile_command(source=source, output=output)
    command_str = " ".join('"%s"'%p if ' ' in p else p for p in command)
    logging.info(command_str)
    try:
        # Need shell=True on windows to keep console box from popping up.
        shell = (os.name == 'nt')
        subprocess.check_output(command, shell=shell, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as exc:
        output = decode(exc.output)
        raise RuntimeError("compile failed.\n%s\n%s"%(command_str, output))
    if not os.path.exists(output):
        raise RuntimeError("compile failed.  File is in %r"%source)


def dll_name(model_file, dtype):
    # type: (str, np.dtype) ->  str
    """
    Name of the dll containing the model.  This is the base file name without
    any path or extension, with a form such as 'sas_sphere32'.
    """
    bits = 8*dtype.itemsize
    basename = "sas%d_%s"%(bits, model_file)
    basename += ARCH + ".so"

    # Hack to find precompiled dlls.
    path = joinpath(generate.DATA_PATH, '..', 'compiled_models', basename)
    if os.path.exists(path):
        return path

    return joinpath(SAS_DLL_PATH, basename)


def dll_path(model_file, dtype):
    # type: (str, np.dtype) -> str
    """
    Complete path to the dll for the model.  Note that the dll may not
    exist yet if it hasn't been compiled.
    """
    return os.path.join(SAS_DLL_PATH, dll_name(model_file, dtype))


def make_dll(source, model_info, dtype=F64, system=False):
    # type: (str, ModelInfo, np.dtype, bool) -> str
    """
    Returns the path to the compiled model defined by *kernel_module*.

    If the model has not been compiled, or if the source file(s) are newer
    than the dll, then *make_dll* will compile the model before returning.
    This routine does not load the resulting dll.

    *dtype* is a numpy floating point precision specifier indicating whether
    the model should be single, double or long double precision.  The default
    is double precision, *np.dtype('d')*.

    Set *sasmodels.ALLOW_SINGLE_PRECISION_DLLS* to False if single precision
    models are not allowed as DLLs.

    Set *sasmodels.kerneldll.SAS_DLL_PATH* to the compiled dll output path.
    Alternatively, set the environment variable *SAS_DLL_PATH*.
    The default is in ~/.sasmodels/compiled_models.

    *system* is a bool that controls whether these are the precompiled DLLs
    that would be shipped with a binary distribution.
    """
    if dtype == F16:
        raise ValueError("16 bit floats not supported")
    if dtype == F32 and not ALLOW_SINGLE_PRECISION_DLLS:
        dtype = F64  # Force 64-bit dll.
    # Note: dtype may be F128 for long double precision.

    # TODO: Deal with ever-growing ~/.sasmodels/compiled_models.
    # TODO: Is the GPU model cache growing without bound?
    # Tag model name with hash so any change to the model or the wrapper
    # will generate a new dll. Need to tag with dtype as well because
    # double has not yet been converted to the target type in the source.
    # Don't use time stamps for caching since they are not reliable, especially
    # when multiple versions of the application are installed.
    model_file = model_info.id + "_" + generate.tag_source(source)
    dll = dll_path(model_file, dtype)
    logging.debug("make_dll: dll located %s as %s in %s",
                  model_info.id, model_file, dll)

    if not os.path.exists(dll):
        # Make sure the DLL path exists. Use abspath since python docs warn
        # that makedirs is not robust against '..' in path.
        os.makedirs(os.path.abspath(SAS_DLL_PATH), exist_ok=True)
        source = generate.convert_type(source, dtype)
        if not system:
            basename = splitext(os.path.basename(dll))[0] + "_"
            system_fd, filename = tempfile.mkstemp(suffix=".c", prefix=basename)
            logging.debug("make_dll: writing C for dll: %s", filename)
            with os.fdopen(system_fd, "w") as file_handle:
                file_handle.write(source)
        else:
            filename = splitext(dll)[0] + ".c"
            logging.debug("make_dll: writing C for system dll: %s", filename)
            with open(filename, 'w') as file_handle:
                file_handle.write(source)
        compile_model(source=filename, output=dll)
        # Comment the following to keep the generated C file.
        # Note: If there is a syntax error then compile raises an error
        # and the source file will not be deleted.
        os.unlink(filename)
        #print("saving compiled file in %r"%filename)
    else:
        logging.debug("make_dll: cache hit for %s", dll)
    return dll


def load_dll(source, model_info, dtype=F64):
    # type: (str, ModelInfo, np.dtype) -> "DllModel"
    """
    Create and load a dll corresponding to the source.

    *model_info* is the info object returned from
    :func:`.modelinfo.make_model_info`.

    *source* is returned from :func:`.generate.make_source`, as
    *make_source(model_info)['dll']*.

    See :func:`make_dll` for details on controlling the dll path and the
    allowed floating point precision.
    """
    filename = make_dll(source, model_info, dtype=dtype)
    return DllModel(filename, model_info, dtype=dtype)


class DllModel(KernelModel):
    """
    ctypes wrapper for a single model.

    *dllpath* is the stored path to the dll.

    *model_info* is the model definition returned from
    :func:`.modelinfo.make_model_info`.

    *dtype* is the desired model precision.  Any numpy dtype for single
    or double precision floats will do, such as 'f', 'float32' or 'single'
    for single and 'd', 'float64' or 'double' for double.  Double precision
    is an optional extension which may not be available on all devices.

    Call :meth:`release` when done with the kernel.
    """
    def __init__(self, dllpath, model_info, dtype=generate.F32):
        # type: (str, ModelInfo, np.dtype) -> None
        self.info = model_info
        self.dllpath = dllpath
        self._dll = None  # type: ct.CDLL
        self._kernels = None  # type: List[Callable, Callable]
        self.dtype = np.dtype(dtype)

    def _load_dll(self):
        # type: () -> None
        try:
            self._dll = ct.CDLL(self.dllpath)
        except:
            annotate_exception("while loading "+self.dllpath)
            raise

        float_type = (ct.c_float if self.dtype == generate.F32
                      else ct.c_double if self.dtype == generate.F64
                      else ct.c_longdouble)

        # int, int, int, int*, double*, double*, double*, double*, double
        argtypes = [ct.c_int32]*3 + [ct.c_void_p]*4 + [float_type, ct.c_int32]
        names = [generate.kernel_name(self.info, variant)
                 for variant in ("Iq", "Iqxy", "Imagnetic")]
        self._kernels = [self._dll[name] for name in names]
        for k in self._kernels:
            k.argtypes = argtypes

    def __getstate__(self):
        # type: () -> Tuple[ModelInfo, str]
        return self.info, self.dllpath

    def __setstate__(self, state):
        # type: (Tuple[ModelInfo, str]) -> None
        self.info, self.dllpath = state
        self._dll = None

    def make_kernel(self, q_vectors):
        # type: (List[np.ndarray]) -> DllKernel
        q_input = PyInput(q_vectors, self.dtype)
        # Note: DLL is lazy loaded.
        if self._dll is None:
            self._load_dll()
        is_2d = len(q_vectors) == 2
        kernel = self._kernels[1:3] if is_2d else [self._kernels[0]]*2
        return DllKernel(kernel, self.info, q_input)

    def release(self):
        # type: () -> None
        """
        Release any resources associated with the model.
        """
        dll_handle = self._dll._handle
        if os.name == 'nt':
            ct.windll.kernel32.FreeLibrary(dll_handle)
        else:
            _ct.dlclose(dll_handle)
        del self._dll
        self._dll = None


class DllKernel(Kernel):
    """
    Callable SAS kernel.

    *kernel* is the c function to call.

    *model_info* is the module information

    *q_input* is the DllInput q vectors at which the kernel should be
    evaluated.

    The resulting call method takes the *pars*, a list of values for
    the fixed parameters to the kernel, and *pd_pars*, a list of (value, weight)
    vectors for the polydisperse parameters.  *cutoff* determines the
    integration limits: any points with combined weight less than *cutoff*
    will not be calculated.

    Call :meth:`release` when done with the kernel instance.
    """
    def __init__(self, kernel, model_info, q_input):
        # type: (Callable[[], np.ndarray], ModelInfo, PyInput) -> None
        dtype = q_input.dtype
        self.q_input = q_input
        self.kernel = kernel

        # Attributes accessed from the outside.
        self.dim = '2d' if q_input.is_2d else '1d'
        self.info = model_info
        self.dtype = dtype

        # Converter to translate input to target type.
        self._as_dtype = (np.float32 if dtype == generate.F32
                          else np.float64 if dtype == generate.F64
                          else np.float128)

        # Holding place for the returned value.
        nout = 2 if self.info.have_Fq else 1
        extra_q = 4  # Total weight, form volume, shell volume and R_eff.
        self.result = np.empty(self.q_input.nq*nout + extra_q, dtype)

    def _call_kernel(self, call_details, values, cutoff, magnetic,
                     radius_effective_mode):
        # type: (CallDetails, np.ndarray, float, bool, int)

        # Setup kernel function and arguments.
        kernel = self.kernel[1 if magnetic else 0]
        kernel_args = [
            self.q_input.nq,  # Number of inputs.
            None,  # Placeholder for pd_start.
            None,  # Placeholder for pd_stop.
            call_details.buffer.ctypes.data,  # Problem definition.
            values.ctypes.data,  # Parameter values.
            self.q_input.q.ctypes.data,  # Q values.
            self.result.ctypes.data,   # Result storage.
            self._as_dtype(cutoff),  # Probability cutoff.
            radius_effective_mode,  # R_eff mode.
        ]

        # Call kernel and retrieve results.
        #print("Calling DLL")
        #call_details.show(values)
        step = 100
        # TODO: Do we need the explicit sleep like the OpenCL and CUDA loops?
        for start in range(0, call_details.num_eval, step):
            stop = min(start + step, call_details.num_eval)
            kernel_args[1:3] = [start, stop]
            kernel(*kernel_args) # type: ignore

    def release(self):
        # type: () -> None
        """
        Release resources associated with the kernel.
        """
        # TODO: OpenCL/CUDA allocate q_input in __init__ and free it in release.
        # Should we be doing the same for DLL?
        #self.q_input.release()
        pass

    def __del__(self):
        # type: () -> None
        self.release()
