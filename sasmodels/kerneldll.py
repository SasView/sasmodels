r"""
DLL driver for C kernels

The global attribute *ALLOW_SINGLE_PRECISION_DLLS* should be set to *True* if
you wish to allow single precision floating point evaluation for the compiled
models, otherwise it defaults to *False*.

The compiler command line is stored in the attribute *COMPILE*, with string
substitutions for %(source)s and %(output)s indicating what to compile and
where to store it.  The actual command is system dependent.

On windows systems, you have a choice of compilers.  *MinGW* is the GNU
compiler toolchain, available in packages such as anaconda and PythonXY,
or available stand alone. This toolchain has had difficulties on some
systems, and may or may not work for you.  In order to build DLLs, *gcc*
must be on your path.  If the environment variable *SAS_OPENMP* is given
then -fopenmp is added to the compiler flags.  This requires a version
of MinGW compiled with OpenMP support.

An alternative toolchain uses the Microsoft Visual C++ compiler, available
free from microsoft:

    `<http://www.microsoft.com/en-us/download/details.aspx?id=44266>`_

Again, this requires that the compiler is available on your path.  This is
done by running vcvarsall.bat in a windows terminal.  Install locations are
system dependent, such as:

    C:\Program Files (x86)\Common Files\Microsoft\Visual C++ for Python\9.0\vcvarsall.bat

or maybe

    C:\Users\yourname\AppData\Local\Programs\Common\Microsoft\Visual C++ for Python\9.0\vcvarsall.bat

And again, the environment variable *SAS_OPENMP* controls whether OpenMP is
used to compile the C code.  This requires the Microsoft vcomp90.dll library,
which doesn't seem to be included with the compiler, nor does there appear
to be a public download location.  There may be one on your machine already
in a location such as:

    C:\Windows\winsxs\x86_microsoft.vc90.openmp*\vcomp90.dll

If you copy this onto your path, such as the python directory or the install
directory for this application, then OpenMP should be supported.
"""
from __future__ import print_function

import sys
import os
from os.path import join as joinpath, split as splitpath, splitext
import subprocess
import tempfile
import ctypes as ct  # type: ignore
from ctypes import c_void_p, c_int32, c_longdouble, c_double, c_float  # type: ignore
import logging

import numpy as np  # type: ignore

from . import generate
from .kernel import KernelModel, Kernel
from .kernelpy import PyInput
from .exception import annotate_exception
from .generate import F16, F32, F64

try:
    from typing import Tuple, Callable, Any
    from .modelinfo import ModelInfo
    from .details import CallDetails
except ImportError:
    pass

if os.name == 'nt':
    ARCH = "" if sys.maxint > 2**32 else "x86"  # maxint=2**31-1 on 32 bit
    # Windows compiler; check if TinyCC is available
    try:
        import tinycc
    except ImportError:
        tinycc = None
    # call vcvarsall.bat before compiling to set path, headers, libs, etc.
    if "VCINSTALLDIR" in os.environ:
        # MSVC compiler is available, so use it.  OpenMP requires a copy of
        # vcomp90.dll on the path.  One may be found here:
        #       C:/Windows/winsxs/x86_microsoft.vc90.openmp*/vcomp90.dll
        # Copy this to the python directory and uncomment the OpenMP COMPILE
        # TODO: remove intermediate OBJ file created in the directory
        # TODO: maybe don't use randomized name for the c file
        # TODO: maybe ask distutils to find MSVC
        CC = "cl /nologo /Ox /MD /W3 /GS- /DNDEBUG".split()
        if "SAS_OPENMP" in os.environ:
            CC.append("/openmp")
        LN = "/link /DLL /INCREMENTAL:NO /MANIFEST".split()
        def compile_command(source, output):
            return CC + ["/Tp%s"%source] + LN + ["/OUT:%s"%output]
    elif tinycc:
        # TinyCC compiler.
        CC = [tinycc.TCC] + "-shared -rdynamic -Wall".split()
        def compile_command(source, output):
            return CC + [source, "-o", output]
    else:
        # MinGW compiler.
        CC = "gcc -shared -std=c99 -O2 -Wall".split()
        if "SAS_OPENMP" in os.environ:
            CC.append("-fopenmp")
        def compile_command(source, output):
            return CC + [source, "-o", output, "-lm"]
else:
    ARCH = ""
    # Generic unix compile
    # On mac users will need the X code command line tools installed
    #COMPILE = "gcc-mp-4.7 -shared -fPIC -std=c99 -fopenmp -O2 -Wall %s -o %s -lm -lgomp"
    CC = "cc -shared -fPIC -std=c99 -O2 -Wall".split()
    # add openmp support if not running on a mac
    if sys.platform != "darwin":
        CC.append("-fopenmp")
    def compile_command(source, output):
        return CC + [source, "-o", output, "-lm"]

# Windows-specific solution
if os.name == 'nt':
    # Assume the default location of module DLLs is in .sasmodels/compiled_models.
    DLL_PATH = os.path.join(os.path.expanduser("~"), ".sasmodels", "compiled_models")
    if not os.path.exists(DLL_PATH):
        os.makedirs(DLL_PATH)
else:
    # Set up the default path for compiled modules.
    DLL_PATH = tempfile.gettempdir()

ALLOW_SINGLE_PRECISION_DLLS = True

def compile(source, output):
    command = compile_command(source=source, output=output)
    command_str = " ".join('"%s"'%p if ' ' in p else p for p in command)
    logging.info(command_str)
    try:
        # need shell=True on windows to keep console box from popping up
        shell = (os.name == 'nt')
        subprocess.check_output(command, shell=shell, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError("compile failed.\n%s\n%s"%(command_str, exc.output))
    if not os.path.exists(output):
        raise RuntimeError("compile failed.  File is in %r"%source)

def dll_name(model_info, dtype):
    # type: (ModelInfo, np.dtype) ->  str
    """
    Name of the dll containing the model.  This is the base file name without
    any path or extension, with a form such as 'sas_sphere32'.
    """
    bits = 8*dtype.itemsize
    basename = "sas%d_%s"%(bits, model_info.id)
    basename += ARCH + ".so"

    # Hack to find precompiled dlls
    path = joinpath(generate.DATA_PATH, '..', 'compiled_models', basename)
    if os.path.exists(path):
        return path

    return joinpath(DLL_PATH, basename)


def dll_path(model_info, dtype):
    # type: (ModelInfo, np.dtype) -> str
    """
    Complete path to the dll for the model.  Note that the dll may not
    exist yet if it hasn't been compiled.
    """
    return os.path.join(DLL_PATH, dll_name(model_info, dtype))


def make_dll(source, model_info, dtype=F64):
    # type: (str, ModelInfo, np.dtype) -> str
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

    Set *sasmodels.kerneldll.DLL_PATH* to the compiled dll output path.
    The default is the system temporary directory.
    """
    if dtype == F16:
        raise ValueError("16 bit floats not supported")
    if dtype == F32 and not ALLOW_SINGLE_PRECISION_DLLS:
        dtype = F64  # Force 64-bit dll
    # Note: dtype may be F128 for long double precision

    dll = dll_path(model_info, dtype)

    if not os.path.exists(dll):
        need_recompile = True
    elif getattr(sys, 'frozen', None) is not None:
        # TODO: don't suppress time stamp
        # Currently suppressing recompile when running in a frozen environment
        need_recompile = False
    else:
        dll_time = os.path.getmtime(dll)
        newest_source = generate.timestamp(model_info)
        need_recompile = dll_time < newest_source
    if need_recompile:
        basename = os.path.splitext(os.path.basename(dll))[0] + "_"
        fd, filename = tempfile.mkstemp(suffix=".c", prefix=basename)
        source = generate.convert_type(source, dtype)
        with os.fdopen(fd, "w") as file:
            file.write(source)
        compile(source=filename, output=dll)
        # comment the following to keep the generated c file
        # Note: if there is a syntax error then compile raises an error
        # and the source file will not be deleted.
        os.unlink(filename)
        #print("saving compiled file in %r"%filename)
    return dll


def load_dll(source, model_info, dtype=F64):
    # type: (str, ModelInfo, np.dtype) -> "DllModel"
    """
    Create and load a dll corresponding to the source, info pair returned
    from :func:`sasmodels.generate.make` compiled for the target precision.

    See :func:`make_dll` for details on controlling the dll path and the
    allowed floating point precision.
    """
    filename = make_dll(source, model_info, dtype=dtype)
    return DllModel(filename, model_info, dtype=dtype)


class DllModel(KernelModel):
    """
    ctypes wrapper for a single model.

    *source* and *model_info* are the model source and interface as returned
    from :func:`gen.make`.

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
        self.dtype = np.dtype(dtype)

    def _load_dll(self):
        # type: () -> None
        #print("dll", self.dllpath)
        try:
            self._dll = ct.CDLL(self.dllpath)
        except:
            annotate_exception("while loading "+self.dllpath)
            raise

        fp = (c_float if self.dtype == generate.F32
              else c_double if self.dtype == generate.F64
              else c_longdouble)

        # int, int, int, int*, double*, double*, double*, double*, double
        argtypes = [c_int32]*3 + [c_void_p]*4 + [fp]
        self._Iq = self._dll[generate.kernel_name(self.info, "Iq")]
        self._Iqxy = self._dll[generate.kernel_name(self.info, "Iqxy")]
        self._Imagnetic = self._dll[generate.kernel_name(self.info, "Imagnetic")]
        self._Iq.argtypes = argtypes
        self._Iqxy.argtypes = argtypes
        self._Imagnetic.argtypes = argtypes

    def __getstate__(self):
        # type: () -> Tuple[ModelInfo, str]
        return self.info, self.dllpath

    def __setstate__(self, state):
        # type: (Tuple[ModelInfo, str]) -> None
        self.info, self.dllpath = state
        self._dll = None

    def make_kernel(self, q_vectors, magnetic=False):
        # type: (List[np.ndarray]) -> DllKernel
        q_input = PyInput(q_vectors, self.dtype)
        # Note: pickle not supported for DllKernel
        if self._dll is None:
            self._load_dll()
        kernel = [self._Iqxy, self._Imagnetic] if q_input.is_2d else self._Iq
        return DllKernel(kernel, self.info, q_input)

    def release(self):
        # type: () -> None
        """
        Release any resources associated with the model.
        """
        if os.name == 'nt':
            #dll = ct.cdll.LoadLibrary(self.dllpath)
            dll = ct.CDLL(self.dllpath)
            libHandle = dll._handle
            #libHandle = ct.c_void_p(dll._handle)
            del dll, self._dll
            self._dll = None
            ct.windll.kernel32.FreeLibrary(libHandle)
        else:    
            pass 


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
        self.kernel = kernel
        self.info = model_info
        self.q_input = q_input
        self.dtype = q_input.dtype
        self.dim = '2d' if q_input.is_2d else '1d'
        self.result = np.empty(q_input.nq+1, q_input.dtype)
        self.real = (np.float32 if self.q_input.dtype == generate.F32
                     else np.float64 if self.q_input.dtype == generate.F64
                     else np.float128)

    def __call__(self, call_details, values, cutoff, magnetic):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool) -> np.ndarray

        #print("in kerneldll")
        #print("values", values)
        start, stop = 0, call_details.pd_prod
        args = [
            self.q_input.nq, # nq
            start, # pd_start
            stop, # pd_stop pd_stride[MAX_PD]
            call_details.buffer.ctypes.data, # problem
            values.ctypes.data,  #pars
            self.q_input.q.ctypes.data, #q
            self.result.ctypes.data,   # results
            self.real(cutoff), # cutoff
            ]
        #print("calling DLL")
        self.kernel[1 if magnetic else 0](*args) # type: ignore
        return self.result[:-1]

    def release(self):
        # type: () -> None
        """
        Release any resources associated with the kernel.
        """
        self.q_input.release()
