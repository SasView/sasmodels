.. _gpu-setup:

********************
GPU Setup
********************

SAS model evaluations can run on your graphics card (GPU) or they can run
on the processor (CPU). In general, calculations performed on the GPU
will run faster.


OpenCL Installation
*******************
*Warning! GPU devices do not in general offer the same level of memory
protection as CPU devices. If your code attempts to write outside allocated
memory buffers unpredicatable behaviour may result (eg, your video display
may freeze, or your system may crash, etc). Do not install OpenCL drivers
without first checking for known issues (eg, some computer manufacturers
install modified graphics drivers so replacing these may not be a good
idea!). If in doubt, seek advice from an IT professional before proceeding
further.*

Check if you have OpenCL already installed
==========================================

**Windows**

The following instructions are based on
http://web.engr.oregonstate.edu/~mjb/cs475/DoIHaveOpenCL.pdf

* Go to: Start -> Control Panel -> System & Security -> Administrative Tools
* Double Click on Computer Managment
* Click on Device Manager
* Click open Display Adapters
* Right-click on available adapter and select Properties
* Click on Driver
* Go to Driver Details
* Scroll down and see if OpenCL is installed (look for OpenCL*.dll files)

**Mac OSX**

For OS X operating systems higher than 10.6 OpenCL is shipped along with
the system.

However, OpenCL has had a rocky history on Macs. Apple provide a useful
compatibility table at https://support.apple.com/en-us/HT202823

Note: some versions of OS/X have an OpenCL driver for the Intel CPU that
does support precompiled OpenCL program.  The first time that PyOpenCL
compiles a program that only runs on the CPU (e.g., because the GPU does
not support double precision programs), the program will run fine, but
the next time it will give the unhelpful error:

    pyopencl._cl.RuntimeError: clBuildProgram failed:
    BUILD_PROGRAM_FAILURE - clBuildProgram failed: BUILD_PROGRAM_FAILURE -
    clBuildProgram failed: BUILD_PROGRAM_FAILURE

To prevent this you will need to set the environment variable
*PYOPENCL_NO_CACHE=1* before reruning the program.  Only do this if needed
since it will increase the time it takes to restart the program when the
cache is disabled.

Installation
============

**Windows**

Depending on the graphic card in your system, drivers
can be obtained from different sources:

* NVIDIA: https://developer.nvidia.com/opencl
* AMD: http://developer.amd.com/tools-and-sdks/opencl-zone/


**Mac OSX**

N/A

You cannot download OpenCL driver updates for your Mac. They are packaged
with the normal quarterly OS X updates from Apple.


.. note::
    Intel provides OpenCL drivers for Intel processors at
    https://software.intel.com/en-us/articles/opencl-drivers
    These can sometimes make use of special vector instructions across multiple
    processors, so it is worth installing if the GPU does not support double
    precision. You can install this driver alongside the GPU driver for NVIDIA
    or AMD.


GPU Selection
*************

The logic for choosing the compute platform is a little bit complicated.
If the model has the line *single=False* then it requires double precision.
If the GPU is single precision only, then it will try running via OpenCL
on the CPU.  If the OpenCL driver is not available for the CPU then
it will run as a normal program on the CPU.

For models with a large number of parameters or with a lot of code,
the GPU may be too small to run the program effectively. In this case, you
should try simplifying the model, maybe breaking it into several different
models so that you don't need *IF* statements in your code. If it is still
too big, you can set *opencl=False* in the model file and the model will
only run as a normal program on the CPU. This will not usually be necessary.

Device Selection
================
**OpenCL drivers**

If you have multiple GPU devices you can tell the program which device to use.
By default, the program looks for one GPU and one CPU device from available
OpenCL platforms. It prefers AMD or NVIDIA drivers for GPU, and
prefers Intel or Apple drivers for CPU. Both GPU and CPU are included on
the assumption that CPU is always available and supports double precision.

The device order is important: GPU is checked before CPU on the assumption that
it will be faster. By examining ~/sasview.log you can see which device
was used to run the model.

If you want to use a specific driver and devices, you can run the following
from the python console::

    import pyopencl as cl
    cl.create_some_context()

This will provide a menu of different OpenCL drivers available.
When one is selected, it will say "set PYOPENCL_CTX=..."
Use that value as the value of *SAS_OPENCL=driver:device*.

To use the default OpenCL device (rather than CUDA or None),
set *SAS_OPENCL=opencl*.

In batch queues, you may need to set *XDG_CACHE_HOME=~/.cache*
(Linux only) to a different directory, depending on how the filesystem
is configured.  You should also set *SAS_DLL_PATH* for CPU-only modules.


**CUDA drivers**

If OpenCL drivers are not available on your system, but NVidia CUDA
drivers are available, then set *SAS_OPENCL=cuda* or
*SAS_OPENCL=cuda:n* for a particular device number *n*.  If no device
number is specified, then the CUDA drivers looks for look for
*CUDA_DEVICE=n* or a file ~/.cuda-device containing n for the device number.

In batch queues, the SLURM command *sbatch --gres=gpu:1 ...* will set
*CUDA_VISIBLE_DEVICES=n*, which ought to set the correct device
number for *SAS_OPENCL=cuda*.  If not, then set
*CUDA_DEVICE=$CUDA_VISIBLE_DEVICES* within the batch script.  You may
need to set the CUDA cache directory to a folder accessible across the
cluster with *PYCUDA_CACHE_DIR* (or *PYCUDA_DISABLE_CACHE* to disable
caching), and you may need to set environment specific compiler flags
with *PYCUDA_DEFAULT_NVCC_FLAGS*.  You should also set *SAS_DLL_PATH*
for CPU-only modules.

**No GPU support**

If you don't want to use OpenCL or CUDA, you can set *SAS_OPENCL=None*
in your environment settings, and it will only use normal programs.

In batch queues, you may need to set *SAS_DLL_PATH* to a directory
accessible on the compute node.


Device Testing
==============
Unfortunately, not all vendors provide working OpenCL implementations
for their GPU devices.  For example, the HD 6000 Intel GPUs with
double precision support fail for some of the double precision models.

The SasView user interface provides a Fitting OpenCL Options dialog
for selecting amongst and testing the available devices.  After a
few minutes of seeming to freeze, the application will return a list
of model tests which have passed.  The same tests can be run directly
from the python console using::

    from sasmodels.model_tests import main as model_tests
    model_tests("-v", "opencl", "all")

Compiler Selection
==================
For models run as normal programs, you may need to specify a compiler.
This is done using the *SAS_COMPILER* environment variable, and the
*SAS_OPENMP* environment variable if OpenMP support is available for
the compiler.

On Windows, set *SAS_COMPILER=tinycc* for the tinycc compiler,
*SAS_COMPILER=msvc* for the Microsoft Visual C compiler,
or *SAS_COMPILER=mingw* for the MinGW compiler. If TinyCC is available
on the python path (it is provided with SasView), that will be the
default. If you want one of the other compilers, be sure to have it
available in your *PATH* so we can find it!

On Mac OS/X and Linux, set *SAS_COMPILER=unix* for the compiler.  This
will use the unix cc command to compile the model, with gcc style
command line options.  For OS/X you will need to install the Xcode
package to make the compiler available.

Environment Variables
=====================

The following environment variables are used by the system::

    SAS_MODELPATH=path - sets directory containing custom models
    SAS_OPENCL=vendor:device|cuda:device|none - sets the target GPU device
    XDG_CACHE_HOME=~/.cache - sets the pyopencl cache root (linux only)
    SAS_COMPILER=tinycc|msvc|mingw|unix - sets the DLL compiler
    SAS_OPENMP=1 - turns on OpenMP for the DLLs
    SAS_DLL_PATH=path - sets the path to the compiled modules
    SAS_NUMBA=1|2 - enables numba and numba.cuda calculations if available
    PYOPENCL_NO_CACHE=1 - turns off caching for PyOpenCL


*Document History*

| 2019-11-06 Paul Kienzle
