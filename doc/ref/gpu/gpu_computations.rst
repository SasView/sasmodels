.. _models-computation:

****************
GPU Computations
****************
SasView model evaluations can run on your graphics card (GPU) or they can run
on the processor (CPU). In general, calculations performed on the GPU will run faster.

To run on the GPU, your computer must have OpenCL drivers installed.
For information about OpenCL installation see this
:ref:`opencl-installation` guidance.

Where the model is evaluated is a little bit complicated.
If the model has the line *single=False* then it requires double precision.
If the GPU is single precision only, then it will try running via OpenCL
on the CPU.  If the OpenCL driver is not available for the CPU then
it will run as a normal program on the CPU.

For models with a large number of parameters or with a lot of code,
the GPU may be too small to run the program effectively.
In this case, you should try simplifying the model, maybe breaking it
into several different modules so that you don't need *IF* statements in your code.
If it is still too big, you can set *opencl=False* in the model file and
the model will only run as a normal program on the CPU.
This will not usually be necessary.

Device Selection
................
If you have multiple GPU devices you can tell SasView which device to use.
By default, SasView looks for one GPU and one CPU device
from available OpenCL platforms.

SasView prefers AMD or NVIDIA drivers for GPU, and prefers Intel or
Apple drivers for CPU. Both GPU and CPU are included on the assumption that CPU 
is always available and supports double precision.

The device order is important: GPU is checked before CPU on the assumption that
it will be faster. By examining ~/sasview.log you can see which device SasView
chose to run the model.

**If you don't want to use OpenCL, you can set** *SAS_OPENCL=None*
**in your environment settings, and it will only use normal programs.**

If you want to use one of the other devices, you can run the following
from the python console in SasView::

    import pyopencl as cl
    cl.create_some_context()

This will provide a menu of different OpenCL drivers available.
When one is selected, it will say "set PYOPENCL_CTX=..."
Use that value as the value of *SAS_OPENCL*.

Compiler Selection
..................
For models run as normal programs, you may need to specify a compiler.
This is done using the SAS_COMPILER environment variable.
Set it to *tinycc* for the tinycc compiler, *msvc* for the
Microsoft Visual C compiler, or *mingw* for the MinGW compiler.
TinyCC is provided with SasView so that is the default.
If you want one of the other compilers, be sure to have it available
in your *PATH* so SasView can find it!


.. note::
    This help document was last changed by Steve King, 08Oct2016
