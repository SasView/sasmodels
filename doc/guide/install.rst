********************
sasmodels Setup
********************


sasmodels Installation
**********************
Sasmodels can be installed using a simple pip installation::

    # Won't work yet! For now follow the developer instructions below.
    $ pip install sasmodels

There are a number of dependent python packages which need to be installed
separately::

    numpy
    scipy
    opencl (for speed if your system has OpenCL drivers)
    tinycc (windows only, if no C compiler and no OpenCL)

When using sasmodels to fit data::

    sasview (for loading data)
    bumps (for fitting)
    matplotlib (for displaying results)

For development and testing it is handy to have some additional packages::

    nose (for testing)
    xmlrunner (for testing)
    matplotlib (for showing models)
    docutils (for showing model docs)
    sphinx (for building docs)
    wx ([optional] for adjusting parameters interactively)

**Developer Installation**

Developers will need to clone the sasmodels reposistory from github::

    $ git clone https://github.com/sasview/sasmodels.git

or if you have a github account::

    $ git clone git@github.com:sasview/sasmodels.git

Sasmodels can be run in-place by including the sasmodels directory on the
python path.  There are numerous ways to do this which will not be outlined
here.  It will be easiest to install the package in "develop" mode using::

    $ cd sasmodels
    $ python setup.py develop

This will allow you to edit the files in the package and have the changes
show up immediately in python the next time you load your program.

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

sasmodels evaluations can run on your graphics card (GPU) or they can run
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
================
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
==================
For models run as normal programs, you may need to specify a compiler.
This is done using the SAS_COMPILER environment variable.
Set it to *tinycc* for the tinycc compiler, *msvc* for the
Microsoft Visual C compiler, or *mingw* for the MinGW compiler.
TinyCC is provided with SasView so that is the default.
If you want one of the other compilers, be sure to have it available
in your *PATH* so SasView can find it!


*Document History*

| 2017-05-17 Paul Kienzle