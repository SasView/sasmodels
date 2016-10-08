.. _opencl-installation:

*******************
OpenCL Installation
*******************
*Warning! GPU devices do not in general offer the same level of memory protection as CPU devices. If your code attempts to write outside allocated memory buffers unpredicatable behaviour may result (eg, your video display may freeze, or your system may crash, etc). Do not install OpenCL drivers without first checking for known issues (eg, some computer manufacturers install modified graphics drivers so replacing these may not be a good idea!). If in doubt, seek advice from an IT professional before proceeding further.*

1. Check if you have OpenCL already installed
=============================================

Windows
.......
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

Mac OSX
.......
    For OS X operating systems higher than 10.6 OpenCL is shipped along with the system.

    However, OpenCL has had a rocky history on Macs. Apple provide a useful compatibility table at https://support.apple.com/en-us/HT202823


2. Installation
===============

Windows
.......
    Depending on the graphic card in your system, drivers
    can be obtained from different sources:

    * NVIDIA: https://developer.nvidia.com/opencl
    * AMD: http://developer.amd.com/tools-and-sdks/opencl-zone/

Mac OSX
.......
    N/A
	
	You cannot download OpenCL driver updates for your Mac. They are packaged with the normal quarterly OS X updates from Apple. 


.. note::
    Intel provides OpenCL drivers for Intel processors at     https://software.intel.com/en-us/articles/opencl-drivers 
    These can sometimes make use of special vector instructions across multiple
    processors, so it is worth installing if the GPU does not support double
    precision. You can install this driver alongside the GPU driver for NVIDIA
    or AMD.

	
.. note::
    This help document was last changed by Steve King, 08Oct2016
