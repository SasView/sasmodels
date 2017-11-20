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
