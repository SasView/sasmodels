Sasmodels
=========

Theory models for small angle scattering.

The models provided are usable directly in the bumps fitting package and
in the sasview analysis package.  If OpenCL is available, the models will
run much faster.  If not, then precompiled versions will be included with
the distributed package.  New models can be added if OpenCL or a C compiler
is available.

Install
-------

The easiest way to use sasmodels is from `SasView <http://www.sasview.org/>`_.

You can also install sasmodels as a standalone package in python. Use
`miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
or `anaconda <https://www.anaconda.com/>`_
to create a python environment with the sasmodels dependencies::

    $ conda create -n sasmodels -c conda-forge numpy scipy matplotlib pyopencl

The option ``-n sasmodels`` names the environment sasmodels, and the option
``-c conda-forge`` selects the conda-forge package channel because pyopencl
is not part of the base anaconda distribution.

Activate the environment and install sasmodels::

    $ conda activate sasmodels
    (sasmodels) $ pip install sasmodels

Install `bumps <https://github.com/bumps/bumps>`_ if you want to use it to fit
your data::

    (sasmodels) $ pip install bumps

Usage
-----

Check that the works::

    (sasmodels) $ python -m sasmodels.compare cylinder

To show the orientation explorer::

    (sasmodels) $ python -m sasmodels.jitter

Documentation is available online as part of the SasView
`fitting perspective <https://www.sasview.org/docs/index.html>`_
as well as separate pages for
`individual models <https://www.sasview.org/docs/user/qtgui/Perspectives/Fitting/models/index.html>`_.
Programming details for sasmodels are available in the
`developer documentation <https://www.sasview.org/docs/dev/dev.html>`_.


Fitting Example
---------------

The example directory contains a radial+tangential data set for an oriented
rod-like shape.

To load the example data, you will need the SAS data loader from the sasview
package. This is not yet available on PyPI, so you will need a copy of the
SasView source code to run it.  Create a directory somewhere to hold the
sasview and sasmodels source code, which we will refer to as $SOURCE.

Use the following to install sasview, and the sasmodels examples::

    (sasmodels) $ cd $SOURCE
    (sasmodels) $ conda install git
    (sasmodels) $ git clone https://github.com/sasview/sasview.git
    (sasmodels) $ git clone https://github.com/sasview/sasmodels.git

Set the path to the sasview source on your python path within the sasmodels
environment.  On Windows, this will be::

    (sasmodels)> set PYTHONPATH="$SOURCE\sasview\src"
    (sasmodels)> cd $SOURCE/sasmodels/example
    (sasmodels)> python -m bumps.cli fit.py cylinder --preview

On Mac/Linux with the standard shell this will be::

    (sasmodels) $ export PYTHONPATH="$SOURCE/sasview/src"
    (sasmodels) $ cd $SOURCE/sasmodels/example
    (sasmodels) $ bumps fit.py cylinder --preview

The fit.py model accepts up to two arguments.  The first argument is the
model type, which has been defined for cylinder, capped_cylinder,
core_shell_cylinder, ellipsoid, triaxial_ellipsoid and lamellar.  The
second argument is view, which can be radial or tangential.  To fit
both radial and tangential simultaneously, use the word "both".

See `bumps documentation <https://bumps.readthedocs.io/>`_ for detailed
instructions on running the fit.

|TestStatus|_ |TravisStatus|_

.. |TestStatus| image:: https://github.com/SasView/sasmodels/workflows/Test/badge.svg?branch=master
.. _TestStatus: https://github.com/SasView/sasmodels/actions?query=workflow%3ATest+branch%3Amaster

.. |TravisStatus| image:: https://travis-ci.org/SasView/sasmodels.svg?branch=master
.. _TravisStatus: https://travis-ci.org/SasView/sasmodels
