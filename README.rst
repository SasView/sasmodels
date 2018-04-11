Sasmodels
=========

Theory models for small angle scattering.

The models provided are usable directly in the bumps fitting package and
in the sasview analysis package.  If OpenCL is available, the models will
run much faster.  If not, then precompiled versions will be included with
the distributed package.  New models can be added if OpenCL or a C compiler
is available.

Example
-------

The example directory contains a radial+tangential data set for an oriented
rod-like shape.

The data is loaded by sas.dataloader from the sasview package, so sasview
is needed to run the example.

To run the example, you need sasview, sasmodels and bumps.  Assuming these
repositories are installed side by side, change to the sasmodels/example
directory and enter::

    PYTHONPATH=..:../../sasview/src ../../bumps/run.py fit.py \
        cylinder --preview

See bumps documentation for instructions on running the fit.  With the
python packages installed, e.g., into a virtual environment, then the
python path need not be set, and the command would be::

    bumps fit.py cylinder --preview

The fit.py model accepts up to two arguments.  The first argument is the
model type, which has been defined for cylinder, capped_cylinder,
core_shell_cylinder, ellipsoid, triaxial_ellipsoid and lamellar.  The
second argument is view, which can be radial or tangential.  To fit
both radial and tangential simultaneously, use the word "both".

Notes
-----

cylinder.c + cylinder.py is the cylinder model with renamed variables and
sld scaled by 1e6 so the numbers are nicer.  The model name is "cylinder"

lamellar.py is an example of a single file model with embedded C code.

|TravisStatus|_

.. |TravisStatus| image:: https://travis-ci.org/SasView/sasmodels.svg?branch=master
.. _TravisStatus: https://travis-ci.org/SasView/sasmodels
