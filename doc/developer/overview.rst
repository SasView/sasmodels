.. py:currentmodule:: sasmodels

***************************
Code Overview
***************************

Computational kernels
---------------------

At the heart of *sasmodels* are the individual computational kernels.  These
functions take a particular $q$ value and a set of parameter values and
return the expected scattering for that $q$. The instructions for writing
a kernel are documented in :ref:`Writing_a_Plugin`.  The source code for
the kernels is stored in :mod:`sasmodels.models`.

The primary interface to the models is through :mod:`sasmodels.core`, which
provides functions for listing available models, loading the model definition
and compiling the model.  Use :func:`sasmodels.core.load_model` to load in
a model definition and compile it.

The computational kernels are difficult to call directly, but an example
was provided in :ref:`Scripting_Interface`.

The code for loading the kernel is in
:mod:`sasmodels.generate`, but a big part of it, which is interpreting the
definition, lives in :mod:`sasmodels.modelinfo`.  This defines the properties
of the model in :class:`sasmodels.modelinfo.ModelInfo`, the available
model parameters in :class:`sasmodels.modelinfo.ParameterTable` and the
individual parameter attributes such as units and hard limits in
:class:`sasmodels.modelinfo.Parameter`.

The kernel functions for the most part do not define polydispersity,
resolution functions or magnetism.  Instead, code for these is generated in
:mod:`sasmodels.generate`



Sasview kernels:

* :mod:`
    ('kernel', 'Evaluator type definitions'),
    ('kernelcl', 'OpenCL model evaluator'),
    ('kerneldll', 'Ctypes model evaluator'),
    ('kernelpy', 'Python model evaluator'),

Bumps fitting routines:

* :mod:`sasmodels.bumps_model` is a wrapper around
    ('bumps_model', 'Bumps interface'),
    ('direct_model', 'Simple interface'),


Utility functions:

* :mod:`sasmodels` provides :func:`sasmodels.__init__.data_files` for the installer.
* :mod:`sasmodels.rst2html` provides tools for converting model docs to html and viewing
  the html.  When run as a main program, it can display arbitrary rst files.
* :mod:`sasmodels.exception` annotates the current exception with a context string,
  such as "while opening myfile.dat" without adjusting the traceback.


And so on::


    ('compare_many', 'Batch compare models on different compute engines'),
    ('compare', 'Compare models on different compute engines'),
    ('convert', 'Sasview to sasmodel converter'),
    ('core', 'Model access'),
    ('data', 'Data layout and plotting routines'),
    ('details', 'Parameter packing for kernel calls'),
    ('generate', 'Model parser'),
    ('list_pars', 'Identify all parameters in all models'),
    ('mixture', 'Mixture model evaluator'),
    ('model_test', 'Unit test support'),
    ('modelinfo', 'Parameter and model definitions'),
    ('product', 'Product model evaluator'),
    ('resolution', '1-D resolution functions'),
    ('resolution2d', '2-D resolution functions'),
    ('sasview_model', 'Sasview interface'),
    ('sesans', 'SESANS calculation routines'),
    ('weights', 'Distribution functions'),
