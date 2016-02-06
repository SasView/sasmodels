"""
sasmodels
=========

**sasmodels** is a package containing models for small angle neutron and xray
scattering.  Models supported are the one dimensional circular average and
two dimensional oriented patterns.  As well as the form factor calculations
for the individual shapes **sasmodels** also provides automatic shape
polydispersity, angular dispersion and resolution convolution.  SESANS
patterns can be computed for any model.

Models can be written in python or in C.  C models can run on the GPU if
OpenCL drivers are available.  See :mod:`generate` for details on
defining new models.
"""

__version__ = "0.9"
