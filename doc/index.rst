sasmodels
=========
Small angle X-ray and Neutron scattering (SAXS and SANS) examines the
scattering patterns produced by a beam travelling through the sample
and scattering at low angles.  The scattering is computed as a function
of reciprocal space $q$, which arises from a combination of beam wavelength
and scattering angles. Each pixel on the detector corresponds to
a different scattering angle, and has a distinct $q_x$ and $q_y$. If the
sample is unoriented, the scattering pattern will appear as rings on the
detector.  In this case, a circular average can be taken with 1-dimension
data at $q = \surd (q_x^2 + q_y^2)$ compared to the orientationally
averaged SAS scattering pattern.

The sasmodels package provides theory functions for small angle scattering
calculations for different shapes, including the effects of resolution,
polydispersity and orientational dispersion.

.. only:: html

   :Release: |version|
   :Date:    |today|

   Download pdf version `here <SASModels.pdf>`_

.. toctree::
   :maxdepth: 4

   guide/index.rst
   guide/models/index.rst
   developer/index.rst
   api/index.rst

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`

.. only:: html

  * :ref:`search`
