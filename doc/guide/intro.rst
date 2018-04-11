
.. _models-intro:

************
Introduction
************

This software provides form factors for various particle shapes and computes
resolution, polydispersity and angular dispersion. After giving a mathematical
definition of each model, we show the list of parameters available to the user.
Validation plots for each model are also presented.

To easily compare to the scattering intensity measured in experiments, we
normalize the form factors by the volume of the particle

.. math::

    P(\vec q) = \frac{P_o(\vec q)}{V} = \frac{1}{V} F(\vec q) F^*(\vec q)

with

.. math::

    F(\vec q) = \int\int\int dV\rho(\vec r) e^{-i\vec q \cdot \vec r}

where $P_0(\vec q)$ is the un-normalized form factor, $\rho(\vec r)$ is
the scattering length density at a given point in space and the integration
is done over the volume $V$ of the scatterer.

For systems without inter-particle interference, the form factors we provide
can be related to the scattering intensity by the particle volume fraction

.. math::

    I(\vec q) = \Phi P(\vec q)

Our so-called 1D scattering intensity functions provide $P(Q)$ for the case
where the scatterer is randomly oriented. In that case, the scattering
intensity only depends on the length of $Q$ . The intensity measured on
the plane of the SAS detector will have an azimuthal symmetry around $Q=0$.

Our so-called 2D scattering intensity functions provide $P(Q,\phi)$ for an
oriented system as a function of a $q$ vector in the plane of the detector.
We define the angle $\phi$ as the angle between the $q$ vector and the
horizontal ($x$) axis of the plane of the detector.

For information about polarised and magnetic scattering, see :ref:`magnetism`.

The models are used within the SasView package.  Instructions on how to
use SasView itself are available separately.

Many of our models use the form factor calculations as implemented in the
C-library provided by the NIST Center for Neutron Research and thus some
content and figures in this document are originated from or shared with
the NIST SANS Igor-based analysis package.

*Document History*

| 2017-05-07 Paul Kienzle