.. currentmodule:: sasmodels
.. Wim Bouwman, DUT, written at codecamp-V, Oct2016

.. _SESANS:

SANS to SESANS conversion
=========================

The conversion from SANS into SESANS in absolute units is a simple Hankel transformation when all the small-angle scattered neutrons are detected.
First we calculate the Hankel transform including the absolute intensities by

.. math:: G(\delta) = 2 \pi \int_0^{\infty} J_0(Q \delta) \frac{d \Sigma}{d \Omega} (Q) Q dQ \!,

in which :math:`J_0` is the zeroth order Bessel function, :math:`\delta` the spin-echo length, :math:`Q` the wave vector transfer
and :math:`\frac{d \Sigma}{d \Omega} (Q)` the scattering cross section in absolute units.
This is a 1-dimensional integral, which can be rather fast.
In the numerical calculation we integrate from :math:`Q_{min} = 0.1 \times 2 \pi / R_{max}` in which :math:`R_{max}` will be model dependent.
We determined the factor 0.1 by varying its value until the value of the integral was stable.
This happened at a value of 0.3.  The have a safety margin of a factor of three we have choosen the value 0.1.
For the solid sphere we took 3 times the radius for :math:`R_{max}`.
The real integration is performed to :math:`Q_{max}` which is an instrumental parameter that is read in from the measurement file.
From the equation above we can calculate the polarisation that we measure in a SESANS experiment:

.. math:: P(\delta) = e^{t \left( \frac{ \lambda}{2 \pi} \right)^2 \left(G(\delta) - G(0) \right)} \!,

in which :math:`t` is the thickness of the sample and :math:`\lambda` is the wavelength of the neutrons.