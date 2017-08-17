.. currentmodule:: sasmodels
.. Wim Bouwman, DUT, written at codecamp-V, Oct2016

.. _SESANS:

SANS to SESANS conversion
=========================

The conversion from SANS into SESANS in absolute units is a simple Hankel
transformation when all the small-angle scattered neutrons are detected.
First we calculate the Hankel transform including the absolute intensities by

.. math:: G(\delta) = 2 \pi \int_0^{\infty} J_0(Q \delta) \frac{d \Sigma}{d \Omega} (Q) Q dQ \!,

in which :math:`J_0` is the zeroth order Bessel function, :math:`\delta`
the spin-echo length, :math:`Q` the wave vector transfer and :math:`\frac{d \Sigma}{d \Omega} (Q)`
the scattering cross section in absolute units.

Out of necessity, a 1-dimensional numerical integral approximates the exact Hankel transform.
The upper bound of the numerical integral is :math:`Q_{max}`, which is calculated from the wavelength and the instrument's maximum acceptance angle, both of which are included in the file.
While the true Hankel transform has a lower bound of zero, most scattering models are undefined at :math: `Q=0`, so the integral requires an effective lower bound.
The lower bound of the integral is :math:`Q_{min} = 0.1 \times 2 \pi / R_{max}`, in which :math:`R_{max}` is the maximum length scale probed by the instrument multiplied by the number of data points.
This lower bound is the minimum expected Q value for the given length scale multiplied by a constant.
The constant, 0.1, was chosen empirically by integrating multiple curves and finding where the value at which the integral was stable.
A constant value of 0.3 gave numerical stability to the integral, so a factor of three safety margin was included to give the final value of 0.1.


From the equation above we can calculate the polarisation that we measure in a SESANS experiment:

.. math:: P(\delta) = e^{t \left( \frac{ \lambda}{2 \pi} \right)^2 \left(G(\delta) - G(0) \right)} \!,

in which :math:`t` is the thickness of the sample and :math:`\lambda` is the wavelength of the neutrons.
