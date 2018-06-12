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

Log Spaced SESANS
-----------------

For computational efficiency, the integral in the Hankel transform is
converted into a Reimann sum


.. math:: G(\delta) \approx
	  2 \pi
	  \sum_{Q=q_{min}}^{q_{max}} J_0(Q \delta)
	  \frac{d \Sigma}{d \Omega} (Q)
	  Q \Delta Q \!

However, this model approximates more than is strictly necessary.
Specifically, it is approximating the entire integral, when it is only
the scattering function that cannot be handled analytically.  A better
approximation might be

.. math:: G(\delta) \approx
	  \sum_{n=0} 2 \pi \frac{d \Sigma}{d \Omega} (q_n)
	  \int_{q_{n-1}}^{q_n} J_0(Q \delta) Q dQ
	  =
	  \sum_{n=0} \frac{2 \pi}{\delta} \frac{d \Sigma}{d \Omega} (q_n)
	  (q_n J_1(q_n \delta) - q_{n-1}J_1(q_{n-1} \delta))\!,

Assume that vectors :math:`q_n` and :math:`I_n` represent the q points
and corresponding intensity data, respectively.  Further assume that
:math:`\delta_m` and :math:`G_m` are the spin echo lengths and
corresponding Hankel transform value.

.. math:: G_m = H_{nm} I_n

where

.. math:: H_{nm} = \frac{2 \pi}{\delta_m}
	  (q_n J_1(q_n \delta_m) - q_{n-1} J_1(q_{n-1} \delta_m))

Also not that, for the limit as :math:`\delta_m` approaches zero,

.. math:: G(0)
	  =
	  \sum_{n=0} \pi \frac{d \Sigma}{d \Omega} (q_n) (q_n^2 - q_{n-1}^2)
