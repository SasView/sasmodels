# ellipsoid model
# Note: model title and parameter table are inserted automatically
r"""
The form factor is normalized by the particle volume.

Definition
----------

The output of the 2D scattering intensity function for oriented ellipsoids
is given by (Feigin, 1987)

.. math::

    P(Q,\alpha) = {\text{scale} \over V} F^2(Q) + \text{background}

where

.. math::

    F(Q) = {3 (\Delta rho)) V (\sin[Qr(R_p,R_e,\alpha)]
                - \cos[Qr(R_p,R_e,\alpha)])
            \over [Qr(R_p,R_e,\alpha)]^3 }

and

.. math::

    r(R_p,R_e,\alpha) = \left[ R_e^2 \sin^2 \alpha
        + R_p^2 \cos^2 \alpha \right]^{1/2}


$\alpha$ is the angle between the axis of the ellipsoid and $\vec q$,
$V$ is the volume of the ellipsoid, $R_p$ is the polar radius along the
rotational axis of the ellipsoid, $R_e$ is the equatorial radius perpendicular
to the rotational axis of the ellipsoid and $\Delta \rho$ (contrast) is the
scattering length density difference between the scatterer and the solvent.

To provide easy access to the orientation of the ellipsoid, we define
the rotation axis of the ellipsoid using two angles $\theta$ and $\phi$.
These angles are defined in the
:ref:`cylinder orientation figure <cylinder-orientation>`.
For the ellipsoid, $\theta$ is the angle between the rotational axis
and the $z$-axis.

NB: The 2nd virial coefficient of the solid ellipsoid is calculated based
on the $R_p$ and $R_e$ values, and used as the effective radius for
$S(Q)$ when $P(Q) \cdot S(Q)$ is applied.

.. _ellipsoid-1d:

.. figure:: img/ellipsoid_1d.JPG

    The output of the 1D scattering intensity function for randomly oriented
    ellipsoids given by the equation above.


The $\theta$ and $\phi$ parameters are not used for the 1D output. Our
implementation of the scattering kernel and the 1D scattering intensity
use the c-library from NIST.

.. _ellipsoid-geometry:

.. figure:: img/ellipsoid_geometry.JPG

    The angles for oriented ellipsoid.

Validation
----------

Validation of our code was done by comparing the output of the 1D model
to the output of the software provided by the NIST (Kline, 2006).
:num:`Figure ellipsoid-comparison-1d` below shows a comparison of
the 1D output of our model and the output of the NIST software.

.. _ellipsoid-comparison-1d:

.. figure:: img/ellipsoid_comparison_1d.jpg

    Comparison of the SasView scattering intensity for an ellipsoid
    with the output of the NIST SANS analysis software.  The parameters
    were set to: *scale* = 1.0, *rpolar* = 20 |Ang|,
    *requatorial* =400 |Ang|, *contrast* = 3e-6 |Ang^-2|,
    and *background* = 0.01 |cm^-1|.

Averaging over a distribution of orientation is done by evaluating the
equation above. Since we have no other software to compare the
implementation of the intensity for fully oriented ellipsoids, we can
compare the result of averaging our 2D output using a uniform distribution
$p(\theta,\phi) = 1.0$.  :num:`Figure #ellipsoid-comparison-2d`
shows the result of such a cross-check.

.. _ellipsoid-comparison-2d:

.. figure:: img/ellipsoid_comparison_2d.jpg

    Comparison of the intensity for uniformly distributed ellipsoids
    calculated from our 2D model and the intensity from the NIST SANS
    analysis software. The parameters used were: *scale* = 1.0,
    *rpolar* = 20 |Ang|, *requatorial* = 400 |Ang|,
    *contrast* = 3e-6 |Ang^-2|, and *background* = 0.0 |cm^-1|.

The discrepancy above *q* = 0.3 |cm^-1| is due to the way the form factors
are calculated in the c-library provided by NIST. A numerical integration
has to be performed to obtain $P(Q)$ for randomly oriented particles.
The NIST software performs that integration with a 76-point Gaussian
quadrature rule, which will become imprecise at high $Q$ where the amplitude
varies quickly as a function of $Q$. The SasView result shown has been
obtained by summing over 501 equidistant points. Our result was found
to be stable over the range of $Q$ shown for a number of points higher
than 500.

REFERENCE

L A Feigin and D I Svergun. *Structure Analysis by Small-Angle X-Ray and Neutron Scattering*, Plenum,
New York, 1987.
"""

from numpy import pi, inf

name = "ellipsoid"
title = "Ellipsoid of revolution with uniform scattering length density."

description = """\
P(q.alpha)= scale*f(q)^2 + background, where f(q)= 3*(sld
		- solvent_sld)*V*[sin(q*r(Rp,Re,alpha))
		-q*r*cos(qr(Rp,Re,alpha))]
		/[qr(Rp,Re,alpha)]^3"

     r(Rp,Re,alpha)= [Re^(2)*(sin(alpha))^2
		+ Rp^(2)*(cos(alpha))^2]^(1/2)

		sld: SLD of the ellipsoid
		solvent_sld: SLD of the solvent
		V: volume of the ellipsoid
		Rp: polar radius of the ellipsoid
		Re: equatorial radius of the ellipsoid
"""

parameters = [
#   [ "name", "units", default, [lower, upper], "type",
#     "description" ],
    [ "sld", "1e-6/Ang^2", 4, [-inf,inf], "",
      "Ellipsoid scattering length density" ],
    [ "solvent_sld", "1e-6/Ang^2", 1, [-inf,inf], "",
      "Solvent scattering length density" ],
    [ "rpolar", "Ang",  20, [0, inf], "volume",
      "Polar radius" ],
    [ "requatorial", "Ang",  400, [0, inf], "volume",
      "Equatorial radius" ],
    [ "theta", "degrees", 60, [-inf, inf], "orientation",
      "In plane angle" ],
    [ "phi", "degrees", 60, [-inf, inf], "orientation",
      "Out of plane angle" ],
    ]

source = [ "lib/J1.c", "lib/gauss76.c", "ellipsoid.c"]

def ER(rpolar, requatorial):
    import numpy as np

    ee = np.empty_like(rpolar)
    idx = rpolar > requatorial
    ee[idx] = (rpolar[idx]**2 - requatorial[idx]**2)/rpolar[idx]**2
    idx = rpolar < requatorial
    ee[idx] = (requatorial[idx]**2 - rpolar[idx]**2)/requatorial[idx]**2
    idx = rpolar == requatorial
    ee[idx] = 2*rpolar[idx]
    valid = (rpolar*requatorial != 0)
    bd = 1.0-ee[valid]
    e1 = np.sqrt(ee[valid])
    b1 = 1.0 + np.arcsin(e1)/(e1*np.sqrt(bd))
    bL = (1.0+e1)/(1.0-e1)
    b2 = 1.0 + bd/2/e1*np.log(bL)
    delta = 0.75*b1*b2

    ddd = np.zeros_like(rpolar)
    ddd[valid] = 2.0*(delta+1.0)*rpolar*requatorial**2
    return 0.5*ddd**(1.0/3.0)
