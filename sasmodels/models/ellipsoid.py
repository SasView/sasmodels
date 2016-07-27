# ellipsoid model
# Note: model title and parameter table are inserted automatically
r"""
The form factor is normalized by the particle volume 

Definition
----------

The output of the 2D scattering intensity function for oriented ellipsoids
is given by (Feigin, 1987)

.. math::

    P(q,\alpha) = \frac{\text{scale}}{V} F^2(q,\alpha) + \text{background}

where

.. math::

    F(q,\alpha) = \frac{3 \Delta \rho V (\sin[qr(R_p,R_e,\alpha)]
                - \cos[qr(R_p,R_e,\alpha)])}
                {[qr(R_p,R_e,\alpha)]^3}

and

.. math::

    r(R_p,R_e,\alpha) = \left[ R_e^2 \sin^2 \alpha
        + R_p^2 \cos^2 \alpha \right]^{1/2}


$\alpha$ is the angle between the axis of the ellipsoid and $\vec q$,
$V = (4/3)\pi R_pR_e^2$ is the volume of the ellipsoid , $R_p$ is the polar radius along the
rotational axis of the ellipsoid, $R_e$ is the equatorial radius perpendicular
to the rotational axis of the ellipsoid and $\Delta \rho$ (contrast) is the
scattering length density difference between the scatterer and the solvent.

To provide easy access to the orientation of the ellipsoid, we define
the rotation axis of the ellipsoid using two angles $\theta$ and $\phi$.
These angles are defined in the
:ref:`cylinder orientation figure <cylinder-angle-definition>`.
For the ellipsoid, $\theta$ is the angle between the rotational axis
and the $z$ -axis.

NB: The 2nd virial coefficient of the solid ellipsoid is calculated based
on the $R_p$ and $R_e$ values, and used as the effective radius for
$S(q)$ when $P(q) \cdot S(q)$ is applied.


The $\theta$ and $\phi$ parameters are not used for the 1D output.

.. _ellipsoid-geometry:

.. figure:: img/ellipsoid_angle_projection.jpg

    The angles for oriented ellipsoid, shown here as oblate, $a$ = $R_p$ and $b$ = $R_e$

Validation
----------

Validation of the code was done by comparing the output of the 1D model
to the output of the software provided by the NIST (Kline, 2006).

The implementation of the intensity for fully oriented ellipsoids was
validated by averaging the 2D output using a uniform distribution
$p(\theta,\phi) = 1.0$ and comparing with the output of the 1D calculation.


.. _ellipsoid-comparison-2d:

.. figure:: img/ellipsoid_comparison_2d.jpg

    Comparison of the intensity for uniformly distributed ellipsoids
    calculated from our 2D model and the intensity from the NIST SANS
    analysis software. The parameters used were: *scale* = 1.0,
    *r_polar* = 20 |Ang|, *r_equatorial* = 400 |Ang|,
    *contrast* = 3e-6 |Ang^-2|, and *background* = 0.0 |cm^-1|.

The discrepancy above $q$ = 0.3 |cm^-1| is due to the way the form factors
are calculated in the c-library provided by NIST. A numerical integration
has to be performed to obtain $P(q)$ for randomly oriented particles.
The NIST software performs that integration with a 76-point Gaussian
quadrature rule, which will become imprecise at high $q$ where the amplitude
varies quickly as a function of $q$. The SasView result shown has been
obtained by summing over 501 equidistant points. Our result was found
to be stable over the range of $q$ shown for a number of points higher
than 500.

References
----------

L A Feigin and D I Svergun.
*Structure Analysis by Small-Angle X-Ray and Neutron Scattering*,
Plenum Press, New York, 1987.
"""

from numpy import inf

name = "ellipsoid"
title = "Ellipsoid of revolution with uniform scattering length density."

description = """\
P(q.alpha)= scale*f(q)^2 + background, where f(q)= 3*(sld
        - sld_solvent)*V*[sin(q*r(Rp,Re,alpha))
        -q*r*cos(qr(Rp,Re,alpha))]
        /[qr(Rp,Re,alpha)]^3"

     r(Rp,Re,alpha)= [Re^(2)*(sin(alpha))^2
        + Rp^(2)*(cos(alpha))^2]^(1/2)

        sld: SLD of the ellipsoid
        sld_solvent: SLD of the solvent
        V: volume of the ellipsoid
        Rp: polar radius of the ellipsoid
        Re: equatorial radius of the ellipsoid
"""
category = "shape:ellipsoid"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 4, [-inf, inf], "sld",
               "Ellipsoid scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["r_polar", "Ang", 20, [0, inf], "volume",
               "Polar radius"],
              ["r_equatorial", "Ang", 400, [0, inf], "volume",
               "Equatorial radius"],
              ["theta", "degrees", 60, [-inf, inf], "orientation",
               "In plane angle"],
              ["phi", "degrees", 60, [-inf, inf], "orientation",
               "Out of plane angle"],
             ]

source = ["lib/sph_j1c.c", "lib/gauss76.c", "ellipsoid.c"]

def ER(r_polar, r_equatorial):
    import numpy as np

    ee = np.empty_like(r_polar)
    idx = r_polar > r_equatorial
    ee[idx] = (r_polar[idx] ** 2 - r_equatorial[idx] ** 2) / r_polar[idx] ** 2
    idx = r_polar < r_equatorial
    ee[idx] = (r_equatorial[idx] ** 2 - r_polar[idx] ** 2) / r_equatorial[idx] ** 2
    idx = r_polar == r_equatorial
    ee[idx] = 2 * r_polar[idx]
    valid = (r_polar * r_equatorial != 0)
    bd = 1.0 - ee[valid]
    e1 = np.sqrt(ee[valid])
    b1 = 1.0 + np.arcsin(e1) / (e1 * np.sqrt(bd))
    bL = (1.0 + e1) / (1.0 - e1)
    b2 = 1.0 + bd / 2 / e1 * np.log(bL)
    delta = 0.75 * b1 * b2

    ddd = np.zeros_like(r_polar)
    ddd[valid] = 2.0 * (delta + 1.0) * r_polar * r_equatorial ** 2
    return 0.5 * ddd ** (1.0 / 3.0)


demo = dict(scale=1, background=0,
            sld=6, sld_solvent=1,
            r_polar=50, r_equatorial=30,
            theta=30, phi=15,
            r_polar_pd=.2, r_polar_pd_n=15,
            r_equatorial_pd=.2, r_equatorial_pd_n=15,
            theta_pd=15, theta_pd_n=45,
            phi_pd=15, phi_pd_n=1)
