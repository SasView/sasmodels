# triaxial ellipsoid model
# Note: model title and parameter table are inserted automatically
r"""
All three axes are of different lengths with $R_a \leq R_b \leq R_c$
**Users should maintain this inequality for all calculations**.

.. math::

    P(q) = \text{scale} V \left< F^2(q) \right> + \text{background}

where the volume $V = 4/3 \pi R_a R_b R_c$, and the averaging
$\left<\ldots\right>$ is applied over all orientations for 1D.

.. figure:: img/triaxial_ellipsoid_geometry.jpg

    Ellipsoid schematic.

Definition
----------

The form factor calculated is

.. math::

    P(q) = \frac{\text{scale}}{V}\int_0^1\int_0^1
        \Phi^2(qR_a^2\cos^2( \pi x/2) + qR_b^2\sin^2(\pi y/2)(1-y^2) + R_c^2y^2)
        dx dy

where

.. math::

    \Phi(u) = 3 u^{-3} (\sin u - u \cos u)

To provide easy access to the orientation of the triaxial ellipsoid,
we define the axis of the cylinder using the angles $\theta$, $\phi$
and $\psi$. These angles are defined on
:numref:`triaxial-ellipsoid-angles`.
The angle $\psi$ is the rotational angle around its own $c$ axis
against the $q$ plane. For example, $\psi = 0$ when the
$a$ axis is parallel to the $x$ axis of the detector.

.. _triaxial-ellipsoid-angles:

.. figure:: img/triaxial_ellipsoid_angle_projection.jpg

    The angles for oriented ellipsoid.

The radius-of-gyration for this system is  $R_g^2 = (R_a R_b R_c)^2/5$.

The contrast is defined as SLD(ellipsoid) - SLD(solvent).  In the
parameters, $R_a$ is the minor equatorial radius, $R_b$ is the major
equatorial radius, and $R_c$ is the polar radius of the ellipsoid.

NB: The 2nd virial coefficient of the triaxial solid ellipsoid is
calculated based on the polar radius $R_p = R_c$ and equatorial
radius $R_e = \sqrt{R_a R_b}$, and used as the effective radius for
$S(q)$ when $P(q) \cdot S(q)$ is applied.

Validation
----------

Validation of our code was done by comparing the output of the
1D calculation to the angular average of the output of 2D calculation
over all possible angles.


References
----------

L A Feigin and D I Svergun, *Structure Analysis by Small-Angle X-Ray
and Neutron Scattering*, Plenum, New York, 1987.
"""

from numpy import inf

name = "triaxial_ellipsoid"
title = "Ellipsoid of uniform scattering length density with three independent axes."

description = """\
Note: During fitting ensure that the inequality ra<rb<rc is not
	violated. Otherwise the calculation will
	not be correct.
"""
category = "shape:ellipsoid"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 4, [-inf, inf], "sld",
               "Ellipsoid scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["radius_equat_minor", "Ang", 20, [0, inf], "volume",
               "Minor equatorial radius"],
              ["radius_equat_major", "Ang", 400, [0, inf], "volume",
               "Major equatorial radius"],
              ["radius_polar", "Ang", 10, [0, inf], "volume",
               "Polar radius"],
              ["theta", "degrees", 60, [-inf, inf], "orientation",
               "In plane angle"],
              ["phi", "degrees", 60, [-inf, inf], "orientation",
               "Out of plane angle"],
              ["psi", "degrees", 60, [-inf, inf], "orientation",
               "Out of plane angle"],
             ]

source = ["lib/sph_j1c.c", "lib/gauss76.c", "triaxial_ellipsoid.c"]

def ER(radius_equat_minor, radius_equat_major, radius_polar):
    """
        Returns the effective radius used in the S*P calculation
    """
    import numpy as np
    from .ellipsoid import ER as ellipsoid_ER
    return ellipsoid_ER(radius_polar, np.sqrt(radius_equat_minor * radius_equat_major))

demo = dict(scale=1, background=0,
            sld=6, sld_solvent=1,
            theta=30, phi=15, psi=5,
            radius_equat_minor=25, radius_equat_major=36, radius_polar=50,
            radius_equat_minor_pd=0, radius_equat_minor_pd_n=1,
            radius_equat_major_pd=0, radius_equat_major_pd_n=1,
            radius_polar_pd=.2, radius_polar_pd_n=30,
            theta_pd=15, theta_pd_n=45,
            phi_pd=15, phi_pd_n=1,
            psi_pd=15, psi_pd_n=1)
