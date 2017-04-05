# triaxial ellipsoid model
# Note: model title and parameter table are inserted automatically
r"""
Definition
----------

.. figure:: img/triaxial_ellipsoid_geometry.jpg

    Ellipsoid with $R_a$ as *radius_equat_minor*, $R_b$ as *radius_equat_major*
    and $R_c$ as *radius_polar*.

Given an ellipsoid

.. math::

    \frac{X^2}{R_a^2} + \frac{Y^2}{R_b^2} + \frac{Z^2}{R_c^2} = 1

the scattering for randomly oriented particles is defined by the average over all orientations $\Omega$ of:

.. math::

    P(q) = \text{scale}(\Delta\rho)^2\frac{V}{4 \pi}\int_\Omega \Phi^2(qr) d\Omega + \text{background}

where

.. math::

    \Phi(qr) &= 3 j_1(qr)/qr = 3 (\sin qr - qr \cos qr)/(qr)^3 \\
    r^2 &= R_a^2e^2 + R_b^2f^2 + R_c^2g^2 \\
    V &= \tfrac{4}{3} \pi R_a R_b R_c

The $e$, $f$ and $g$ terms are the projections of the orientation vector on $X$,
$Y$ and $Z$ respectively.  Keeping the orientation fixed at the canonical
axes, we can integrate over the incident direction using polar angle
$-\pi/2 \le \gamma \le \pi/2$ and equatorial angle $0 \le \phi \le 2\pi$
(as defined in ref [1]),

 .. math::

     \langle\Phi^2\rangle = \int_0^{2\pi} \int_{-\pi/2}^{\pi/2} \Phi^2(qr) \cos \gamma\,d\gamma d\phi

with $e = \cos\gamma \sin\phi$, $f = \cos\gamma \cos\phi$ and $g = \sin\gamma$.
A little algebra yields

.. math::

    r^2 = b^2(p_a \sin^2 \phi \cos^2 \gamma + 1 + p_c \sin^2 \gamma)

for

.. math::

    p_a = \frac{a^2}{b^2} - 1 \text{ and } p_c = \frac{c^2}{b^2} - 1

Due to symmetry, the ranges can be restricted to a single quadrant
$0 \le \gamma \le \pi/2$ and $0 \le \phi \le \pi/2$, scaling the resulting
integral by 8. The computation is done using the substitution $u = \sin\gamma$,
$du = \cos\gamma\,d\gamma$, giving

.. math::

    \langle\Phi^2\rangle &= 8 \int_0^{\pi/2} \int_0^1 \Phi^2(qr) du d\phi \\
    r^2 &= b^2(p_a \sin^2(\phi)(1 - u^2) + 1 + p_c u^2)

To provide easy access to the orientation of the triaxial ellipsoid,
we define the axis of the cylinder using the angles $\theta$, $\phi$
and $\psi$. These angles are defined on
:numref:`triaxial-ellipsoid-angles` .
The angle $\psi$ is the rotational angle around its own $c$ axis
against the $q$ plane. For example, $\psi = 0$ when the
$a$ axis is parallel to the $x$ axis of the detector.

.. _triaxial-ellipsoid-angles:

.. figure:: img/triaxial_ellipsoid_angle_projection.jpg

    The angles for oriented ellipsoid.

The radius-of-gyration for this system is  $R_g^2 = (R_a R_b R_c)^2/5$.

The contrast $\Delta\rho$ is defined as SLD(ellipsoid) - SLD(solvent).  In the
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

[1] Finnigan, J.A., Jacobs, D.J., 1971.
*Light scattering by ellipsoidal particles in solution*,
J. Phys. D: Appl. Phys. 4, 72-77. doi:10.1088/0022-3727/4/1/310

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Paul Kienzle (improved calculation) **Date:** April 4, 2017
* **Last Reviewed by:** Paul Kienzle &Richard Heenan **Date:**  April 4, 2017

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
               "Minor equatorial radius, Ra"],
              ["radius_equat_major", "Ang", 400, [0, inf], "volume",
               "Major equatorial radius, Rb"],
              ["radius_polar", "Ang", 10, [0, inf], "volume",
               "Polar radius, Rc"],
              ["theta", "degrees", 60, [-inf, inf], "orientation",
               "In plane angle"],
              ["phi", "degrees", 60, [-inf, inf], "orientation",
               "Out of plane angle"],
              ["psi", "degrees", 60, [-inf, inf], "orientation",
               "Out of plane angle"],
             ]

source = ["lib/sas_3j1x_x.c", "lib/gauss76.c", "triaxial_ellipsoid.c"]

def ER(radius_equat_minor, radius_equat_major, radius_polar):
    """
        Returns the effective radius used in the S*P calculation
    """
    import numpy as np
    from .ellipsoid import ER as ellipsoid_ER
     # now that radii can be in any size order, radii need sorting a,b,c where a~b and c is either much smaller or much larger
     # also need some unit tests!
    
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
