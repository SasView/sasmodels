# triaxial ellipsoid model
# Note: model title and parameter table are inserted automatically
r"""
All three axes are of different lengths with $R_a \le R_b <= R_c$
**Users should maintain this inequality for all calculations**.

.. math::

    P(Q) = \text{scale} V \left< F^2(Q) \right> + \text{background}

where the volume $V = 4/3 \pi R_a R_b R_c$, and the averaging
$\left< \cdots \right>$ is applied over all orientations for 1D.

.. figure:: img/triaxial_ellipsoid_geometry.jpg

    Ellipsoid schematic.

The returned value is in units of |cm^-1|, on absolute scale.

Definition
----------

The form factor calculated is

.. math::

    P(Q) = \frac{\text{scale}}{V}\int_0^1\int_0^1
        \Phi^2(QR_a^2\cos^2( \pi x/2) + QR_b^2\sin^2(\pi y/2)(1-y^2) + c^2y^2)
        dx dy

where

.. math::

    \Phi(u) = 3 u^{-3} (\sin u - u \cos u)

To provide easy access to the orientation of the triaxial ellipsoid,
we define the axis of the cylinder using the angles $\theta$, $\phi$
and $\psi$. These angles are defined on
:num:`figure #triaxial-ellipsoid-angles`.
The angle $\psi$ is the rotational angle around its own $c$ axis
against the $Q$ plane. For example, $\psi = 0$ when the
$a$ axis is parallel to the $x$ axis of the detector.

.. _triaxial-ellipsoid-angles:

.. figure:: img/triaxial_ellipsoid_angles.jpg

    The angles for oriented ellipsoid.

The radius-of-gyration for this system is  $R_g^2 = (R_a R_b R_c)^2/5$.

The contrast is defined as SLD(ellipsoid) - SLD(solvent).  In the
parameters, *a* is the minor equatorial radius, *b* is the major
equatorial radius, and c is the polar radius of the ellipsoid.

NB: The 2nd virial coefficient of the triaxial solid ellipsoid is
calculated based on the polar radius $R_p = R_c$ and equatorial
radius $R_e = \sqrt{R_a R_b}$, and used as the effective radius for
$S(Q)$ when $P(Q) \cdot S(Q)$ is applied.

.. figure:: img/triaxial_ellipsoid_1d.jpg

    1D plot using the default values (w/1000 data point).

Validation
----------

Validation of our code was done by comparing the output of the
1D calculation to the angular average of the output of 2D calculation
over all possible angles.
:num:`Figure #triaxial-ellipsoid-comparison` shows the comparison where
the solid dot refers to averaged 2D while the line represents the
result of 1D calculation (for 2D averaging, 76, 180, and 76 points
are taken for the angles of $\theta$, $\phi$, and $\psi$ respectively).

.. _triaxial-ellipsoid-comparison:

.. figure:: img/triaxial_ellipsoid_comparison.png

    Comparison between 1D and averaged 2D.

Our model uses the form factor calculations implemented in a c-library provided by the NIST Center for Neutron Research
(Kline, 2006)

REFERENCE

L A Feigin and D I Svergun, *Structure Analysis by Small-Angle X-Ray and Neutron Scattering*, Plenum,
New York, 1987.
"""

from numpy import pi, inf

name = "triaxial_ellipsoid"
title = "Ellipsoid of uniform scattering length density with three independent axes."

description = """\
Note: During fitting ensure that the inequality ra<rb<rc is not
	violated. Otherwise the calculation will
	not be correct.
"""

parameters = [
#   [ "name", "units", default, [lower, upper], "type",
#     "description" ],
    [ "sld", "1e-6/Ang^2", 4, [-inf,inf], "",
      "Ellipsoid scattering length density" ],
    [ "solvent_sld", "1e-6/Ang^2", 1, [-inf,inf], "",
      "Solvent scattering length density" ],
    [ "req_minor", "Ang",  20, [0, inf], "volume",
      "Minor equitorial radius" ],
    [ "req_major", "Ang",  400, [0, inf], "volume",
      "Major equatorial radius" ],
    [ "rpolar", "Ang",  10, [0, inf], "volume",
      "Polar radius" ],
    [ "theta", "degrees", 60, [-inf, inf], "orientation",
      "In plane angle" ],
    [ "phi", "degrees", 60, [-inf, inf], "orientation",
      "Out of plane angle" ],
    [ "psi", "degrees", 60, [-inf, inf], "orientation",
      "Out of plane angle" ],
    ]

source = [ "lib/J1.c", "lib/gauss76.c", "triaxial_ellipsoid.c"]

def ER(req_minor, req_major, rpolar):
    import numpy as np
    from .ellipsoid import ER as ellipsoid_ER
    return ellipsoid_ER(rpolar, np.sqrt(req_minor*req_major))

