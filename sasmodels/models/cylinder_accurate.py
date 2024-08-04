# cylinder model
# Note: model title and parameter table are inserted automatically
r"""

For information about polarised and magnetic scattering, see
the :ref:`magnetism` documentation.

Definition
----------

The output of the 2D scattering intensity function for oriented cylinders is
given by (Guinier, 1955)

.. math::

    I(q,\alpha) = \frac{\text{scale}}{V} F^2(q,\alpha) + \text{background}

where

.. math::

    F(q,\alpha) = 2 (\Delta \rho) V
           \frac{\sin \left(\tfrac12 qL\cos\alpha \right)}
                {\tfrac12 qL \cos \alpha}
           \frac{J_1 \left(q R \sin \alpha\right)}{q R \sin \alpha}

and $\alpha$ is the angle between the axis of the cylinder and $\vec q$,
$V =\pi R^2L$ is the volume of the cylinder, $L$ is the length of the cylinder,
$R$ is the radius of the cylinder, and $\Delta\rho$ (contrast) is the
scattering length density difference between the scatterer and the solvent.
$J_1$ is the first order Bessel function.

For randomly oriented particles:

.. math::

    P(q)=F^2(q)=\int_{0}^{\pi/2}{F^2(q,\alpha)\sin(\alpha)d\alpha}


The output of the 1D scattering intensity function for randomly oriented
cylinders is thus given by

.. math::

    I(q) = \frac{\text{scale}}{V}
        \int_0^{\pi/2} F^2(q,\alpha) \sin \alpha\ d\alpha + \text{background}


NB: The 2nd virial coefficient of the cylinder is calculated based on the
radius and length values, and used as the effective radius for $S(q)$
when $P(q) \cdot S(q)$ is applied.

For 2d scattering from oriented cylinders, we define the direction of the
axis of the cylinder using two angles $\theta$ (note this is not the same as
the scattering angle used in q) and $\phi$. Those angles are defined in
:numref:`cylinder-angle-definition` , for further details see
:ref:`orientation`.

.. _cylinder-angle-definition:

.. figure:: img/cylinder_angle_definition.png

    Angles $\theta$ and $\phi$ orient the cylinder relative to the beam line
    coordinates, where the beam is along the $z$ axis. Rotation $\theta$,
    initially in the $xz$ plane, is carried out first, then rotation $\phi$
    about the $z$ axis. Orientation distributions are described as rotations
    about two perpendicular axes $\delta_1$ and $\delta_2$ in the frame of
    the cylinder itself, which when $\theta = \phi = 0$ are parallel to the
    $Y$ and $X$ axes.

.. figure:: img/cylinder_angle_projection.png

    Examples for oriented cylinders.

The $\theta$ and $\phi$ parameters to orient the cylinder only appear in the
model when fitting 2d data.

Validation
----------

Validation of the code was done by comparing the output of the 1D model
to the output of the software provided by the NIST (Kline, 2006).
The implementation of the intensity for fully oriented cylinders was done
by averaging over a uniform distribution of orientations using

.. math::

    P(q) = \int_0^{\pi/2} d\phi
        \int_0^\pi p(\theta) P_0(q,\theta) \sin \theta\ d\theta


where $p(\theta,\phi) = 1$ is the probability distribution for the orientation
and $P_0(q,\theta)$ is the scattering intensity for the fully oriented
system, and then comparing to the 1D result.

References
----------

#.  J. Pedersen, *Adv. Colloid Interface Sci.*, 70 (1997) 171-210
#.  G. Fournet, *Bull. Soc. Fr. Mineral. Cristallogr.*, 74 (1951) 39-113
#.  L. Onsager, *Ann. New York Acad. Sci.*, 51 (1949) 627-659

Authorship and Verification
----------------------------

* **Author:**
* **Last Modified by:** Paul Butler (docs only) November 10, 2022
* **Last Reviewed by:**
"""

import numpy as np  # type: ignore
from numpy import pi, inf  # type: ignore
from scipy.special import j1  # type: ignore
import os
from .libpy import tuned_quad_integrate, load_tuned_quad_h5, kronrod_points_dict

name = "cylinder_accurate"
title = "Right circular cylinder with uniform scattering length density."
description = """
    Identical to cylinder model but with a more accurate integration scheme.
     f(q,alpha) = 2*(sld - sld_solvent)*V*sin(qLcos(alpha)/2))
                /[qLcos(alpha)/2]*J1(qRsin(alpha))/[qRsin(alpha)]

            P(q,alpha)= scale/V*f(q,alpha)^(2)+background
            V: Volume of the cylinder
            R: Radius of the cylinder
            L: Length of the cylinder
            J1: The Bessel function
            alpha: angle between the axis of the
            cylinder and the q-vector for 1D
            :the ouput is P(q)=scale/V*integral
            from pi/2 to zero of...
            f(q,alpha)^(2)*sin(alpha)*dalpha + background
"""
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#             [ "name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["sld", "1e-6/Ang^2", 4, [-inf, inf], "sld", "Cylinder scattering length density"],
    ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld", "Solvent scattering length density"],
    ["radius", "Ang", 20, [0, inf], "volume", "Cylinder radius"],
    ["length", "Ang", 400, [0, inf], "volume", "Cylinder length"],
    ["rtol", "", 1e-4, [0, inf], "", "Relative tolerance for integration"],
    # ["theta", "degrees", 60, [-360, 360], "orientation", "cylinder axis to beam angle"],
    # ["phi", "degrees", 60, [-360, 360], "orientation", "rotation about beam"],
]
cylinder_accurate_quad = load_tuned_quad_h5("cylinder_accurate", kronrod_points_dict)

# pylint: enable=bad-whitespace, line-too-long

def J1x(x):
    """
    Calculate the Bessel function J1(x) for scalar x.
    """
    return np.where(x != 0, j1(x) / x, 0.5)


def Fq(alpha, params):
    A = params["A"]
    B = params["B"]
    print("Working Direction from Fq", os.getcwd())
    return (np.sinc(A * np.cos(alpha) / np.pi) * J1x(B * np.sin(alpha))) ** 2 * np.sin(alpha)


def Iq(q, sld, sld_solvent, radius, length, rtol):

    # , theta, phi
    V = pi * radius ** 2 * length
    # fq_coeff = 2 * (sld - sld_solvent) * V
    s = (sld-sld_solvent) * V * 2
    fq_params = {"rtol": rtol, "A": q * length / 2, "B": q * radius}
    # print("Working Direction from Iq", os.getcwd())
    IntFq = tuned_quad_integrate(cylinder_accurate_quad, Fq, 0, pi / 2, fq_params)
    IntFq *= s**2*1e-4
    return IntFq
    # return 1
