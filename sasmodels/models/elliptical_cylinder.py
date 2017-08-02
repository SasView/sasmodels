# pylint: disable=line-too-long
r"""
Definition for 2D (orientated system)
-------------------------------------

The angles $\theta$ and $\phi$ define the orientation of the axis of the
cylinder. The angle $\Psi$ is defined as the orientation of the major
axis of the ellipse with respect to the vector $Q$. A gaussian polydispersity
can be added to any of the orientation angles, and also for the minor
radius and the ratio of the ellipse radii.

.. figure:: img/elliptical_cylinder_geometry.png

   Elliptical cylinder geometry $a = r_\text{minor}$
   and $\nu = r_\text{major} / r_\text{minor}$ is the *axis_ratio*.

The function calculated is

.. math::

    I(\vec q)=\frac{1}{V_\text{cyl}}\int{d\psi}\int{d\phi}\int{
        p(\theta,\phi,\psi)F^2(\vec q,\alpha,\psi)\sin(\alpha)d\alpha}

with the functions

.. math::

    F(q,\alpha,\psi) = 2\frac{J_1(a)\sin(b)}{ab}

where

.. math::

    a = qr'\sin(\alpha)

    b = q\frac{L}{2}\cos(\alpha)

    r'=\frac{r_{minor}}{\sqrt{2}}\sqrt{(1+\nu^{2}) + (1-\nu^{2})cos(\psi)}


and the angle $\psi$ is defined as the orientation of the major axis of the
ellipse with respect to the vector $\vec q$. The angle $\alpha$ is the angle
between the axis of the cylinder and $\vec q$.


Definition for 1D (no preferred orientation)
--------------------------------------------

The form factor is averaged over all possible orientation before normalized
by the particle volume

.. math::

    P(q) = \text{scale}  <F^2> / V

To provide easy access to the orientation of the elliptical cylinder, we
define the axis of the cylinder using two angles $\theta$, $\phi$ and $\Psi$
(see :ref:`cylinder orientation <cylinder-angle-definition>`). The angle
$\Psi$ is the rotational angle around its own long_c axis.

All angle parameters are valid and given only for 2D calculation; ie, an
oriented system.

.. figure:: img/elliptical_cylinder_angle_definition.png

    Definition of angles for oriented elliptical cylinder, where axis_ratio is drawn >1,
    and angle $\Psi$ is now a rotation around the axis of the cylinder.

.. figure:: img/elliptical_cylinder_angle_projection.png

    Examples of the angles for oriented elliptical cylinders against the
    detector plane, with $\Psi$ = 0.

The $\theta$ and $\phi$ parameters to orient the cylinder only appear in the model when fitting 2d data.
On introducing "Orientational Distribution" in the angles, "distribution of theta" and "distribution of phi" parameters will
appear. These are actually rotations about the axes $\delta_1$ and $\delta_2$ of the cylinder, the $b$ and $a$ axes of the
cylinder cross section. (When $\theta = \phi = 0$ these are parallel to the $Y$ and $X$ axes of the instrument.)
The third orientation distribution, in $\psi$, is about the $c$ axis of the particle. Some experimentation may be required to
understand the 2d patterns fully. (Earlier implementations had numerical integration issues in some circumstances when orientation
distributions passed through 90 degrees, such situations, with very broad distributions, should still be approached with care.)

NB: The 2nd virial coefficient of the cylinder is calculated based on the
averaged radius $(=\sqrt{r_\text{minor}^2 * \text{axis ratio}})$ and length
values, and used as the effective radius for $S(Q)$ when $P(Q)*S(Q)$ is applied.


Validation
----------

Validation of our code was done by comparing the output of the 1D calculation
to the angular average of the output of the 2D calculation over all possible
angles.

In the 2D average, more binning in the angle $\phi$ is necessary to get the
proper result. The following figure shows the results of the averaging by
varying the number of angular bins.

.. figure:: img/elliptical_cylinder_averaging.png

    The intensities averaged from 2D over different numbers of bins and angles.

References
----------

L A Feigin and D I Svergun, *Structure Analysis by Small-Angle X-Ray and
Neutron Scattering*, Plenum, New York, (1987) [see table 3.4]

Authorship and Verification
----------------------------

* **Author:**
* **Last Modified by:**
* **Last Reviewed by:**  Richard Heenan - corrected equation in docs **Date:** December 21, 2016

"""

from numpy import pi, inf, sqrt, sin, cos

name = "elliptical_cylinder"
title = "Form factor for an elliptical cylinder."
description = """
    Form factor for an elliptical cylinder.
    See L A Feigin and D I Svergun, Structure Analysis by Small-Angle X-Ray and Neutron Scattering, Plenum, New York, (1987).
"""
category = "shape:cylinder"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["radius_minor",     "Ang",        20.0,  [0, inf],    "volume",      "Ellipse minor radius"],
              ["axis_ratio",   "",          1.5,   [1, inf],    "volume",      "Ratio of major radius over minor radius"],
              ["length",      "Ang",        400.0, [1, inf],    "volume",      "Length of the cylinder"],
              ["sld",         "1e-6/Ang^2", 4.0,   [-inf, inf], "sld",         "Cylinder scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1.0,   [-inf, inf], "sld",         "Solvent scattering length density"],
              ["theta",       "degrees",    90.0,  [-360, 360], "orientation", "cylinder axis to beam angle"],
              ["phi",         "degrees",    0,     [-360, 360], "orientation", "rotation about beam"],
              ["psi",         "degrees",    0,     [-360, 360], "orientation", "rotation about cylinder axis"]]

# pylint: enable=bad-whitespace, line-too-long

source = ["lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c", "lib/gauss20.c",
          "elliptical_cylinder.c"]

demo = dict(scale=1, background=0, radius_minor=100, axis_ratio=1.5, length=400.0,
            sld=4.0, sld_solvent=1.0, theta=10.0, phi=20, psi=30,
            theta_pd=10, phi_pd=2, psi_pd=3)

def ER(radius_minor, axis_ratio, length):
    """
        Equivalent radius
        @param radius_minor: Ellipse minor radius
        @param axis_ratio: Ratio of major radius over minor radius
        @param length: Length of the cylinder
    """
    radius = sqrt(radius_minor * radius_minor * axis_ratio)
    ddd = 0.75 * radius * (2 * radius * length
                           + (length + radius) * (length + pi * radius))
    return 0.5 * (ddd) ** (1. / 3.)

def random():
    import numpy as np
    # V = pi * radius_major * radius_minor * length;
    V = 10**np.random.uniform(3, 9)
    length = 10**np.random.uniform(1, 3)
    axis_ratio = 10**np.random.uniform(-2, 2)
    radius_minor = np.sqrt(V/length/axis_ratio)
    Vf = 10**np.random.uniform(-4, -2)
    pars = dict(
        #background=0, sld=0, sld_solvent=1,
        scale=1e9*Vf/V,
        length=length,
        radius_minor=radius_minor,
        axis_ratio=axis_ratio,
    )
    return pars

q = 0.1
# april 6 2017, rkh added a 2d unit test, NOT READY YET pull #890 branch assume correct!
qx = q*cos(pi/6.0)
qy = q*sin(pi/6.0)

tests = [
    [{'radius_minor': 20.0, 'axis_ratio': 1.5, 'length':400.0}, 'ER', 79.89245454155024],
    [{'radius_minor': 20.0, 'axis_ratio': 1.2, 'length':300.0}, 'VR', 1],

    # The SasView test result was 0.00169, with a background of 0.001
    [{'radius_minor': 20.0, 'axis_ratio': 1.5, 'sld': 4.0, 'length':400.0,
      'sld_solvent':1.0, 'background':0.0},
     0.001, 675.504402],
#    [{'theta':80., 'phi':10.}, (qx, qy), 7.88866563001 ],
]
