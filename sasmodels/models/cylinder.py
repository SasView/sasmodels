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

    P(q,\alpha) = \frac{\text{scale}}{V} F^2(q,\alpha).sin(\alpha) + \text{background}

where

.. math::

    F(q,\alpha) = 2 (\Delta \rho) V
           \frac{\sin \left(\tfrac12 qL\cos\alpha \right)}
                {\tfrac12 qL \cos \alpha}
           \frac{J_1 \left(q R \sin \alpha\right)}{q R \sin \alpha}

and $\alpha$ is the angle between the axis of the cylinder and $\vec q$, $V =\pi R^2L$
is the volume of the cylinder, $L$ is the length of the cylinder, $R$ is the
radius of the cylinder, and $\Delta\rho$ (contrast) is the scattering length
density difference between the scatterer and the solvent. $J_1$ is the
first order Bessel function.

For randomly oriented particles:

.. math::

    F^2(q)=\int_{0}^{\pi/2}{F^2(q,\alpha)\sin(\alpha)d\alpha}=\int_{0}^{1}{F^2(q,u)du}


Numerical integration is simplified by a change of variable to $u = cos(\alpha)$ with
$sin(\alpha)=\sqrt{1-u^2}$.

The output of the 1D scattering intensity function for randomly oriented
cylinders is thus given by

.. math::

    P(q) = \frac{\text{scale}}{V}
        \int_0^{\pi/2} F^2(q,\alpha) \sin \alpha\ d\alpha + \text{background}


NB: The 2nd virial coefficient of the cylinder is calculated based on the
radius and length values, and used as the effective radius for $S(q)$
when $P(q) \cdot S(q)$ is applied.

For 2d scattering from oriented cylinders, we define the direction of the
axis of the cylinder using two angles $\theta$ (note this is not the
same as the scattering angle used in q) and $\phi$. Those angles
are defined in :numref:`cylinder-angle-definition` , for further details see :ref:`orientation` .

.. _cylinder-angle-definition:

.. figure:: img/cylinder_angle_definition.png

    Angles $\theta$ and $\phi$ orient the cylinder relative
    to the beam line coordinates, where the beam is along the $z$ axis. Rotation $\theta$, initially
    in the $xz$ plane, is carried out first, then rotation $\phi$ about the $z$ axis. Orientation distributions
    are described as rotations about two perpendicular axes $\delta_1$ and $\delta_2$
    in the frame of the cylinder itself, which when $\theta = \phi = 0$ are parallel to the $Y$ and $X$ axes.

.. figure:: img/cylinder_angle_projection.png

    Examples for oriented cylinders.

The $\theta$ and $\phi$ parameters to orient the cylinder only appear in the model when fitting 2d data.

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

J. S. Pedersen, Adv. Colloid Interface Sci. 70, 171-210 (1997).
G. Fournet, Bull. Soc. Fr. Mineral. Cristallogr. 74, 39-113 (1951).
"""

import numpy as np  # type: ignore
from numpy import pi, inf  # type: ignore

name = "cylinder"
title = "Right circular cylinder with uniform scattering length density."
description = """
     f(q,alpha) = 2*(sld - sld_solvent)*V*sin(qLcos(alpha)/2))
                /[qLcos(alpha)/2]*J1(qRsin(alpha))/[qRsin(alpha)]

            P(q,alpha)= scale/V*f(q,alpha)^(2)+background
            V: Volume of the cylinder
            R: Radius of the cylinder
            L: Length of the cylinder
            J1: The bessel function
            alpha: angle between the axis of the
            cylinder and the q-vector for 1D
            :the ouput is P(q)=scale/V*integral
            from pi/2 to zero of...
            f(q,alpha)^(2)*sin(alpha)*dalpha + background
"""
category = "shape:cylinder"

#             [ "name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld", "1e-6/Ang^2", 4, [-inf, inf], "sld",
               "Cylinder scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["radius", "Ang", 20, [0, inf], "volume",
               "Cylinder radius"],
              ["length", "Ang", 400, [0, inf], "volume",
               "Cylinder length"],
              ["theta", "degrees", 60, [-360, 360], "orientation",
               "cylinder axis to beam angle"],
              ["phi", "degrees", 60, [-360, 360], "orientation",
               "rotation about beam"],
             ]

source = ["lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c", "cylinder.c"]

def ER(radius, length):
    """
        Return equivalent radius (ER)
    """
    ddd = 0.75 * radius * (2 * radius * length + (length + radius) * (length + pi * radius))
    return 0.5 * (ddd) ** (1. / 3.)

def random():
    volume = 10**np.random.uniform(5, 12)
    length = 10**np.random.uniform(-2, 2)*volume**0.333
    radius = np.sqrt(volume/length/np.pi)
    pars = dict(
        #scale=1,
        #background=0,
        length=length,
        radius=radius,
    )
    return pars


# parameters for demo
demo = dict(scale=1, background=0,
            sld=6, sld_solvent=1,
            radius=20, length=300,
            theta=60, phi=60,
            radius_pd=.2, radius_pd_n=9,
            length_pd=.2, length_pd_n=10,
            theta_pd=10, theta_pd_n=5,
            phi_pd=10, phi_pd_n=5)

qx, qy = 0.2 * np.cos(2.5), 0.2 * np.sin(2.5)
# After redefinition of angles, find new tests values.  Was 10 10 in old coords
tests = [
    [{}, 0.2, 0.042761386790780453],
    [{}, [0.2], [0.042761386790780453]],
    #  new coords
    [{'theta':80.1534480601659, 'phi':10.1510817110481}, (qx, qy), 0.03514647218513852],
    [{'theta':80.1534480601659, 'phi':10.1510817110481}, [(qx, qy)], [0.03514647218513852]],
    # old coords
    #[{'theta':10.0, 'phi':10.0}, (qx, qy), 0.03514647218513852],
    #[{'theta':10.0, 'phi':10.0}, [(qx, qy)], [0.03514647218513852]],
]
del qx, qy  # not necessary to delete, but cleaner
# ADDED by:  RKH  ON: 18Mar2016 renamed sld's etc
