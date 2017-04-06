# parallelepiped model
# Note: model title and parameter table are inserted automatically
r"""
The form factor is normalized by the particle volume.
For information about polarised and magnetic scattering, see
the :ref:`magnetism` documentation.

Definition
----------

| This model calculates the scattering from a rectangular parallelepiped
| (\:numref:`parallelepiped-image`\).
| If you need to apply polydispersity, see also :ref:`rectangular-prism`.

.. _parallelepiped-image:

.. figure:: img/parallelepiped_geometry.jpg

   Parallelepiped with the corresponding definition of sides.

.. note::

   The edge of the solid must satisfy the condition that $A < B < C$.
   This requirement is not enforced in the model, so it is up to the
   user to check this during the analysis.

The 1D scattering intensity $I(q)$ is calculated as:

.. Comment by Miguel Gonzalez:
   I am modifying the original text because I find the notation a little bit
   confusing. I think that in most textbooks/papers, the notation P(Q) is
   used for the form factor (adim, P(Q=0)=1), although F(q) seems also to
   be used. But here (as for many other models), P(q) is used to represent
   the scattering intensity (in cm-1 normally). It would be good to agree on
   a common notation.

.. math::

    I(q) = \frac{\text{scale}}{V} (\Delta\rho \cdot V)^2
           \left< P(q, \alpha) \right> + \text{background}

where the volume $V = A B C$, the contrast is defined as
$\Delta\rho = \rho_\text{p} - \rho_\text{solvent}$,
$P(q, \alpha)$ is the form factor corresponding to a parallelepiped oriented
at an angle $\alpha$ (angle between the long axis C and $\vec q$),
and the averaging $\left<\ldots\right>$ is applied over all orientations.

Assuming $a = A/B < 1$, $b = B /B = 1$, and $c = C/B > 1$, the
form factor is given by (Mittelbach and Porod, 1961)

.. math::

    P(q, \alpha) = \int_0^1 \phi_Q\left(\mu \sqrt{1-\sigma^2},a\right)
        \left[S(\mu c \sigma/2)\right]^2 d\sigma

with

.. math::

    \phi_Q(\mu,a) &= \int_0^1
        \left\{S\left[\frac{\mu}{2}\cos\left(\frac{\pi}{2}u\right)\right]
               S\left[\frac{\mu a}{2}\sin\left(\frac{\pi}{2}u\right)\right]
               \right\}^2 du

    S(x) &= \frac{\sin x}{x}

    \mu &= qB


The scattering intensity per unit volume is returned in units of |cm^-1|.

NB: The 2nd virial coefficient of the parallelepiped is calculated based on
the averaged effective radius $(=\sqrt{A B / \pi})$ and
length $(= C)$ values, and used as the effective radius for
$S(q)$ when $P(q) \cdot S(q)$ is applied.

To provide easy access to the orientation of the parallelepiped, we define
three angles $\theta$, $\phi$ and $\Psi$. The definition of $\theta$ and
$\phi$ is the same as for the cylinder model (see also figures below).

.. Comment by Miguel Gonzalez:
   The following text has been commented because I think there are two
   mistakes. Psi is the rotational angle around C (but I cannot understand
   what it means against the q plane) and psi=0 corresponds to a||x and b||y.

   The angle $\Psi$ is the rotational angle around the $C$ axis against
   the $q$ plane. For example, $\Psi = 0$ when the $B$ axis is parallel
   to the $x$-axis of the detector.

The angle $\Psi$ is the rotational angle around the $C$ axis.
For $\theta = 0$ and $\phi = 0$, $\Psi = 0$ corresponds to the $B$ axis
oriented parallel to the y-axis of the detector with $A$ along the z-axis.
For other $\theta$, $\phi$ values, the parallelepiped has to be first rotated
$\theta$ degrees around $z$ and $\phi$ degrees around $y$,
before doing a final rotation of $\Psi$ degrees around the resulting $C$ to
obtain the final orientation of the parallelepiped.
For example, for $\theta = 0$ and $\phi = 90$, we have that $\Psi = 0$
corresponds to $A$ along $x$ and $B$ along $y$,
while for $\theta = 90$ and $\phi = 0$, $\Psi = 0$ corresponds to
$A$ along $z$ and $B$ along $x$.

.. _parallelepiped-orientation:

.. figure:: img/parallelepiped_angle_definition.jpg

    Definition of the angles for oriented parallelepipeds.

.. figure:: img/parallelepiped_angle_projection.jpg

    Examples of the angles for oriented parallelepipeds against the
    detector plane.

For a given orientation of the parallelepiped, the 2D form factor is
calculated as

.. math::

    P(q_x, q_y) = \left[\frac{\sin(qA\cos\alpha/2)}{(qA\cos\alpha/2)}\right]^2
                  \left[\frac{\sin(qB\cos\beta/2)}{(qB\cos\beta/2)}\right]^2
                  \left[\frac{\sin(qC\cos\gamma/2)}{(qC\cos\gamma/2)}\right]^2

with

.. math::

    \cos\alpha &= \hat A \cdot \hat q,

    \cos\beta  &= \hat B \cdot \hat q,

    \cos\gamma &= \hat C \cdot \hat q

and the scattering intensity as:

.. math::

    I(q_x, q_y) = \frac{\text{scale}}{V} V^2 \Delta\rho^2 P(q_x, q_y)
            + \text{background}

.. Comment by Miguel Gonzalez:
   This reflects the logic of the code, as in parallelepiped.c the call
   to _pkernel returns $P(q_x, q_y)$ and then this is multiplied by
   $V^2 * (\Delta \rho)^2$. And finally outside parallelepiped.c it will be
   multiplied by scale, normalized by $V$ and the background added. But
   mathematically it makes more sense to write
   $I(q_x, q_y) = \text{scale} V \Delta\rho^2 P(q_x, q_y) + \text{background}$,
   with scale being the volume fraction.


Validation
----------

Validation of the code was done by comparing the output of the 1D calculation
to the angular average of the output of a 2D calculation over all possible
angles.

This model is based on form factor calculations implemented in a c-library
provided by the NIST Center for Neutron Research (Kline, 2006).

References
----------

P Mittelbach and G Porod, *Acta Physica Austriaca*, 14 (1961) 185-211

R Nayuk and K Huber, *Z. Phys. Chem.*, 226 (2012) 837-854
"""

import numpy as np
from numpy import pi, inf, sqrt

name = "parallelepiped"
title = "Rectangular parallelepiped with uniform scattering length density."
description = """
    I(q)= scale*V*(sld - sld_solvent)^2*P(q,alpha)+background
        P(q,alpha) = integral from 0 to 1 of ...
           phi(mu*sqrt(1-sigma^2),a) * S(mu*c*sigma/2)^2 * dsigma
        with
            phi(mu,a) = integral from 0 to 1 of ..
            (S((mu/2)*cos(pi*u/2))*S((mu*a/2)*sin(pi*u/2)))^2 * du
            S(x) = sin(x)/x
            mu = q*B
        V: Volume of the rectangular parallelepiped
        alpha: angle between the long axis of the 
            parallelepiped and the q-vector for 1D
"""
category = "shape:parallelepiped"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 4, [-inf, inf], "sld",
               "Parallelepiped scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["length_a", "Ang", 35, [0, inf], "volume",
               "Shorter side of the parallelepiped"],
              ["length_b", "Ang", 75, [0, inf], "volume",
               "Second side of the parallelepiped"],
              ["length_c", "Ang", 400, [0, inf], "volume",
               "Larger side of the parallelepiped"],
              ["theta", "degrees", 60, [-inf, inf], "orientation",
               "In plane angle"],
              ["phi", "degrees", 60, [-inf, inf], "orientation",
               "Out of plane angle"],
              ["psi", "degrees", 60, [-inf, inf], "orientation",
               "Rotation angle around its own c axis against q plane"],
             ]

source = ["lib/gauss76.c", "parallelepiped.c"]

def ER(length_a, length_b, length_c):
    """
    Return effective radius (ER) for P(q)*S(q)
    """
    # now that axes can be in any size order, need to sort a,b,c where a~b and c is either much smaller 
    # or much larger
    abc = np.vstack((length_a, length_b, length_c))
    abc = np.sort(abc, axis=0)
    selector = (abc[1] - abc[0]) > (abc[2] - abc[1])
    length = np.where(selector, abc[0], abc[2])
    # surface average radius (rough approximation)
    radius = np.sqrt(np.where(~selector, abc[0]*abc[1], abc[1]*abc[2]) / pi)

    ddd = 0.75 * radius * (2*radius*length + (length + radius)*(length + pi*radius))
    return 0.5 * (ddd) ** (1. / 3.)

# VR defaults to 1.0

# parameters for demo
demo = dict(scale=1, background=0,
            sld=6.3, sld_solvent=1.0,
            length_a=35, length_b=75, length_c=400,
            theta=45, phi=30, psi=15,
            length_a_pd=0.1, length_a_pd_n=10,
            length_b_pd=0.1, length_b_pd_n=1,
            length_c_pd=0.1, length_c_pd_n=1,
            theta_pd=10, theta_pd_n=1,
            phi_pd=10, phi_pd_n=1,
            psi_pd=10, psi_pd_n=10)

qx, qy = 0.2 * np.cos(2.5), 0.2 * np.sin(2.5)
tests = [[{}, 0.2, 0.17758004974],
         [{}, [0.2], [0.17758004974]],
         [{'theta':10.0, 'phi':10.0}, (qx, qy), 0.00560296014],
         [{'theta':10.0, 'phi':10.0}, [(qx, qy)], [0.00560296014]],
        ]
del qx, qy  # not necessary to delete, but cleaner
