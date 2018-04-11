r"""
This model provides the form factor, $P(q)$, for a monodisperse hollow right
angle circular cylinder (rigid tube) where the form factor is normalized by the
volume of the tube (i.e. not by the external volume).

.. math::

    P(q) = \text{scale} \left<F^2\right>/V_\text{shell} + \text{background}

where the averaging $\left<\ldots\right>$ is applied only for the 1D calculation.

The inside and outside of the hollow cylinder are assumed have the same SLD.

Definition
----------

The 1D scattering intensity is calculated in the following way (Guinier, 1955)

.. math::

    P(q)           &= (\text{scale})V_\text{shell}\Delta\rho^2
            \int_0^{1}\Psi^2
            \left[q_z, R_\text{outer}(1-x^2)^{1/2},
                       R_\text{core}(1-x^2)^{1/2}\right]
            \left[\frac{\sin(qHx)}{qHx}\right]^2 dx \\
    \Psi[q,y,z]    &= \frac{1}{1-\gamma^2}
            \left[ \Lambda(qy) - \gamma^2\Lambda(qz) \right] \\
    \Lambda(a)     &= 2 J_1(a) / a \\
    \gamma         &= R_\text{core} / R_\text{outer} \\
    V_\text{shell} &= \pi \left(R_\text{outer}^2 - R_\text{core}^2 \right)L \\
    J_1(x)         &= (\sin(x)-x\cdot \cos(x)) / x^2

where *scale* is a scale factor, $H = L/2$ and $J_1$ is the 1st order
Bessel function.

**NB**: The 2nd virial coefficient of the cylinder is calculated
based on the outer radius and full length, which give an the effective radius
for structure factor $S(q)$ when $P(q) \cdot S(q)$ is applied.

In the parameters,the *radius* is $R_\text{core}$ while *thickness*
is $R_\text{outer} - R_\text{core}$.

To provide easy access to the orientation of the core-shell cylinder, we define
the axis of the cylinder using two angles $\theta$ and $\phi$
(see :ref:`cylinder model <cylinder-angle-definition>`).

References
----------

L A Feigin and D I Svergun, *Structure Analysis by Small-Angle X-Ray and
Neutron Scattering*, Plenum Press, New York, (1987)

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Richard Heenan **Date:** October 06, 2016
   (reparametrised to use thickness, not outer radius)
* **Last Reviewed by:** Richard Heenan **Date:** October 06, 2016
"""

import numpy as np
from numpy import pi, inf, sin, cos

name = "hollow_cylinder"
title = ""
description = """
P(q) = scale*<f*f>/Vol + background, where f is the scattering amplitude.
radius = the radius of core
thickness = the thickness of shell
length = the total length of the cylinder
sld = SLD of the shell
sld_solvent = SLD of the solvent
background = incoherent background
"""
category = "shape:cylinder"
# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["radius",      "Ang",     20.0, [0, inf],    "volume",      "Cylinder core radius"],
    ["thickness",   "Ang",     10.0, [0, inf],    "volume",      "Cylinder wall thickness"],
    ["length",      "Ang",    400.0, [0, inf],    "volume",      "Cylinder total length"],
    ["sld",         "1e-6/Ang^2",  6.3, [-inf, inf], "sld",         "Cylinder sld"],
    ["sld_solvent", "1e-6/Ang^2",  1,   [-inf, inf], "sld",         "Solvent sld"],
    ["theta",       "degrees", 90,   [-360, 360], "orientation", "Cylinder axis to beam angle"],
    ["phi",         "degrees",  0,   [-360, 360], "orientation", "Rotation about beam"],
    ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c", "hollow_cylinder.c"]

# pylint: disable=W0613
def ER(radius, thickness, length):
    """
    :param radius:      Cylinder core radius
    :param thickness:   Cylinder wall thickness
    :param length:      Cylinder length
    :return:            Effective radius
    """
    router = radius + thickness
    if router == 0 or length == 0:
        return 0.0
    len1 = router
    len2 = length/2.0
    term1 = len1*len1*2.0*len2/2.0
    term2 = 1.0 + (len2/len1)*(1.0 + 1/len2/2.0)*(1.0 + pi*len1/len2/2.0)
    ddd = 3.0*term1*term2
    diam = pow(ddd, (1.0/3.0))
    return diam

def VR(radius, thickness, length):
    """
    :param radius:      Cylinder radius
    :param thickness:   Cylinder wall thickness
    :param length:      Cylinder length
    :return:            Volf ratio for P(q)*S(q)
    """
    router = radius + thickness
    vol_core = pi*radius*radius*length
    vol_total = pi*router*router*length
    vol_shell = vol_total - vol_core
    return vol_shell, vol_total

def random():
    length = 10**np.random.uniform(1, 4.7)
    outer_radius = 10**np.random.uniform(1, 4.7)
    # Use a distribution with a preference for thin shell or thin core
    # Avoid core,shell radii < 1
    thickness = np.random.beta(0.5, 0.5)*(outer_radius-2) + 1
    radius = outer_radius - thickness
    pars = dict(
        length=length,
        radius=radius,
        thickness=thickness,
    )
    return pars

# parameters for demo
demo = dict(scale=1.0, background=0.0, length=400.0, radius=20.0,
            thickness=10, sld=6.3, sld_solvent=1, theta=90, phi=0,
            thickness_pd=0.2, thickness_pd_n=9,
            length_pd=.2, length_pd_n=10,
            radius_pd=.2, radius_pd_n=9,
            theta_pd=10, theta_pd_n=5,
           )
q = 0.1
# april 6 2017, rkh added a 2d unit test, assume correct!
qx = q*cos(pi/6.0)
qy = q*sin(pi/6.0)
# Parameters for unit tests
tests = [
    [{}, 0.00005, 1764.926],
    [{}, 'VR', 1.8],
    [{}, 0.001, 1756.76],
    [{}, (qx, qy), 2.36885476192],
]
del qx, qy  # not necessary to delete, but cleaner
