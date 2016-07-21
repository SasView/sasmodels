r"""
This model provides the form factor, $P(q)$, for a monodisperse hollow right
angle circular cylinder (tube) where the form factor is normalized by the
volume of the tube

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
            \left[q_z, R_\text{shell}(1-x^2)^{1/2},
                       R_\text{core}(1-x^2)^{1/2}\right]
            \left[\frac{\sin(qHx)}{qHx}\right]^2 dx \\
    \Psi[q,y,z]    &= \frac{1}{1-\gamma^2}
            \left[ \Lambda(qy) - \gamma^2\Lambda(qz) \right] \\
    \Lambda(a)     &= 2 J_1(a) / a \\
    \gamma         &= R_\text{core} / R_\text{shell} \\
    V_\text{shell} &= \pi \left(R_\text{shell}^2 - R_\text{core}^2 \right)L \\
    J_1(x)         &= (\sin(x)-x\cdot \cos(x)) / x^2

where *scale* is a scale factor, $H = L/2$ and $J_1$ is the 1st order
Bessel function.

**NB**: The 2nd virial coefficient of the cylinder is calculated
based on the radius and 2 length values, and used as the effective radius
for $S(q)$ when $P(q) \cdot S(q)$ is applied.

In the parameters, the contrast represents SLD :sub:`shell` - SLD :sub:`solvent`
and the *radius* is $R_\text{shell}$ while *core_radius* is $R_\text{core}$.

To provide easy access to the orientation of the core-shell cylinder, we define
the axis of the cylinder using two angles $\theta$ and $\phi$
(see :ref:`cylinder model <cylinder-angle-definition>`).

References
----------

L A Feigin and D I Svergun, *Structure Analysis by Small-Angle X-Ray and
Neutron Scattering*, Plenum Press, New York, (1987)
"""

from numpy import pi, inf

name = "hollow_cylinder"
title = ""
description = """
P(q) = scale*<f*f>/Vol + background, where f is the scattering amplitude.
core_radius = the radius of core
radius = the radius of shell
length = the total length of the cylinder
sld = SLD of the shell
sld_solvent = SLD of the solvent
background = incoherent background
"""
category = "shape:cylinder"
# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["radius",      "Ang",     30.0, [0, inf],    "volume",      "Cylinder radius"],
    ["core_radius", "Ang",     20.0, [0, inf],    "volume",      "Hollow core radius"],
    ["length",      "Ang",    400.0, [0, inf],    "volume",      "Cylinder length"],
    ["sld",         "1/Ang^2",  6.3, [-inf, inf], "sld",         "Cylinder sld"],
    ["sld_solvent", "1/Ang^2",  1,   [-inf, inf], "sld",         "Solvent sld"],
    ["theta",       "degrees", 90,   [-360, 360], "orientation", "Theta angle"],
    ["phi",         "degrees",  0,   [-360, 360], "orientation", "Phi angle"],
    ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/polevl.c","lib/sas_J1.c", "lib/gauss76.c", "hollow_cylinder.c"]

# pylint: disable=W0613
def ER(radius, core_radius, length):
    """
    :param radius:      Cylinder radius
    :param core_radius: Hollow core radius, UNUSED
    :param length:      Cylinder length
    :return:            Effective radius
    """
    if radius == 0 or length == 0:
        return 0.0
    len1 = radius
    len2 = length/2.0
    term1 = len1*len1*2.0*len2/2.0
    term2 = 1.0 + (len2/len1)*(1.0 + 1/len2/2.0)*(1.0 + pi*len1/len2/2.0)
    ddd = 3.0*term1*term2
    diam = pow(ddd, (1.0/3.0))
    return diam

def VR(radius, core_radius, length):
    """
    :param radius:      Cylinder radius
    :param core_radius: Hollow core radius
    :param length:      Cylinder length
    :return:            Volf ratio for P(q)*S(q)
    """
    vol_core = pi*core_radius*core_radius*length
    vol_total = pi*radius*radius*length
    vol_shell = vol_total - vol_core
    return vol_shell, vol_total

# parameters for demo
demo = dict(scale=1.0, background=0.0, length=400.0, radius=30.0,
            core_radius=20.0, sld=6.3, sld_solvent=1, theta=90, phi=0,
            radius_pd=.2, radius_pd_n=9,
            length_pd=.2, length_pd_n=10,
            core_radius_pd=.2, core_radius_pd_n=9,
            theta_pd=10, theta_pd_n=5,
           )

# Parameters for unit tests
tests = [
    [{"radius": 30.0}, 0.00005, 1764.926],
    [{}, 'VR', 1.8],
    [{}, 0.001, 1756.76]
    ]
