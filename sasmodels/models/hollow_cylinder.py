r"""
This model provides the form factor, *P(q)*, for a monodisperse hollow right 
angle circular cylinder (tube) where the form factor is normalized by the
volume of the tube

*P(q)* = *scale* \* *<F*\ :sup:`2`\ *>* / *V*\ :sub:`shell` + *background*

where the averaging < > is applied only for the 1D calculation.

The inside and outside of the hollow cylinder are assumed have the same SLD.

Definition
----------

The 1D scattering intensity is calculated in the following way (Guinier, 1955)

.. math::

    \begin{eqnarray}
    P(q)&=&(\text{scale})V_{shell}(\Delta\rho)^2\int_0^{1}\Psi^2[q_z,
    R_{shell}(1-x^2)^{1/2},R_{core}(1-x^2)^{1/2}][\frac{sin(qHx)}{qHx}]^2dx\\
    \Psi[q,y,z]&=&\frac{1}{1-\gamma^2}[\Lambda(qy)-\gamma^2\Lambda(qz)]\\
    \Lambda(a)&=&2J_1(a)/a\\
    \gamma&=&R_{core}/R_{shell}\\
    V_{shell}&=&\pi(R_{shell}^2-R_{core}^2)L\\
    J_1(x)&=&\frac{(sin(x)-x\cdot cos(x))}{x^2}\\
    \end{eqnarray}

where *scale* is a scale factor and *J1* is the 1st order Bessel function.

To provide easy access to the orientation of the core-shell cylinder, we define
the axis of the cylinder using two angles |theta| and |phi|\ . As for the case 
of the cylinder, those angles are defined in Figure 2 of the CylinderModel.

NB: The 2nd virial coefficient of the cylinder is calculated based on the radius
and 2 length values, and used as the effective radius for *S(Q)* when 
*P(Q)* \* *S(Q)* is applied.

In the parameters, the contrast represents SLD :sub:`shell` - SLD :sub:`solvent`
and the *radius* = *R*\ :sub:`shell` while *core_radius* = *R*\ :sub:`core`.

.. image:: img/image074.jpg

*Figure. 1D plot using the default values (w/1000 data point).*

Our model uses the form factor calculations implemented in a c-library provided
by the NIST Center for Neutron Research (Kline, 2006).

.. image:: img/image061.jpg

*Figure. Definition of the angles for the oriented HollowCylinderModel.*

.. image:: img/image062.jpg

*Figure. Examples of the angles for oriented pp against the detector plane.*

REFERENCE

L A Feigin and D I Svergun, *Structure Analysis by Small-Angle X-Ray and
Neutron Scattering*, Plenum Press, New York, (1987)
"""

from numpy import inf

name = "hollow_cylinder"
title = ""
description = """
P(q) = scale*<f*f>/Vol + background, where f is the scattering amplitude.
core_radius = the radius of core
radius = the radius of shell
length = the total length of the cylinder
sld = SLD of the shell
solvent_sld = SLD of the solvent
background = incoherent background
"""
category = "shape:cylinder"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [
              ["radius", "Ang", 30.0, [0, inf], "", "Cylinder radius"],
              ["core_radius", "Ang", 20.0, [0, inf], "", "Hollow core radius"],
              ["length", "Ang", 400.0, [0, inf], "", "Cylinder length"],
              ["sld", "1/Ang^2", 6.3, [-inf, inf], "", "Cylinder sld"],
              ["solvent_sld", "1/Ang^2", 1, [-inf, inf], "", "Solvent sld"],
              ["axis_theta", "[deg]", 90, [-360, 360], "", "Theta angle"],
              ["axis_phi", "[deg]", 0, [-360, 360], "", "Phi angle"],
              ]

source = ["lib/J1.c", "lib/gauss76.c", "hollow_cylinder.c"]

# parameters for demo
demo = dict(scale=1.0,background=0.0,length=400.0,radius=30.0,core_radius=20.0,
            sld=6.3,solvent_sld=1,axis_theta=90,axis_phi=0)

# For testing against the old sasview models, include the converted parameter
# names and the target sasview model name.
oldname = 'HollowCylinderModel'
oldpars = dict(scale='scale',background='background',radius='radius',
               core_radius='core_radius',sld='sldCyl',length='length',
               solvent_sld='sldSolv',axis_phi='axis_phi',axis_theta='axis_theta')