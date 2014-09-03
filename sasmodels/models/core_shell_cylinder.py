# core shell cylinder model
# Note: model title and parameter table are inserted automatically
r"""
The form factor is normalized by the particle volume.

Definition
----------

The output of the 2D scattering intensity function for oriented core-shell
cylinders is given by (Kline, 2006)

.. math::

    P(Q,\alpha) = {\text{scale} \over V_s} F^2(Q) + \text{background}

where

.. math::

    F(Q) = &\ (\rho_c - \rho_s) V_c
           {\sin \left( Q \tfrac12 L\cos\alpha \right)
               \over Q \tfrac12 L\cos\alpha }
           {2 J_1 \left( QR\sin\alpha \right)
               \over QR\sin\alpha } \\
         &\ + (\rho_s - \rho_\text{solv}) V_s
           {\sin \left( Q \left(\tfrac12 L+T\right) \cos\alpha \right)
               \over Q \left(\tfrac12 L +T \right) \cos\alpha }
           { 2 J_1 \left( Q(R+T)\sin\alpha \right)
               \over Q(R+T)\sin\alpha }

and

.. math::

    V_s = \pi (R + T)^2 (L + 2T)

and $\alpha$ is the angle between the axis of the cylinder and $\vec q$,
$V_s$ is the volume of the outer shell (i.e. the total volume, including
the shell), $V_c$ is the volume of the core, $L$ is the length of the core,
$R$ is the radius of the core, $T$ is the thickness of the shell, $\rho_c$
is the scattering length density of the core, $\rho_s$ is the scattering
length density of the shell, $\rho_\text{solv}$ is the scattering length
density of the solvent, and *background* is the background level.  The outer
radius of the shell is given by $R+T$ and the total length of the outer
shell is given by $L+2T$. $J1$ is the first order Bessel function.

.. _core-shell-cylinder-geometry:

.. figure:: img/core_shell_cylinder_geometry.jpg

    Core shell cylinder schematic.

To provide easy access to the orientation of the core-shell cylinder, we
define the axis of the cylinder using two angles $\theta$ and $\phi$. As
for the case of the cylinder, those angles are defined in
:num:`figure #cylinder-orientation`.

NB: The 2nd virial coefficient of the cylinder is calculated based on
the radius and 2 length values, and used as the effective radius for
$S(Q)$ when $P(Q) \cdot S(Q)$ is applied.

The $\theta$ and $\phi$ parameters are not used for the 1D output. Our
implementation of the scattering kernel and the 1D scattering intensity
use the c-library from NIST.

Validation
----------

Validation of our code was done by comparing the output of the 1D model to
the output of the software provided by the NIST (Kline, 2006).
:num:`Figure #core-shell-cylinder-1d` shows a comparison
of the 1D output of our model and the output of the NIST software.

.. _core-shell-cylinder-1d:

.. figure:: img/core_shell_cylinder_1d.jpg

    Comparison of the SasView scattering intensity for a core-shell cylinder
    with the output of the NIST SANS analysis software. The parameters were
    set to: *scale* = 1.0 |Ang|, *radius* = 20 |Ang|, *thickness* = 10 |Ang|,
    *length* =400 |Ang|, *core_sld* =1e-6 |Ang^-2|, *shell_sld* = 4e-6 |Ang^-2|,
    *solvent_sld* = 1e-6 |Ang^-2|, and *background* = 0.01 |cm^-1|.

Averaging over a distribution of orientation is done by evaluating the
equation above. Since we have no other software to compare the
implementation of the intensity for fully oriented cylinders, we can
compare the result of averaging our 2D output using a uniform
distribution $p(\theta,\phi) = 1.0$.
:num:`Figure #core-shell-cylinder-2d` shows the result
of such a cross-check.

.. _core-shell-cylinder-2d:

.. figure:: img/core_shell_cylinder_2d.jpg

    Comparison of the intensity for uniformly distributed core-shell
    cylinders calculated from our 2D model and the intensity from the
    NIST SANS analysis software. The parameters used were: *scale* = 1.0,
    *radius* = 20 |Ang|, *thickness* = 10 |Ang|, *length* = 400 |Ang|,
    *core_sld* = 1e-6 |Ang^-2|, *shell_sld* = 4e-6 |Ang^-2|,
    *solvent_sld* = 1e-6 |Ang^-2|, and *background* = 0.0 |cm^-1|.

2013/11/26 - Description reviewed by Heenan, R.
"""

from numpy import pi, inf

name = "core_shell_cylinder"
title = "Right circular cylinder with a core-shell scattering length density profile."
description = """
P(q,alpha)= scale/Vs*f(q)^(2) + background,
      where: f(q)= 2(core_sld - solvant_sld)
		* Vc*sin[qLcos(alpha/2)]
		/[qLcos(alpha/2)]*J1(qRsin(alpha))
		/[qRsin(alpha)]+2(shell_sld-solvent_sld)
		*Vs*sin[q(L+T)cos(alpha/2)][[q(L+T)
		*cos(alpha/2)]*J1(q(R+T)sin(alpha))
		/q(R+T)sin(alpha)]

	alpha:is the angle between the axis of
		the cylinder and the q-vector
	Vs: the volume of the outer shell
	Vc: the volume of the core
	L: the length of the core
    	shell_sld: the scattering length density of the shell
	solvent_sld: the scattering length density of the solvent
	background: the background
	T: the thickness
    	R+T: is the outer radius
 	L+2T: The total length of the outershell
	J1: the first order Bessel function
 	theta: axis_theta of the cylinder
 	phi: the axis_phi of the cylinder
"""

parameters = [
#   [ "name", "units", default, [lower, upper], "type",
#     "description" ],
    [ "core_sld", "1e-6/Ang^2", 4, [-inf,inf], "",
      "Cylinder core scattering length density" ],
    [ "shell_sld", "1e-6/Ang^2", 4, [-inf,inf], "",
      "Cylinder shell scattering length density" ],
    [ "solvent_sld", "1e-6/Ang^2", 1, [-inf,inf], "",
      "Solvent scattering length density" ],
    [ "radius", "Ang",  20, [0, inf], "volume",
      "Cylinder core radius" ],
    [ "thickness", "Ang",  20, [0, inf], "volume",
      "Cylinder shell thickness" ],
    [ "length", "Ang",  400, [0, inf], "volume",
      "Cylinder length" ],
    [ "theta", "degrees", 60, [-inf, inf], "orientation",
      "In plane angle" ],
    [ "phi", "degrees", 60, [-inf, inf], "orientation",
      "Out of plane angle" ],
    ]

source = [ "lib/J1.c", "lib/gauss76.c", "core_shell_cylinder.c"]

def ER(radius, thickness, length):
    radius = radius + thickness
    length = length + 2*thickness
    ddd = 0.75*radius*(2*radius*length + (length+radius)*(length+pi*radius))
    return 0.5 * (ddd)**(1./3.)

def VR(radius, thickness, length):
    whole = pi * (radius+thickness)**2 * (length+2*thickness)
    core = pi * radius**2 * length
    return whole, whole-core

