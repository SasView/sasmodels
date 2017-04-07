r"""
Definition
----------

This model provides the form factor for an elliptical cylinder with a
core-shell scattering length density profile. Thus this is a variation
of the core-shell bicelle model, but with an elliptical cylinder for the core.
Outer shells on the rims and flat ends may be of different thicknesses and 
scattering length densities. The form factor is normalized by the total particle volume.


.. figure:: img/core_shell_bicelle_geometry.png

    Schematic cross-section of bicelle. Note however that the model here
    calculates for rectangular, not curved, rims as shown below.

.. figure:: img/core_shell_bicelle_parameters.png

   Cross section of model used here. Users will have 
   to decide how to distribute "heads" and "tails" between the rim, face 
   and core regions in order to estimate appropriate starting parameters.

Given the scattering length densities (sld) $\rho_c$, the core sld, $\rho_f$,
the face sld, $\rho_r$, the rim sld and $\rho_s$ the solvent sld, the
scattering length density variation along the bicelle axis is:

.. math::

    \rho(r) = 
      \begin{cases} 
      &\rho_c \text{ for } 0 \lt r \lt R; -L \lt z\lt L \\[1.5ex]
      &\rho_f \text{ for } 0 \lt r \lt R; -(L+2t) \lt z\lt -L;
      L \lt z\lt (L+2t) \\[1.5ex]
      &\rho_r\text{ for } 0 \lt r \lt R; -(L+2t) \lt z\lt -L; L \lt z\lt (L+2t)
      \end{cases}

The form factor for the bicelle is calculated in cylindrical coordinates, where
$\alpha$ is the angle between the $Q$ vector and the cylinder axis, and $\psi$ is the angle
for the ellipsoidal cross section core, to give:

.. math::

    I(Q,\alpha,\psi) = \frac{\text{scale}}{V_t} \cdot
        F(Q,\alpha, \psi)^2.sin(\alpha) + \text{background}

where a numerical integration of $F(Q,\alpha, \psi)^2.sin(\alpha)$ is carried out over \alpha and \psi for:

.. math::

        \begin{align}    
    F(Q,\alpha,\psi) = &\bigg[ 
    (\rho_c - \rho_f) V_c \frac{2J_1(QR'sin \alpha)}{QR'sin\alpha}\frac{sin(QLcos\alpha/2)}{Q(L/2)cos\alpha} \\
    &+(\rho_f - \rho_r) V_{c+f} \frac{2J_1(QR'sin\alpha)}{QR'sin\alpha}\frac{sin(Q(L/2+t_f)cos\alpha)}{Q(L/2+t_f)cos\alpha} \\
    &+(\rho_r - \rho_s) V_t \frac{2J_1(Q(R'+t_r)sin\alpha)}{Q(R'+t_r)sin\alpha}\frac{sin(Q(L/2+t_f)cos\alpha)}{Q(L/2+t_f)cos\alpha}
    \bigg]
    \end{align} 

where

.. math::

    R'=\frac{R}{\sqrt{2}}\sqrt{(1+X_{core}^{2}) + (1-X_{core}^{2})cos(\psi)}
    
    
and $V_t = \pi.(R+t_r)(Xcore.R+t_r)^2.(L+2.t_f)$ is the total volume of the bicelle, 
$V_c = \pi.Xcore.R^2.L$ the volume of the core, $V_{c+f} = \pi.Xcore.R^2.(L+2.t_f)$ 
the volume of the core plus the volume of the faces, $R$ is the radius
of the core, $Xcore$ is the axial ratio of the core, $L$ the length of the core, 
$t_f$ the thickness of the face, $t_r$ the thickness of the rim and $J_1$ the usual 
first order bessel function. The core has radii $R$ and $Xcore.R$ so is circular, 
as for the core_shell_bicelle model, for $Xcore$ =1. Note that you may need to 
limit the range of $Xcore$, especially if using the Monte-Carlo algorithm, as 
setting radius to $R/Xcore$ and axial ratio to $1/Xcore$ gives an equivalent solution!

The output of the 1D scattering intensity function for randomly oriented
bicelles is then given by integrating over all possible $\alpha$ and $\psi$.

For oriented bicellles the *theta*, *phi* and *psi* orientation parameters only appear when fitting 2D data, 
see the :ref:`elliptical-cylinder` model for further information.


.. figure:: img/elliptical_cylinder_angle_definition.png

    Definition of the angles for the oriented core_shell_bicelle_elliptical particles.   



References
----------

.. [#]

Authorship and Verification
----------------------------

* **Author:** Richard Heenan **Date:** December 14, 2016
* **Last Modified by:**  Richard Heenan **Date:** December 14, 2016
* **Last Reviewed by:**  Richard Heenan BEWARE 2d data yet to be checked **Date:** December 14, 2016
"""

from numpy import inf, sin, cos, pi

name = "core_shell_bicelle_elliptical"
title = "Elliptical cylinder with a core-shell scattering length density profile.."
description = """
    core_shell_bicelle_elliptical
    Elliptical cylinder core, optional shell on the two flat faces, and shell of 
    uniform thickness on its rim (extending around the end faces).    
    Please see full documentation for equations and further details.
    Involves a double numerical integral around the ellipsoid diameter
    and the angle of the cylinder axis to Q.
    Compare also the core_shell_bicelle and elliptical_cylinder models.
      """
category = "shape:cylinder"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["radius",         "Ang",       30, [0, inf],    "volume",      "Cylinder core radius"],
    ["x_core",        "None",       3,  [0, inf],    "volume",      "axial ratio of core, X = r_polar/r_equatorial"],
    ["thick_rim",  "Ang",        8, [0, inf],    "volume",      "Rim shell thickness"],
    ["thick_face", "Ang",       14, [0, inf],    "volume",      "Cylinder face thickness"],
    ["length",         "Ang",      50, [0, inf],    "volume",      "Cylinder length"],
    ["sld_core",       "1e-6/Ang^2", 4, [-inf, inf], "sld",         "Cylinder core scattering length density"],
    ["sld_face",       "1e-6/Ang^2", 7, [-inf, inf], "sld",         "Cylinder face scattering length density"],
    ["sld_rim",        "1e-6/Ang^2", 1, [-inf, inf], "sld",         "Cylinder rim scattering length density"],
    ["sld_solvent",    "1e-6/Ang^2", 6, [-inf, inf], "sld",         "Solvent scattering length density"],
    ["theta",          "degrees",   90, [-360, 360], "orientation", "In plane angle"],
    ["phi",            "degrees",    0, [-360, 360], "orientation", "Out of plane angle"],
    ["psi",            "degrees",    0, [-360, 360], "orientation", "Major axis angle relative to Q"],
    ]

# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sas_Si.c", "lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c",
          "core_shell_bicelle_elliptical.c"]

demo = dict(scale=1, background=0,
            radius=30.0,
            x_core=3.0,
            thick_rim=8.0,
            thick_face=14.0,
            length=50.0,
            sld_core=4.0,
            sld_face=7.0,
            sld_rim=1.0,
            sld_solvent=6.0,
            theta=90,
            phi=0,
            psi=0)

q = 0.1
# april 6 2017, rkh added a 2d unit test, NOT READY YET pull #890 branch assume correct!
qx = q*cos(pi/6.0)
qy = q*sin(pi/6.0)

tests = [
    [{'radius': 30.0, 'x_core': 3.0, 'thick_rim':8.0, 'thick_face':14.0, 'length':50.0}, 'ER', 1],
    [{'radius': 30.0, 'x_core': 3.0, 'thick_rim':8.0, 'thick_face':14.0, 'length':50.0}, 'VR', 1],

    [{'radius': 30.0, 'x_core': 3.0, 'thick_rim':8.0, 'thick_face':14.0, 'length':50.0,
    'sld_core':4.0, 'sld_face':7.0, 'sld_rim':1.0, 'sld_solvent':6.0, 'background':0.0},
    0.015, 286.540286],
#    [{'theta':80., 'phi':10.}, (qx, qy), 7.88866563001 ],
        ]

del qx, qy  # not necessary to delete, but cleaner
