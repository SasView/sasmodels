r"""
Definition
----------

Adds an elliptical rod shaped "protein" through the center of 
core_shell_bicelle_elliptical_belt_rough
The protein sticks out beyond the face of the bicelle, on BOTH faces, 
by distance thick_protein. The whole particle remains centrosymmetric.
The assymetric case, where the protein sticks out further on one face 
than the other is not yet possible in sasview,as that will require a 
completely different I(Q) calculation involving complex algebra.

NOTE it is up to the user to keep the elliptical protein cross section
inside the elliptical core of the bicelle, else non-physical results will
be produced!


.. figure:: img/core_shell_bicelle_elliptical_protein.png

   Schematic cross section of core_shell_bicelle_elliptical_protein. Users
   will have to decide how to distribute "heads" and "tails" between the
   rim, face and core regions in order to estimate appropriate starting
   parameters.

The scattering length densities (sld) are: $\rho_c$, the core, $\rho_f$,
the face, $\rho_r$, the rim, $\rho_p$, the protein,  and $\rho_s$ the solvent.

The form factor for the protein containing bicelle is calculated in cylindrical
coordinates, where
$\alpha$ is the angle between the $Q$ vector and the cylinder axis, and $\psi$
is the angle for the ellipsoidal cross section core, to give:

.. math::

    I(Q,\alpha,\psi) = \frac{\text{scale}}{V_t}
        \cdot F(Q,\alpha, \psi)^2 \cdot \sin(\alpha)
        \cdot\exp\left\{ -\frac{1}{2}Q^2\sigma^2 \right\} + \text{background}

where a numerical integration of $F(Q,\alpha, \psi)^2\sin(\alpha)$ is
carried out over $\alpha$ and $\psi$ for:

.. math::

    F(Q,\alpha,\psi) = &\bigg[
      (\rho_c -\rho_r - \rho_f + \rho_s) V_c
      \frac{2J_1(QR'\sin \alpha)}{QR'\sin\alpha}
      \frac{\sin(Q(L/2)\cos\alpha)}{Q(L/2)\cos\alpha} \\
    &+(\rho_f - \rho_s) V_{c+f}
      \frac{2J_1(QR'\sin\alpha)}{QR'\sin\alpha}
      \frac{\sin(Q(L/2+t_f)\cos\alpha)}{Q(L/2+t_f)\cos\alpha} \\
    &+(\rho_r - \rho_s) V_{c+r}
      \frac{2J_1(Q(R'+t_r)\sin\alpha)}{Q(R'+t_r)\sin\alpha} 
      \frac{\sin(Q(L/2)\cos\alpha)}{Q(L/2)\cos\alpha} \\
    &+\bigg((\rho_p - \rho_s) V_p\frac{\sin(Q(L/2 +t_f +t_p)\cos\alpha)}{Q(L/2 +t_f +t_p)\cos\alpha} \\
    &+(\rho_s - \rho_f) V_d\frac{\sin(Q(L/2 +t_f)\cos\alpha)}{Q(L/2 +t_f)\cos\alpha} \\
    &+(\rho_f - \rho_c) V_e\frac{\sin(Q(L/2)\cos\alpha)}{Q(L/2)\cos\alpha} \bigg)
    \frac{2J_1(QR'_p\sin\alpha)}{QR'_p\sin\alpha}\bigg]

where

.. math::

    R' = \frac{R}{\sqrt{2}}
        \sqrt{(1+X_\text{core}^{2}) + (1-X_\text{core}^{2})\cos(\psi)} \\
    R'_p = \frac{R_p}{\sqrt{2}}
        \sqrt{(1+X_\text{protein}^{2}) + (1-X_\text{protein}^{2})\cos(\psi)}


and $V_t = \pi (R+t_r)(X_\text{core} R+t_r) L + 2 \pi X_\text{core} R^2 t_f+2\pi X_\text{protein}R_p^2 thick_\text{protein}$ is
the total volume of the bicelle plus any protein sticking out from the faces,
$V_c = \pi X_\text{core} R^2 L$ the volume of
the core, $V_{c+f} = \pi X_\text{core} R^2 (L+2 t_f)$ the volume of the core
plus the volume of the faces, $V_{c+r} = \pi (R+t_r)(X_\text{core} R+t_r) L$
the volume of the core plus the rim, $V_d = \pi X_\text{protein} R_p^2 (L+2t_f)$ is the
volume of protein intersecting the bicelle, $V_e = \pi X_\text{protein} R_p^2 L$ is the
volume of protein intersecting the core of the bicelle.

$R$ is the radius of the core,
$X_\text{core}$ is the axial ratio of the core, $L$ the length of the core,
$t_f$ the thickness of the face, $t_r$ the thickness of the rim and $J_1$ the
usual first order bessel function. The bicelle core has radii $R$ and $X_\text{core} R$
so is circular, as for the core_shell_bicelle model, for $X_\text{core}=1$.
Note that you may need to limit the range of $X_\text{core}$, especially if
using the Monte-Carlo algorithm, as setting radius to $R/X_\text{core}$ and
axial ratio to $1/X_\text{core}$ gives an equivalent solution!

$R_p$ is the radius of the protein,
$X_\text{protein}$ is the axial ratio of the protein, which has total length $L+2t_f+2thick_\text{protein}$.
The axes of the elliptical cross section protein cylinder are common with the
bicelle, so making $X_\text{protein} < 1$ whilst $X_\text{core} > 1$ can effectively
rotate the protein by 90 degrees relative to the bicelle.

An approximation for the effects of "Gaussian interfacial roughness" $\sigma$
is included, by multiplying $I(Q)$ by
$\exp\left \{ -\frac{1}{2}Q^2\sigma^2 \right \}$. This applies, in some way, to
all interfaces in the model not just the external ones. It will only work well 
when $\sigma$ is very small compared to the particle size.(Note that for a one
dimensional system convolution of the scattering length density profile with
a Gaussian of standard deviation $\sigma$ does exactly this multiplication.)
Leave $\sigma$ set to zero for the usual sharp interfaces.

The output of the 1D scattering intensity function for randomly oriented
bicelles is then given by integrating over all possible $\alpha$ and $\psi$.

For oriented bicelles the *theta*, *phi* and *psi* orientation parameters
will appear when fitting 2D data, for further details of the calculation
and angular dispersions  see :ref:`orientation` .


.. figure:: img/elliptical_cylinder_angle_definition.png


    Definition of the angles for the oriented core_shell_bicelle_elliptical
    particles.



References
----------

.. [#] Solution structure of bilayer membrane-embedded proton-translocating 
pyrophosphatase revealed via small-angle X-ray scattering”, Orion Shih, Yi-Qi 
Yeh, Kuei-Fen Liao, Kun-Mou Li, Jia-Yin Tsai, Chieh-Chin Li, Yun-Wei Chiang, 
Richard K. Heenan, Yuh-Ju Sun & U-Ser Jeng. Materials Chemistry and Physics, 
308(2023)128253.  https://doi.org/10.1016/j.matchemphys.2023.128253

Source
------

Authorship and Verification
----------------------------

* **Author:** Richard Heenan **Date:** April 24, 2019
* **Last Modified by:Richard Heenan, Sept 2023**
* **Last Reviewed by:**
"""

from numpy import inf, sin, cos, pi

name = "core_shell_bicelle_elliptical_protein"
title = "Elliptical cylinder with a core-shell and protein through the center"
description = """
    core_shell_bicelle_elliptical_protein
    Elliptical cylinder core, optional shell on the two flat faces, and "belt" shell of
    uniform thickness on its rim (in this case NOT extending around the end faces),
    with approximate interfacial roughness, all with an elliptical cylinder of protein
    inserted through the center of the bicelle.
    Please see full documentation for equations and further details.
    Involves a double numerical integral around the ellipsoid diameter
    and the angle of the cylinder axis to Q.
    Compare also the core_shell_bicelle and elliptical_cylinder models.
      """
category = "shape:cylinder"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["radius",         "Ang",       30, [0, inf],    "volume",      "Cylinder core radius r_minor"],
    ["x_core",        "None",       3,  [0, inf],    "volume",      "Axial ratio of core, X = r_major/r_minor"],
    ["thick_rim",     "Ang",         8, [0, inf],    "volume",      "Rim or belt shell thickness"],
    ["thick_face",    "Ang",        14, [0, inf],    "volume",      "Cylinder face thickness"],
    ["length",         "Ang",       50, [0, inf],    "volume",      "Cylinder length"],
    ["sld_core",       "1e-6/Ang^2", 4, [-inf, inf], "sld",         "Cylinder core scattering length density"],
    ["sld_face",       "1e-6/Ang^2", 7, [-inf, inf], "sld",         "Cylinder face scattering length density"],
    ["sld_rim",        "1e-6/Ang^2", 1, [-inf, inf], "sld",         "Cylinder rim scattering length density"],
    ["sld_solvent",    "1e-6/Ang^2", 6, [-inf, inf], "sld",         "Solvent scattering length density"],
    ["sigma",          "Ang",     0,    [0, inf],    "",            "Interfacial roughness"],
    ["radius_protein", "Ang",       18, [0, inf],    "volume",      "Protein radius r_minor"],
    ["x_protein",     "None",     1.5,  [0, inf],    "volume",      "Axial ratio of protein, X = r_major/r_minor"],
    ["thick_protein",  "Ang",       12, [0, inf],    "volume",      "Extension of protein beyond faces"],
    ["sld_protein",    "1e-6/Ang^2", 2, [-inf, inf], "sld",         "Protein scattering length density"],
    ["theta",       "degrees",    90.0, [-360, 360], "orientation", "Cylinder axis to beam angle"],
    ["phi",         "degrees",    0,    [-360, 360], "orientation", "Rotation about beam"],
    ["psi",         "degrees",    0,    [-360, 360], "orientation", "Rotation about cylinder axis"],
    ]

# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sas_Si.c", "lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c",
          "core_shell_bicelle_elliptical_protein.c"]
have_Fq = True
radius_effective_modes = [
    "equivalent cylinder excluded volume", "equivalent volume sphere",
    "outer rim average radius", "outer rim min radius",
    "outer max radius", "half outer thickness", "half diagonal",
    ]

# TODO: No random() for core-shell bicelle elliptical belt rough

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
            radius_protein=18.0,
            x_protein=1.5,
            thick_protein=12.0,
            sld_protein=2.0,
            theta=90,
            phi=0,
            psi=0,
            sigma=0)

q = 0.1
# april 6 2017, rkh added a 2d unit test, NOT READY YET pull #890 branch assume correct!
qx = q*cos(pi/6.0)
qy = q*sin(pi/6.0)

tests = [
    #[{'radius': 30.0, 'x_core': 3.0, 'thick_rim':8.0, 'thick_face':14.0, 'length':50.0}, 'ER', 1],
    #[{'radius': 30.0, 'x_core': 3.0, 'thick_rim':8.0, 'thick_face':14.0, 'length':50.0}, 'VR', 1],
#
     # this test as per the original core_shell_bicelle_elliptical_belt_rough
    [{'radius': 30.0, 'x_core': 3.0, 'thick_rim': 8.0, 'thick_face': 14.0,
      'length': 50.0, 'sld_core': 4.0, 'sld_face': 7.0, 'sld_rim': 1.0,
      'radius_protein':0.0, 'x_protein':1.5, 'thick_protein':0.0, 'sld_protein':2.0,
      'sld_solvent': 6.0, 'background': 0.0},
     0.015, 189.328],
    #[{'theta':80., 'phi':10.}, (qx, qy), 7.88866563001 ],
]

del qx, qy  # not necessary to delete, but cleaner
