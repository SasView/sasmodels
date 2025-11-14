# octahedron_truncated model
# Note: model title and parameter table are inserted automatically
r"""

This model provides the form factor, $P(q)$, for a general octahedron.
It can be a regular octahedron shape with all edges of the same length.
Or a general shape with different elongations along the three perpendicular two-fold axes.
It includes the possibility to add an adjustable square truncation at all the six vertices.
This model includes the general cuboctahedron shape for the maximum value of truncation.
The form factor expression is obtained by analytical integration over the volume of the shape.
This model is constructed in a similar way as the rectangular prism model.
It contains both the form factor for a fixed orientation and the 1D form factor after orientation average (Gauss-Legendre).

Definition
----------

The general octahedron is defined by its dimensions along its three perpendicular two-fold axis along x, y and z directions.
length_a, length_b and length_c are the distances from the center of the general octahedron to its 6 vertices.

Coordinates of the six vertices are:
    (length_a, 0, 0)
    (-length_a, 0, 0)
    (0, length_b, 0)
    (0, -length_b, 0)
    (0, 0, length_c)
    (0, 0, -length_c)

t is the truncation parameter.
Truncation adds a square facet for each vertice that is perpendicular to a 2-fold axis.
The resulting shape is made of 6 squares and 8 hexagons, that may be non regular depending on the three distances.
A square facet crosses the x, y, z directions at distances equal to t*length_a, t*length_b and t*length_c.

A regular octahedron corresponds to:
    length_a = length_b = length_c
    t = 1

A regular cuboctahedron shape with 6 squares and 8 triangles corresponds to:
    length_a = length_b = length_c
    t = 1/2

The model contains 4 parameters: length_a, the two ratios b2a_ratio and c2a_ratio and t:
    b2a_ratio = length_b/length_a
    c2a_ratio = length_c/length_a
    1/2 < t < 1

For a regular shape:
    b2a_ratio = c2a_ratio = 1

Volume of the general shape including truncation is:
V = (4/3) * length_a * (length_a*b2a_ratio) * (length_a*c2a_ratio)*(1-3(1-t)^3)

The general octahedron is made of eight triangular faces.
The three lengths A_edge, B_edge and C_edge of their edges are equal to:

    A_edge^2 = length_a^2+length_b^2
    B_edge^2 = length_a^2+length_c^2
    C_edge^2 = length_b^2+length_c^2

For a regular shape:
    b2a_ratio = c2a_ratio = 1
    A_edge = B_edge = C_edge = length_a*sqrt(2)
    length_a = length_b = length_c = A_edge/sqrt(2)
    V = (4/3) * length_a^3 = (sqrt(2)/3) * A_edge^3 * (1-3(1-t)^3)

Amplitude of the form factor AP for the fixed orientation of the shape reads:
    AP = 6./(1.-3(1.-t)*(1.-t)*(1.-t))*(AA+BB+CC);

    with:

    AA = 1./(2*(qy*qy-qz*qz)*(qy*qy-qx*qx))*((qy-qx)*sin(qy*(1.-t)-qx*t)+(qy+qx)*sin(qy*(1.-t)+qx*t))+ \
                1./(2*(qz*qz-qx*qx)*(qz*qz-qy*qy))*((qz-qx)*sin(qz*(1.-t)-qx*t)+(qz+qx)*sin(qz*(1.-t)+qx*t));

    BB = 1./(2*(qz*qz-qx*qx)*(qz*qz-qy*qy))*((qz-qy)*sin(qz*(1.-t)-qy*t)+(qz+qy)*sin(qz*(1.-t)+qy*t))+ \
                1./(2*(qx*qx-qy*qy)*(qx*qx-qz*qz))*((qx-qy)*sin(qx*(1.-t)-qy*t)+(qx+qy)*sin(qx*(1.-t)+qy*t));

    CC = 1./(2*(qx*qx-qy*qy)*(qx*qx-qz*qz))*((qx-qz)*sin(qx*(1.-t)-qz*t)+(qx+qz)*sin(qx*(1.-t)+qz*t))+ \
                1./(2*(qy*qy-qz*qz)*(qy*qy-qx*qx))*((qy-qz)*sin(qy*(1.-t)-qz*t)+(qy+qz)*sin(qy*(1.-t)+qz*t));

    and:
            Qx = q * sin_theta * cos_phi;
            Qy = q * sin_theta * sin_phi
            Qz = q * cos_theta
            qx = Qx * length_a
            qy = Qy * length_b
            qz = Qz * length_c

$\theta$ is the angle between the scattering vector and the $z$ axis of the octahedron ($length_c$).
$\phi$ is the angle between the scattering vector (lying in the $xy$ plane) and the $x$ axis ($length_a$).

The normalized form factor in 1D is obtained after averaging over all possible
orientations and this part is implemented in the same way as in the rectangular prism model.

.. math::
    P(q) =  \frac{2}{\pi} \int_0^{\frac{\pi}{2}} \,
    \int_0^{\frac{\pi}{2}} A_P^2(q) \, \sin\theta \, d\theta \, d\phi

And the 1D scattering intensity is calculated as

.. math::
    I(q) = \text{scale} \times V \times (\rho_\text{p} -
    \rho_\text{solvent})^2 \times P(q)

where $V$ is the volume of the truncated octahedron, $\rho_\text{p}$
is the scattering length inside the volume, $\rho_\text{solvent}$
is the scattering length of the solvent, and (if the data are in absolute
units) *scale* represents the volume fraction (which is unitless).


.. figure:: img/octa-truncated.png

    Truncated octahedron shape.

.. figure:: img/octahedrons_intensity_plot.png

    Scattering intensity of a cuboctahedron (t=0.5) and a regular octahedron (t=1) of a = 300 Angstroms.

Validation
----------

Validation of the code is made using numerical checks.
Comparisons with Debye formula calculations were made using DebyeCalcultor library (https://github.com/FrederikLizakJohansen/DebyeCalculator).
Good agreement was found at q < 0.1 1/Angstrom.

References
----------
1. Wei-Ren Chen et al. "Scattering functions of Platonic solids".
   In: Journal of Applied Crystallography - J APPL CRYST 44 (June 2011).
   DOI: 10.1107/S0021889811011691
2. Croset, Bernard, "Form factor of any polyhedron: a general compact
   formula and its singularities" In: J. Appl. Cryst. (2017). 50, 1245–1255
   https://doi.org/10.1107/S1600576717010147
3. Wuttke, J. Numerically stable form factor of any polygon and polyhedron
   J Appl Cryst 54, 580-587 (2021)
   https://doi.org/10.1107/S160057672100171

Authorship and Verification
----------------------------

* **Authors:** Marianne Imperor-Clerc (marianne.imperor@universite-paris-saclay.fr)
             Helen Ibrahim (helenibrahim1@outlook.com)
             Sara Mokhtari (smokhtari@insa-toulouse.fr)

* **Last Modified by:** SM **Date:** 13 November 2025

* **Last Reviewed by:** SM **Date:** November 2025

"""

from numpy import inf

name = "octahedron_truncated"
title = "Truncated Octahedron."
description = """
            AP = 6./(1.-3(1.-t)*(1.-t)*(1.-t))*(AA+BB+CC);
            AA = 1./(2*(qy*qy-qz*qz)*(qy*qy-qx*qx))*((qy-qx)*sin(qy*(1.-t)-qx*t)+(qy+qx)*sin(qy*(1.-t)+qx*t))+
                                1./(2*(qz*qz-qx*qx)*(qz*qz-qy*qy))*((qz-qx)*sin(qz*(1.-t)-qx*t)+(qz+qx)*sin(qz*(1.-t)+qx*t));

            BB = 1./(2*(qz*qz-qx*qx)*(qz*qz-qy*qy))*((qz-qy)*sin(qz*(1.-t)-qy*t)+(qz+qy)*sin(qz*(1.-t)+qy*t))+
                                1./(2*(qx*qx-qy*qy)*(qx*qx-qz*qz))*((qx-qy)*sin(qx*(1.-t)-qy*t)+(qx+qy)*sin(qx*(1.-t)+qy*t));

            CC = 1./(2*(qx*qx-qy*qy)*(qx*qx-qz*qz))*((qx-qz)*sin(qx*(1.-t)-qz*t)+(qx+qz)*sin(qx*(1.-t)+qz*t))+
                                1./(2*(qy*qy-qz*qz)*(qy*qy-qx*qx))*((qy-qz)*sin(qy*(1.-t)-qz*t)+(qy+qz)*sin(qy*(1.-t)+qz*t));

normalisation to 1. of AP at q = 0. Division by a Factor 4/3.
Qx = q * sin_theta * cos_phi;
Qy = q * sin_theta * sin_phi;
Qz = q * cos_theta;
qx = Qx * length_a;
qy = Qy * length_b;
qz = Qz * length_c;
0 < t < 1

"""
category = "shape:parallelepiped"

#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld", "1e-6/Ang^2", 126., [-inf, inf], "sld",
               "Octahedron scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 9.4, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["length_a", "Ang", 400, [0, inf], "volume",
               "half height along a axis"],
              ["b2a_ratio", "", 1, [0, inf], "volume",
               "Ratio b/a"],
              ["c2a_ratio", "", 1, [0, inf], "volume",
               "Ratio c/a"],
              ["t", "", 0.89, [0.5, 1.0], "volume",
               "truncation"],
              ["theta", "degrees", 0, [-360, 360], "orientation",
               "c axis to beam angle"],
              ["phi", "degrees", 0, [-360, 360], "orientation",
               "rotation about beam"],
              ["psi", "degrees", 0, [-360, 360], "orientation",
               "rotation about c axis"],
              ]

source = ["lib/gauss20.c", "octahedron_truncated.c"]
# change to "lib/gauss76.c" or "lib/gauss150.c" to increase the number of integration points
# for the orientational average. Note that it will increase calculation times.
have_Fq = True

tests = [
    [{"background": 0, "scale": 1, "length_a": 100, "t": 1, "sld": 1., "sld_solvent": 0.},
     0.01, 120.57218749910827],
    [{"background": 0, "scale": 1, "length_a": 100, "t": 1, "sld": 1., "sld_solvent": 0.},
     [0.01, 0.1], [120.57218749910827, 0.3741758252985113]],
]
