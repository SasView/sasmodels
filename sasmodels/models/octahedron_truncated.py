# octahedron_truncated model
# Note: model title and parameter table are inserted automatically
r"""

This model provides the form factor P(q) for a general octahedron.
It can be a regular octahedron shape with all edges of the same length.
Or a general shape with different elongations along the three perpendicular two-fold axes.
It includes the possibility to add an adjustable square truncation at each of the six vertices.
This model includes the general cuboctahedron shape for the maximum value of truncation.
The form factor expression is obtained by analytical integration over the volume of the shape.
This model is constructed in a similar way as the rectangular prism model.
It contains both the form factor for a reference orientation and the 1D form factor after orientation average (Gauss-Legendre).

Definition
----------

The general octahedron is defined by its dimensions along its three perpendicular two-fold axes along x, y and z directions.
length_a, length_b and length_c are the distances from the center of the general octahedron to its 6 vertices.

Coordinates of the six vertices are:
    (length_a, 0, 0),
    (-length_a, 0, 0),
    (0, length_b, 0),
    (0, -length_b, 0),
    (0, 0, length_c),
    (0, 0, -length_c)

t is the truncation parameter.
Truncation adds a square facet for each vertex that is perpendicular to a 2-fold axis.
The resulting shape consists of six squares and eight hexagons, which may be irregular depending on the three dimensions
A square facet crosses the x, y, z directions at distances equal to t length_a, t length_b and t length_c.

A regular octahedron corresponds to:

.. math::

    length_a = length_b = length_c, \quad t = 1

A regular cuboctahedron shape with 6 squares and 8 triangles corresponds to:

.. math:: 

    length_a = length_b = length_c, \quad t = \frac{1}{2}

The model contains 4 parameters: length_a, the two ratios b2a_ratio and c2a_ratio and t:

.. math::

    b2a_{\text{ratio}} = \frac{length_b}{length_a}, \quad
    c2a_{\text{ratio}} = \frac{length_c}{length_a}, \quad
    \frac{1}{2} < t < 1

For a regular shape:

.. math::

    b2a_{\text{ratio}} = c2a_{\text{ratio}} = 1

Volume of the general shape including truncation is given by:

.. math::

    V = \frac{4}{3}\, \text{length\_a}^{3}\, \text{b2a\_ratio}\, \text{c2a\_ratio}\,\bigl(1 - 3(1 - t)^{3}\bigr)

The general octahedron is made of eight triangular faces. The three edge lengths
are:

.. math::

    A_{\text{edge}}^{2} = \text{length\_a}^{2} + \text{length\_b}^{2},\qquad
    B_{\text{edge}}^{2} = \text{length\_a}^{2} + \text{length\_c}^{2},\qquad
    C_{\text{edge}}^{2} = \text{length\_b}^{2} + \text{length\_c}^{2}

For a regular shape (no elongation):

.. math::

    b2a_{\text{ratio}} = c2a_{\text{ratio}} = 1,\qquad
    A_{\text{edge}} = B_{\text{edge}} = C_{\text{edge}} = \text{length\_a} \sqrt{2},\qquad
    length_{\text{a}} = length_{\text{b}} = length_{\text{c}} = A_{\text{edge}} / \sqrt{2}

.. math::

    V = \frac{4}{3} \, \text{length\_a}^{3} \, \bigl(1 - 3 (1 - t)^3 \bigr)

The reference orientation of the shape is: a along x, b along y and c along z.
Amplitude of the form factor AP for the reference orientation of the shape reads 

.. math::

    AP(q,\theta,\phi) = \frac{6}{1 - 3(1 - t)^3}\,(AA + BB + CC)

.. math::

    AA = \frac{1}{2\,(q_y^2 - q_z^2)\,(q_y^2 - q_x^2)}\Big[(q_y - q_x)\sin\big(q_y(1 - t) - q_x t\big)
    + (q_y + q_x)\sin\big(q_y(1 - t) + q_x t\big)\Big]
    + \frac{1}{2\,(q_z^2 - q_x^2)\,(q_z^2 - q_y^2)}\Big[(q_z - q_x)\sin\big(q_z(1 - t) - q_x t\big)
    + (q_z + q_x)\sin\big(q_z(1 - t) + q_x t\big)\Big]

.. math::

    BB = \frac{1}{2\,(q_z^2 - q_x^2)\,(q_z^2 - q_y^2)}\Big[(q_z - q_y)\sin\big(q_z(1 - t) - q_y t\big)
    + (q_z + q_y)\sin\big(q_z(1 - t) + q_y t\big)\Big]
    + \frac{1}{2\,(q_x^2 - q_y^2)\,(q_x^2 - q_z^2)}\Big[(q_x - q_y)\sin\big(q_x(1 - t) - q_y t\big)
    + (q_x + q_y)\sin\big(q_x(1 - t) + q_y t\big)\Big]

.. math::

    CC = \frac{1}{2\,(q_x^2 - q_y^2)\,(q_x^2 - q_z^2)}\Big[(q_x - q_z)\sin\big(q_x(1 - t) - q_z t\big)
    + (q_x + q_z)\sin\big(q_x(1 - t) + q_z t\big)\Big]
    + \frac{1}{2\,(q_y^2 - q_z^2)\,(q_y^2 - q_x^2)}\Big[(q_y - q_z)\sin\big(q_y(1 - t) - q_z t\big)
    + (q_y + q_z)\sin\big(q_y(1 - t) + q_z t\big)\Big]

Capital Qx Qy Qz are the three components in [A-1] of the scattering vector.
qx qy qz are rescaled components (no unit) for computing AA, BB and CC terms.

.. math::

    Q_x = q\,\sin\theta\,\cos\phi, \qquad
    Q_y = q\,\sin\theta\,\sin\phi, \qquad
    Q_z = q\,\cos\theta

.. math::

    q_x = Q_x \, \text{length\_a},\qquad
    q_y = Q_y \, \text{length\_b},\qquad
    q_z = Q_z \, \text{length\_c}


θ is the angle between the scattering vector and the z axis.
ϕ is the rotation angle in the xy plane.

The octahedron is in its reference orientation, with the c-axis aligned along z and the a-axis aligned along x.

The 1D form factor P(q) corresponds to the orientation average with all the possible orientations having the same probability.
Instead of rotating the shape through all the possible orientations in the integral,
it is equivalent to integrate the 3D scattering vector over a sphere of radius q with the shape in its reference orientation.

.. math::

    P(q) =  \frac{2}{\pi} \int_0^{\frac{\pi}{2}} \,
    \int_0^{\frac{\pi}{2}} A_P^2(q) \, \sin\theta \, d\theta \, d\phi

And the 1D scattering intensity is calculated as

.. math::

    I(q) = \text{scale} \times V \times (\rho_\text{p} -
    \rho_\text{solvent})^2 \times P(q)

where V is the volume of the truncated octahedron, ρ
is the scattering length inside the volume, ρ *solvent*
is the scattering length of the solvent, and (if the data are in absolute
units) *scale* represents the volume fraction (which is unitless).


.. figure:: img/octa-truncated.png

    Truncated octahedron shape for different truncation.

.. figure:: img/octahedrons_intensity_plot.png

    Scattering intensity of a cuboctahedron (t=0.5) and a regular octahedron (t=1) of a = 300 Angstroms.

Validation
----------

Validation of the code is made using numerical checks.
Comparisons with Debye formula calculations were made using DebyeCalculator library (https://github.com/FrederikLizakJohansen/DebyeCalculator).
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

* **Authors:** Marianne Imperor-Clerc (marianne.imperor@cnrs.fr)
             Helen Ibrahim (helenibrahim1@outlook.com)
             Sara Mokhtari (smokhtari@insa-toulouse.fr)

* **Last Modified by:** MIC **Date:** 11 December 2025

* **Last Reviewed by:** SM **Date:** 10 December 2025

"""

from numpy import inf

name = "octahedron_truncated"
title = "Truncated Octahedron."
description = """
The amplitude AP is defined as follows.

AP = 6./(1.-3*(1.-t)*(1.-t)*(1.-t))*(AA+BB+CC)

AA = 1./(2*(qy*qy-qz*qz)*(qy*qy-qx*qx)) * ((qy-qx)*sin(qy*(1.-t)-qx*t) + (qy+qx)*sin(qy*(1.-t)+qx*t)) + 1./(2*(qz*qz-qx*qx)*(qz*qz-qy*qy)) * ((qz-qx)*sin(qz*(1.-t)-qx*t) + (qz+qx)*sin(qz*(1.-t)+qx*t))

BB = 1./(2*(qz*qz-qx*qx)*(qz*qz-qy*qy)) * ((qz-qy)*sin(qz*(1.-t)-qy*t) + (qz+qy)*sin(qz*(1.-t)+qy*t)) + 1./(2*(qx*qx-qy*qy)*(qx*qx-qz*qz)) * ((qx-qy)*sin(qx*(1.-t)-qy*t) + (qx+qy)*sin(qx*(1.-t)+qy*t))

CC = 1./(2*(qx*qx-qy*qy)*(qx*qx-qz*qz)) * ((qx-qz)*sin(qx*(1.-t)-qz*t) + (qx+qz)*sin(qx*(1.-t)+qz*t)) + 1./(2*(qy*qy-qz*qz)*(qy*qy-qx*qx)) * ((qy-qz)*sin(qy*(1.-t)-qz*t) + (qy+qz)*sin(qy*(1.-t)+qz*t))

With capital QX QY QZ are the three components in [A-1] of the scattering vector
and qx qy qz are the rescaled components (no unit) for computing AP term.

Qx = q * sin_theta * cos_phi
Qy = q * sin_theta * sin_phi
Qz = q * cos_theta
qx = Qx * length_a
qy = Qy * length_b
qz = Qz * length_c

Reference orientation is with a along x axis, b along y axis and c along z axis

Valid truncation parameter range: 0.5 < t < 1.

t=1 is for octahedron
t=0.5 is for cuboctahedron
"""

category = "shape:polyhedron"

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

# Fq() function is used in the .c code
have_Fq = True

tests = [
    [{"background": 0, "scale": 1, "length_a": 100, "t": 1, "sld": 1., "sld_solvent": 0.},
     0.01, 120.57219],
    [{"background": 0, "scale": 1, "length_a": 100, "t": 1, "sld": 1., "sld_solvent": 0.},
     [0.01, 0.1], [120.57219, 0.37418]],
]
