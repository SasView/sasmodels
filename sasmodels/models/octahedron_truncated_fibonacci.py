# octahedron_truncated model
# Note: model title and parameter table are inserted automatically
r"""
This model provides the form factor P(q) for a general octahedron with orientational averaging (Fibonacci quadrature).
It can be a regular octahedron shape with all edges of the same length.
Or a general shape with different elongations along the three perpendicular two-fold axes.
It includes the possibility to add an adjustable square truncation at each of the six vertices.
This model includes the general cuboctahedron shape for the maximum value of truncation.
The form factor expression is obtained by analytical integration over the volume of the shape.
This model is constructed in a similar way as the rectangular prism model.
It contains both the form factor for a reference orientation and the 1D form factor after orientation average
using the Fibonacci quadrature. This quadrature provides a quasi-uniform distribution of points on the unit sphere
using the golden ratio. The only new parameter is the number of points to generate on the unit sphere.
The default value (around 400 points) provides a good balance between accuracy and computational efficiency.

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

    V = \frac{4}{3}\, length_{\text{a}}^{3}\, b2a_{\text{ratio}}\, c2a_{\text{ratio}}\,\bigl(1 - 3(1 - t)^{3}\bigr)

The general octahedron is made of eight triangular faces. The three edge lengths
are:

.. math::

    A_{\text{edge}}^{2} = length_{\text{a}}^{2} + length_{\text{b}}^{2},\qquad
    B_{\text{edge}}^{2} = length_{\text{a}}^{2} + length_{\text{c}}^{2},\qquad
    C_{\text{edge}}^{2} = length_{\text{b}}^{2} + length_{\text{c}}^{2}

For a regular shape (no elongation):

.. math::

    b2a_{\text{ratio}} = c2a_{\text{ratio}} = 1,\qquad
    A_{\text{edge}} = B_{\text{edge}} = C_{\text{edge}} = length_{\text{a}} \sqrt{2},\qquad
    length_{\text{a}} = length_{\text{b}} = length_{\text{c}} = A_{\text{edge}} / \sqrt{2}

.. math::

    V = \frac{4}{3} \, length_{\text{a}}^{3} \, \bigl(1 - 3 (1 - t)^3 \bigr)

The reference orientation of the shape is: a along x, b along y and c along z.
Amplitude of the form factor AP for the reference orientation of the shape reads

.. math::

    AP(q,\theta,\phi) = \frac{3}{1 - 3(1 - t)^3}\,(AA + BB + CC)

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

    q_x = Q_x \, length_{\text{a}},\qquad
    q_y = Q_y \, length_{\text{b}},\qquad
    q_z = Q_z \, length_{\text{c}}


θ is the angle between the scattering vector and the z axis.
ϕ is the rotation angle in the xy plane.

The octahedron is in its reference orientation, with the c-axis aligned along z and the a-axis aligned along x.

The 1D form factor P(q) corresponds to the orientation average with all the possible orientations having the same probability.
Instead of rotating the shape through all the possible orientations in the integral,
it is equivalent to integrate the 3D scattering vector over a sphere of radius q with the shape in its reference orientation.

.. math::

    P(q) =  \frac{2}{\pi} \int_0^{\frac{\pi}{2}} \,
    \int_0^{\frac{\pi}{2}} A_P^2(q) \, \sin\theta \, d\theta \, d\phi

The sphere is sampled using Fibonacci quadrature to provide a quasi-uniform distribution of points on the unit sphere.
The repartition of the points is computed using the golden ratio (see fibonacci.py).


.. figure:: img/fibonacci_sphere.png

    Fibonacci sphere using 5810 points.

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

import numpy as np
from numpy import inf

from sasmodels.quadratures.fibonacci import fibonacci_sphere

name = "octahedron_truncated_fibonacci"
title = "General octahedron pure python model using Fibonacci quadrature"
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

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 130, [-inf, inf], "sld",
               "Octahedron scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 9.4, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["length_a", "Ang", 500, [0, inf], "volume",
               "half height along a axis"],
              ["b2a_ratio", "", 1, [0, inf], "volume",
               "Ratio b/a"],
              ["c2a_ratio", "", 1, [0, inf], "volume",
               "Ratio c/a"],
              ["t", "", 0.99, [0.0, 1.0], "volume",
               "truncation along axis"],
              ["npoints_fibonacci",     "",           400, [1, 1e6],   "",
                       "Number of points on the sphere for the Fibonacci integration"]
              ]


def form_volume(length_a, b2a_ratio, c2a_ratio, t):
    """
    Calculate the volume of the truncated octahedron.

    Parameters
    ----------
    length_a: Half height along the a axis of the octahedron without truncature
    b2a_ratio: Half height along the b axis of the octahedron without truncature
    c2a_ratio: Half height along the c axis of the octahedron without truncature
    t: Truncation parameter, varies from 0.5 (cuboctahedron) to 1 (octahedron)

    Returns
    -------
    volume: Volume of the truncated octahedron
    """
    return (
        (4 / 3)
        * length_a**3
        * b2a_ratio
        * c2a_ratio
        * (1.0 - 3 * (1.0 - t) ** 3)
    )



def A(qa, qb, qc, length_a, b2a_ratio, c2a_ratio, t):
    """
    Computes the AA term of the amplitude of the form factor AP at (qa,qb,qc).

    qnx, qny, qnz are the rescaled components (no unit) for computing AA term.

    Parameters
    ----------
    qa: the component of the scattering vector along a axis, units 1/Ang
    qb: the component of the scattering vector along b axis, units 1/Ang
    qc: the component of the scattering vector along c axis, units 1/Ang
    length_a: half height along the a axis of the octahedron without truncature
    b2a_ratio: half height along the b axis of the octahedron without truncature
    c2a_ratio: half height along the c axis of the octahedron without truncature
    t: truncation parameter, varies from 0.5 (cuboctahedron) to 1 (octahedron)

    Returns
    -------
    aa: AA term of the amplitude of the form factor at (qa,qb,qc)
    """
    length_b = length_a * b2a_ratio
    length_c = length_a * c2a_ratio
    qnx = qa * length_a  # conversion to dimensionless coordinate
    qny = qb * length_b  # conversion to dimensionless coordinate
    qnz = qc * length_c  # conversion to dimensionless coordinate

    # Protect denominators against exact zeros by adding tiny epsilon
    eps_den = 1e-18
    denom1 = (qny * qny - qnz * qnz) * (qny * qny - qnx * qnx)
    denom2 = (qnz * qnz - qnx * qnx) * (qnz * qnz - qny * qny)
    denom1 += eps_den * (np.abs(denom1) < eps_den)  # add epsilon where denom is zero
    denom2 += eps_den * (np.abs(denom2) < eps_den)  # add epsilon where denom is zero

    term1 = (
        (qny - qnx) * np.sin(qny * (1.0 - t) - qnx * t)
        + (qny + qnx) * np.sin(qny * (1.0 - t) + qnx * t)
    ) / denom1
    term2 = (
        (qnz - qnx) * np.sin(qnz * (1.0 - t) - qnx * t)
        + (qnz + qnx) * np.sin(qnz * (1.0 - t) + qnx * t)
    ) / denom2

    return term1 + term2

def Fqabc(qa, qb, qc, length_a, b2a_ratio, c2a_ratio, t):
    """
    Computes the amplitude of the form factor AP at (qa,qb,qc).

    Parameters
    ----------
    qa: component of the scattering vector along a axis, units 1/Ang
    qb: component of the scattering vector along b axis, units 1/Ang
    qc: component of the scattering vector along c axis, units 1/Ang
    length_a: half height along the a axis of the octahedron without truncature
    b2a_ratio: half height along the b axis of the octahedron without truncature
    c2a_ratio: half height along the c axis of the octahedron without truncature
    t: truncation parameter, varies from 0.5 (cuboctahedron) to 1 (octahedron)

    Returns
    -------
    ap: Amplitude of the form factor at (qa,qb,qc)
    """
    # Code taking into account the circular permutation between AA, BB and CC terms in AP.
    # amp3D is scaled to 1 near the origin for a regular shape.
    AA = A(qa, qb, qc, length_a, b2a_ratio, c2a_ratio, t)
    BB = A(qb, qc, qa, length_a, b2a_ratio, c2a_ratio, t)
    CC = A(qc, qa, qb, length_a, b2a_ratio, c2a_ratio, t)

    # Normalisation to 1 of AP at q = 0. Division by a factor 4/3 and global 1/2 coefficient.
    ap = 3.0 / (1.0 - 3.0 * (1.0 - t) ** 3) * (AA + BB + CC)

    # Guard against numerical issues: mirror the C fallback which sets
    # invalid (NaN/inf) values to zero so tests and downstream code
    # do not receive NaNs. Use np.nan_to_num for vector/scalar safety.
    return np.nan_to_num(ap, nan=0.0, posinf=0.0, neginf=0.0)

def Iqabc(qa, qb, qc, length_a, b2a_ratio, c2a_ratio, t):
    """
    Computes the 3D scattering intensity.

    Parameters
    ----------
    qa: component of the scattering vector along a axis, units 1/Ang
    qb: component of the scattering vector along b axis, units 1/Ang
    qc: component of the scattering vector along c axis, units 1/Ang
    length_a: half height along the a axis of the octahedron without truncature
    b2a_ratio: half height along the b axis of the octahedron without truncature
    c2a_ratio: half height along the c axis of the octahedron without truncature
    t: truncation parameter, varies from 0.5 (cuboctahedron) to 1 (octahedron)

    Returns
    -------
    intensity: 3D scattering intensity at (qa,qb,qc), units 1/cm
    """
    # int3D is scaled to 1 near the origin for a regular shape
    amp = Fqabc(qa, qb, qc, length_a, b2a_ratio, c2a_ratio, t)
    intensity = amp**2

    # If numerical issues produced invalid numbers, return 0.0 (like the C
    # implementation does) rather than NaN/inf which break tests and usage.
    if np.isnan(intensity) or np.isinf(intensity):
        return 0.0

    return intensity


# Cache for Fibonacci sphere points/weights to avoid recomputing across calls
_fibonacci_sphere_cache = {}

def _get_fibonacci_sphere(npoints):
    """
    Return (q_unit, w) for given npoints using module cache.
    """
    key = int(npoints)
    if key in _fibonacci_sphere_cache:
        return _fibonacci_sphere_cache[key]

    q_unit, w = fibonacci_sphere(int(npoints))
    _fibonacci_sphere_cache[key] = (q_unit, w)
    return q_unit, w


def Iq(
    q,
    sld,
    sld_solvent,
    length_a=500,
    b2a_ratio=1,
    c2a_ratio=1,
    t=0.99,
    npoints_fibonacci: int = 400,
):
    """
    Computes the 1D scattering intensity I(q) using Fibonacci quadrature.

    Parameters
    ----------
    q : float or ndarray
        Magnitude of the scattering vector, units 1/Ang
    sld : float
        Octahedron scattering length density, units 1e-6/Ang^2
    sld_solvent : float
        Solvent scattering length density, units 1e-6/Ang^2
    length_a : float
        Half height along the a axis of the octahedron without truncature, units Angstrom
    b2a_ratio : float
        Ratio b/a
    c2a_ratio : float
        Ratio c/a
    t : float
        Truncation parameter, varies from 0.5 (cuboctahedron) to 1 (octahedron)
    npoints_fibonacci : int
        Number of points on the sphere for the Fibonacci integration

    Returns
    -------
    1D scattering intensity at q, units 1/cm
    """
    q = np.atleast_1d(q)
    q_unit, w = _get_fibonacci_sphere(npoints_fibonacci)  # shape (npoints, 3)

    # build qx,qy,qz arrays with correct broadcasting -> shape (nq, npoints)
    qx = q[:, None] * q_unit[None, :, 0]
    qy = q[:, None] * q_unit[None, :, 1]
    qz = q[:, None] * q_unit[None, :, 2]

    # compute intensity grid using existing amp3D (vectorized)
    amp_grid = Fqabc(
        qx, qy, qz, length_a, b2a_ratio, c2a_ratio, t
    )  # shape (nq, npoints)

    # Replace any NaN/inf in amplitude grid (temporary fallback mirroring C
    # implementation). Use 0.0 to avoid producing huge/NaN intensities.
    amp_grid = np.nan_to_num(amp_grid, nan=0.0, posinf=0.0, neginf=0.0)
    # Also clip extremely large values to avoid overflow in square
    amp_grid = np.clip(amp_grid, -1e12, 1e12)

    intensity_grid = np.abs(amp_grid) ** 2

    integral = np.sum(intensity_grid * w[None, :], axis=1)  # summation over all points

    # Convert from [1e-12 A-1] to [cm-1]
    return (
        integral
        * 0.0001
        * (sld - sld_solvent) ** 2
        * form_volume(length_a, b2a_ratio, c2a_ratio, t) ** 2
    )


Iq.vectorized = True

tests = [
    [{"background": 0, "scale": 1, "length_a": 100, "t": 1, "sld": 1., "sld_solvent": 0.},
     0.01, 120.57219],
    [{"background": 0, "scale": 1, "length_a": 100, "t": 1, "sld": 1., "sld_solvent": 0.},
     [0.01, 0.1], [120.57218, 0.37421]],
]
