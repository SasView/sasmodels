r"""
This model provides the form factor for a bifrustum whose cross-section is a regular polygon. 
Orientation averaging is done by using the Fibonacci quadrature.
This quadrature provides a quasi-uniform distribution of points on the unit sphere
using the golden ratio. The number of points to generate on the unit sphere is set to 500 points, it usually provides
a good balance between accuracy and computational efficiency.

Definition
----------

**Geometry.** The bifrustum is a bipyramid with a regular polygonal cross-section, truncated at both ends.
It is defined by the number of sides :math:`n` (n_sides), the circumradius :math:`R` of the central polygon
(radius_polygon), the half-height :math:`H` of the bipyramid before truncation (half_height) and the truncation
ratio :math:`t` (truncation) with :math:`t = 0` giving a sharp bipyramid and :math:`t = 1` a fully truncated bipyramid to 0 height. The bifrustum has :math:`3n` vertices.

.. figure:: img/bifrustum-doc.png

    Bifrustum geometry: pentagonal cross-section and truncation ratio of 0 and 0.25.


The scattering amplitude of a homogeneous particle is

.. math::

    F(\mathbf{q}) = \int_{V} e^{i\mathbf{q}\cdot\mathbf{r}}\, d^3\mathbf{r}

**Simplex method: decomposition into tetrahedra.** This integral is additive, so the amplitude of the bifrustum is the sum of the amplitudes of non-overlapping
tetrahedra filling its volume:

.. math::

    F_{\text{bifrustum}}(\mathbf{q}) = \sum_{\text{tetrahedra}} F_{\text{tetra}}(\mathbf{q})

The tetrahedra are obtained by a Delaunay triangulation of the :math:`3n` vertices (scipy.spatial.Delaunay).
Delaunay only connects the existing vertices: since the bifrustum is convex, every tetrahedron has its four
corners among the vertices of the bifrustum, even when it crosses the inside of the particle.
Degenerate tetrahedra (volume below :math:`10^{-6}` |Ang^3|) are discarded.
This approach applies to any convex polyhedron defined by its vertices.


.. figure:: img/simplex_method.png

    Bifrustum divided into tetrahedra.

**Form factor of a tetrahedron.** Let :math:`\mathbf{r}_0, \mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_3` be the
vertices of a tetrahedron and :math:`\mathbf{V}_k = \mathbf{r}_k - \mathbf{r}_0` (:math:`k = 1, 2, 3`) its edges.
Any point inside is written :math:`\mathbf{r} = \mathbf{r}_0 + u\mathbf{V}_1 + v\mathbf{V}_2 + w\mathbf{V}_3` with
:math:`u, v, w \geq 0` and :math:`u + v + w \leq 1`. This change of variables maps the integral onto a reference
tetrahedron, with Jacobian :math:`|\det T|`, where :math:`T = (\mathbf{V}_1\ \mathbf{V}_2\ \mathbf{V}_3)`
(:math:`|\det T|` is six times the volume of the tetrahedron). With the projections of :math:`\mathbf{q}` on the edges

.. math::

    Q_k = \mathbf{q}\cdot\mathbf{V}_k

The form factor amplitude of a regular tetrahedron is derived using the projection method described by Yang and al. [1].

.. math::

    F_{\text{tetra}}(\mathbf{q}) = i\,|\det T|\, e^{i\mathbf{q}\cdot\mathbf{r}_0} \left[
        \frac{e^{iQ_1}}{Q_1(Q_1-Q_2)(Q_1-Q_3)}
      + \frac{e^{iQ_2}}{Q_2(Q_2-Q_1)(Q_2-Q_3)}
      + \frac{e^{iQ_3}}{Q_3(Q_3-Q_1)(Q_3-Q_2)}
      - \frac{1}{Q_1 Q_2 Q_3} \right]

where :math:`e^{i\mathbf{q}\cdot\mathbf{r}_0}` is the translation phase of the tetrahedron.
For :math:`\mathbf{q} \to 0`, :math:`F_{\text{tetra}} \to |\det T|/6`, the volume of the tetrahedron.

**Expression used in the code.** Since
:math:`e^{i\mathbf{q}\cdot\mathbf{r}_0} e^{iQ_k} = e^{i\mathbf{q}\cdot\mathbf{r}_k}`, the translation phase can be
moved inside the brackets:

.. math::

    F_{\text{tetra}}(\mathbf{q}) = i\,|\det T| \left[
        \frac{e^{i\mathbf{q}\cdot\mathbf{r}_1}}{Q_1(Q_1-Q_2)(Q_1-Q_3)}
      + \frac{e^{i\mathbf{q}\cdot\mathbf{r}_2}}{Q_2(Q_2-Q_1)(Q_2-Q_3)}
      + \frac{e^{i\mathbf{q}\cdot\mathbf{r}_3}}{Q_3(Q_3-Q_1)(Q_3-Q_2)}
      - \frac{e^{i\mathbf{q}\cdot\mathbf{r}_0}}{Q_1 Q_2 Q_3} \right]

The numerators only depend on the positions of the vertices, and the denominators are obtained from the same
projections since :math:`Q_k = \mathbf{q}\cdot\mathbf{r}_k - \mathbf{q}\cdot\mathbf{r}_0`.
Each vertex of the bifrustum belongs to several tetrahedra. :math:`\mathbf{q}\cdot\mathbf{r}_k`
are therefore computed once per vertex here, and reused by every tetrahedron that contains it
(instead of being recomputed for each tetrahedron).

**Singularities**: the expression of the regular tetreahedron form factor presents numerical singularities in four distinct cases, each handled analytically.

**Case 1** — :math:`\mathbf{q}` perpendicular to vertex :math:`\mathbf{v}_i`
(:math:`Q_i = 0`, with :math:`Q_j, Q_k \neq 0` and :math:`Q_j \neq Q_k`):

.. math::

    F_P(\mathbf{q}) = |\det(\mathbf{T})| \left\{
        \frac{i}{Q_j^2 (Q_j - Q_k)} \exp(i Q_j)
      + \frac{i}{Q_k^2 (Q_k - Q_j)} \exp(i Q_k)
      + \frac{i(Q_j + Q_k + i Q_j Q_k)}{Q_j^2 \, Q_k^2}
    \right\}

**Case 2** — :math:`\mathbf{q}` perpendicular to edge :math:`\mathbf{v}_i \mathbf{v}_j`
(:math:`Q_i = Q_j`, with :math:`Q_k \neq Q_i`):

.. math::

    F_E(\mathbf{q}) = |\det(\mathbf{T})| \left\{
        \frac{i}{Q_k (Q_k - Q_i)^2} \exp(i Q_k)
      + \frac{i Q_k - 2i Q_i - Q_i^2 + Q_i Q_k}{Q_i^2 (Q_k - Q_i)^2} \exp(i Q_i)
      - \frac{i}{Q_i^2 \, Q_k}
    \right\}

**Case 3** — :math:`\mathbf{q}` perpendicular to a face
(:math:`Q_i = Q_j = Q_k  \neq 0`):

.. math::

    F_F(\mathbf{q}) = |\det(\mathbf{T})| \left\{
        \frac{-i}{Q_i^3}
      + \frac{2i + 2Q_i - i Q_i^2}{2 Q_i^3} \exp(i Q_i)
    \right\}

The formula can be expressed by three alternatives
(i.e. :math:`Q_i = q · v_1` or :math:`Q_i = q · v_2`  or :math:`Q_i = q · v_3` ).

**Case 4** — :math:`\mathbf{q} \to 0` (first-order Taylor expansion):

.. math::

    F_V(\mathbf{q}) = V_\mathrm{T}
      + |\det(\mathbf{T})| \, \frac{i(Q_1 + Q_2 + Q_3)}{24}

where the tetrahedral volume is :math:`V_\mathrm{T} = |\det(\mathbf{T})| / 6`.

**Orientation average**: The 1D form factor corresponds to the orientation
average with all the possible orientations having the same probability.
Instead of rotating the shape through all the possible orientations,
it is equivalent to integrate the 3D scattering vector over a sphere
of radius q with the shape in its reference orientation.

The sphere is sampled using Fibonacci quadrature to provide a quasi-uniform
distribution of points on the unit sphere.
The distribution of the N points is computed using the golden ratio
(see fibonacci.py). Each point of the quadrature on the unit sphere
corresponds to a vector :math:`\mathbf{u}_{j}`. In the sum,
all weights :math:`w_j` are taken identical and equal to :math:`\frac{1}{N}`.

.. math::

    P(q) =  \sum_{j=1}^{N} w_j I(q\mathbf{u}_{j}, R, t)

.. figure:: img/fibonacci_sphere.png

    Fibonacci sphere using N=5810 points.

And the 1D scattering intensity is calculated as:

.. math::

    I(q) = \text{scale} \times V \times (\rho_\text{p} -
    \rho_\text{solvent})^2 \times P(q)

where V is the volume of the tetrahedron, ρ
is the scattering length inside the volume, ρ *solvent*
is the scattering length of the solvent, and (if the data are in absolute
units) *scale* represents the volume fraction (which is unitless).

Validation
----------

The model has been tested against experimental data obtained on truncated gold bipyramids with pentagonal cross-section synthetized by Claire Goldmann at LPS, Orsay.
Moreover, comparisons with Debye formula calculations were made using DebyeCalculator library (https://github.com/FrederikLizakJohansen/DebyeCalculator).
Good agreement was found at q < 0.1 1/Angstrom.

References
----------

1. Yang, T., Chen, X., Zhang, J., Ma, J., & Liu, S. (2023).
   Form factor of any polyhedron and its singularities derived from a projection method.
   *Journal of Applied Crystallography*, 56(1), 167–177. https://doi.org/10.1107/s160057672201130x


2. Li, X., Shew, C., He, L., Meilleur, F., Myles, D. A. A., Liu, E., Zhang, Y., Smith, G. S.,
   Herwig, K. W., Pynn, R., & Chen, W. (2011). Scattering functions of Platonic solids.
   *Journal Of Applied Crystallography*, 44(3), 545‑557. https://doi.org/10.1107/s0021889811011691

3. Croset, B. (2017). Form factor of any polyhedron : a general compact formula and its singularities.
   *Journal Of Applied Crystallography*, 50(5), 1245‑1255. https://doi.org/10.1107/s1600576717010147



Authorship and Verification
----------------------------

* **Authors:** Chi-Huan Tung (tungc@ornl.gov)
             Marianne Imperor-Clerc (marianne.imperor@cnrs.fr)
             Sara Mokhtari (smokhtari@insa-toulouse.fr)

* **Last Modified by:** SM **Date:** 04 November 2026

* **Last Reviewed by:** MIC **Date:** 04 November 2026

"""
import math

import numpy as np
from numpy import inf
from scipy.spatial import Delaunay

from sasmodels.special.fibonacci import fibonacci_sphere

name = "bifrustum"
title = "Bifrustum"
description = """
        Model for a bifrustum with orientation average using the Fibonacci quadrature"""
category = "shape:polyhedron"
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld", "1e-6/Ang^2", 126., [-inf, inf], "sld",
               "Prism scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 9.4, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["n_sides", "", 5, [3, 50], "volume",
               "Number of sides"],
              ["radius_polygon", "Ang", 50, [0., inf], "volume",
               "Circumradius of the polygonal cross-section"],
              ["half_height", "Ang", 250, [0., inf], "volume",
               "Half-height of the bifrustum (before truncation)"],
              ["truncation", "", 0.5, [0., 1.], "volume",
               "Truncation ratio, t = 0 → sharp bipyramid and t = 1 → fully truncated to zero height"],
              ["theta", "degrees", 0, [-180, 180], "orientation", "c axis to beam angle"],
              ["phi", "degrees", 0, [-360, 360], "orientation", "rotation about beam"],
              ["psi", "degrees", 0, [-360, 360], "orientation", "rotation about c axis"],
               ]

### Functions for geometrical calculations: volume and vertices of the bifrustum
def form_volume(n_sides, radius_polygon, half_height, truncation):
    """
    Computes the volume of a nanoprism given its number of sides, average radius and length.
    Parameters
    ----------
    n_sides : int
        Number of sides of the regular polygon cross-section.
    radius_polygon : float
        Circumradius of the regular polygon cross-section.
    half_height : float
        Half-height of the bipyramid (before truncation).
    truncation : float
        Truncation ratio of the bipyramid.
    Returns
    -------
    volume : float
        Volume of the nanoprism.
    """
    angle    = math.pi / n_sides
    r_top    = radius_polygon * truncation  # circumradius of the top polygon

    a_base   = 2 * radius_polygon * math.sin(angle)   # edge length of the base polygon
    a_top    = 2 * r_top * math.sin(angle)   # edge length of the top polygon
    apo_base = radius_polygon * math.cos(angle)       # apothem base
    apo_top  = r_top * math.cos(angle)       # apothem top

    A_base = 0.5 * n_sides * a_base * apo_base     # base polygon area
    A_top  = 0.5 * n_sides * a_top  * apo_top      # area of the top polygon

    V_bipyramid = (2/3) * A_base * half_height                  # complete bipyramid
    V_cap       = (1/3) * A_top  * (half_height * truncation)            # each cut pyramid
    volume      = V_bipyramid - 2 * V_cap

    return volume

def bifrustum_verts(n_sides, radius_polygon, half_height, truncation):
    """
    Generate the 3D coordinates of a bifrustum shape (truncated bipyramid).

    Parameters:
    - n_sides: int
        Number of sides in the base polygon (e.g., 5 for a pentagon).
    - radius_polygon: float
        Circumradius of the base polygon.
    - half_height: float
        Half-height from center to apex (before truncation).
    - truncation: float
        Truncation ratio:
          t = 0 → sharp bipyramid,
          t = 1 → fully truncated to zero height.

    Returns:
    - verts: (3n, 3) np.ndarray
        Coordinates of vertices: base (z=0), upper (z=+), and lower (z=-).
    """
    n_sides = int(n_sides)
    angles = np.linspace(0, 2 * np.pi, n_sides, endpoint=False)

    base = np.stack([radius_polygon * np.cos(angles), radius_polygon * np.sin(angles), np.zeros(n_sides)], axis=1)
    radius_truncated = radius_polygon * truncation
    z_up = half_height * (1 - truncation)
    up = np.stack([radius_truncated * np.cos(angles), radius_truncated * np.sin(angles), np.full(n_sides, z_up)], axis=1)
    down = np.stack([radius_truncated * np.cos(angles), radius_truncated * np.sin(angles), np.full(n_sides, -z_up)], axis=1)

    verts = np.vstack([base, up, down])
    return verts

# Functions for the form factor calculation of the bifrustum using the simplex method (Delaunay triangulation)
def delaunay(verts):
    """
    Compute the Delaunay triangulation of a set of vertices.

    Parameters:
    - verts: (3*n_sides, 3) np.ndarray of vertex coordinates.

    Returns:
    - tetrahedra: (M, 4) np.ndarray of indices of vertices forming each tetrahedron.
    """

    delaunay_result = Delaunay(verts)
    return delaunay_result.simplices

def compute_tetrahedra_volumes(verts):
    """
    Compute volume of each tetrahedron.

    Parameters:
    - verts: (3*n_sides, 3) np.ndarray of vertex coordinates.

    Returns:
    - volumes: (M,) np.ndarray of positive volumes.
    """
    tetrahedra = delaunay(verts)
    v0 = verts[tetrahedra[:, 0]]
    v1 = verts[tetrahedra[:, 1]]
    v2 = verts[tetrahedra[:, 2]]
    v3 = verts[tetrahedra[:, 3]]

    # Compute the signed volume of each tetrahedron using the scalar triple product
    cross = np.cross(v1 - v0, v2 - v0)
    dot = np.einsum('ij,ij->i', cross, v3 - v0)
    volumes = np.abs(dot) / 6.0

    volume_threshold = 1e-6

    # Filter out tetrahedra with volumes below the threshold
    valid_indices = np.where(volumes > volume_threshold)[0]

    # Filter tetrahedra and volumes
    filtered_tetrahedra = tetrahedra[valid_indices]
    filtered_volumes = volumes[valid_indices]

    total_volume = np.sum(filtered_volumes)

    return total_volume, filtered_tetrahedra, filtered_volumes


def fix_tetrahedron_orientations(verts, tetrahedra):
    """
    Ensure all tetrahedra are right-handed (positive volume).

    Parameters:
    - verts: (N, 3)
    - tetrahedra: (M, 4)

    Returns:
    - corrected_tetrahedra: (M, 4)
    """
    corrected = tetrahedra.copy()
    # Check the orientation of each tetrahedron
    for i, tet in enumerate(tetrahedra):
        v0, v1, v2, v3 = verts[tet]
        signed_vol = np.dot(np.cross(v1 - v0, v2 - v0), v3 - v0)
        if signed_vol < 0:
            corrected[i] = [tet[0], tet[2], tet[1], tet[3]]
    return corrected


def Fq_ortho(qa, qb, qc, V1, V2, V3):
    """
    General function that computes the form factor amplitude
    of the tetrahedra.
    Takes in account the singularities:
     - 1.1, 1.2, 1.3 q perpendicular to a vertex,
     - 2.1, 2.2, 2.3 q perpendicular to a an edge
     - 3. q perpendicular to a to a face (face is opposite to V0)
     - 4. q tending to 0
     Some cases are not covered in the original paper and are added here:
     - 5. q perpendicular to a face containing V0 (two Q_i = 0)
     - 6. q perpendicular to two opposite edges (Q_i = 0 and Q_j = Q_k).
    Code written by Tianjuan Yang, adapted by Sara Mokhtari.
    https://doi.org/10.1107/S160057672201130X
    Parameters
    ----------
    qa, qb, qc : float
        Components of the scattering vector q
    V1, V2, V3 : array-like
        Coordinates of the three vertices of the tetrahedron
        (V0 is at the origin and not used).
    Returns
    -------
    scattering_amplitude :
        Form factor amplitude of the tetrahedron at the
        specific three dimensional q.
    """

    T = np.array([V1, V2, V3]).T
    Q1 = qa * V1[0] + qb * V1[1] + qc * V1[2]
    Q2 = qa* V2[0] + qb * V2[1] + qc * V2[2]
    Q3 = qa * V3[0] + qb * V3[1] + qc * V3[2]
    det_T = np.linalg.det(T)
    scattering_amplitude = np.zeros_like(Q1, dtype=complex)

    # Calculation at singularities
    # 1.1 q perpendicular to vertex V1 (Q_i = 0, here Q_i is Q_1)
    row2 = (
        (np.abs(Q1) <= 1e-9) & (np.abs(Q2) >= 1e-9) &
        (np.abs(Q3) >= 1e-9) & (np.abs(Q2 - Q3) >= 1e-9))
    scattering_amplitude[row2] = det_T * (
        1j * np.exp(1j * Q2[row2]) /
        (Q2[row2]**2 * (Q2[row2] - Q3[row2])) +
        1j * np.exp(1j * Q3[row2]) /
        (Q3[row2]**2 * (Q3[row2] - Q2[row2])) +
        1j * (Q2[row2] + Q3[row2] +
              1j * Q2[row2] * Q3[row2]) /
        (Q2[row2]**2 * Q3[row2]**2))

    # 1.2 q perpendicular to vertex V2 (Q_i = 0, here Q_i is Q_2)
    row3 = (
        (np.abs(Q2) <= 1e-9) & (np.abs(Q1) >= 1e-9) &
        (np.abs(Q3) >= 1e-9) & (np.abs(Q1 - Q3) >= 1e-9))
    scattering_amplitude[row3] = det_T * (
        1j * np.exp(1j * Q1[row3]) /
        (Q1[row3]**2 * (Q1[row3] - Q3[row3])) +
        1j * np.exp(1j * Q3[row3]) /
        (Q3[row3]**2 * (Q3[row3] - Q1[row3])) +
        1j * (Q1[row3] + Q3[row3] +
              1j * Q1[row3] * Q3[row3]) /
        (Q1[row3]**2 * Q3[row3]**2))

    # 1.3 q perpendicular to vertex V3 (Q_i = 0, here Q_i is Q_3)
    row4 = (
        (np.abs(Q3) <= 1e-9) & (np.abs(Q1) >= 1e-9) &
        (np.abs(Q2) >= 1e-9) & (np.abs(Q1 - Q2) >= 1e-9))
    scattering_amplitude[row4] = det_T * (
        1j * np.exp(1j * Q1[row4]) /
        (Q1[row4]**2 * (Q1[row4] - Q2[row4])) +
        1j * np.exp(1j * Q2[row4]) /
        (Q2[row4]**2 * (Q2[row4] - Q1[row4])) +
        1j * (Q1[row4] + Q2[row4] +
              1j * Q1[row4] * Q2[row4]) /
        (Q1[row4]**2 * Q2[row4]**2))

    # 2.1 q perpendicular to an edge (Q_i = Q_j, here Q_i is Q_1 and Q_j is Q_3)
    row5 = (
        (np.abs(Q1 - Q3) < 1e-9) & (np.abs(Q1) >= 1e-9) &
        (np.abs(Q1 - Q2) >= 1e-9) & (np.abs(Q3) >= 1e-9))
    scattering_amplitude[row5] = det_T * (
        1j * np.exp(1j * Q2[row5]) /
        (Q2[row5] * (Q2[row5] - Q1[row5])**2) +
        np.exp(1j * Q1[row5]) *
        (1j * Q2[row5] - 2j * Q1[row5] -
         Q1[row5]**2 + Q1[row5] * Q2[row5]) /
        ((Q2[row5] - Q1[row5])**2 * Q1[row5]**2) -
        1j / (Q1[row5]**2 * Q2[row5]))

    # 2.2 q perpendicular to an edge (Q_i = Q_j, here Q_i is Q_1 and Q_j is Q_2)
    row6 = (
        (np.abs(Q1 - Q2) < 1e-9) & (np.abs(Q1) >= 1e-9) &
        (np.abs(Q2 - Q3) >= 1e-9) & (np.abs(Q2) >= 1e-9))
    scattering_amplitude[row6] = det_T * (
        1j * np.exp(1j * Q3[row6]) /
        (Q3[row6] * (Q3[row6] - Q1[row6])**2) +
        np.exp(1j * Q1[row6]) *
        (1j * Q3[row6] - 2j * Q1[row6] -
         Q1[row6]**2 + Q1[row6] * Q3[row6]) /
        ((Q3[row6] - Q1[row6])**2 * Q1[row6]**2) -
        1j / (Q1[row6]**2 * Q3[row6]))

    # 2.3 q perpendicular to an edge (Q_i = Q_j, here Q_i is Q_2 and Q_j is Q_3)
    row7 = (
        (np.abs(Q2 - Q3) < 1e-9) & (np.abs(Q2) >= 1e-9) &
        (np.abs(Q1 - Q3) >= 1e-9) & (np.abs(Q3) >= 1e-9))
    scattering_amplitude[row7] = det_T * (
        1j * np.exp(1j * Q1[row7]) /
        (Q1[row7] * (Q1[row7] - Q2[row7])**2) +
        np.exp(1j * Q2[row7]) *
        (1j * Q1[row7] - 2j * Q2[row7] -
         Q2[row7]**2 + Q2[row7] * Q1[row7]) /
        ((Q1[row7] - Q2[row7])**2 * Q2[row7]**2) -
        1j / (Q2[row7]**2 * Q1[row7]))

    # 3. q perpendicular to a face (Q_i = Q_j = Q_k, here Q_i is Q_1, Q_j is Q_2 and Q_k is Q_3)
    # Here the face is the opposite to the vertex V0, the one at the origin
    row8 = (
        (np.abs(Q1 - Q2) <= 1e-5) &
        (np.abs(Q1 - Q3) <= 1e-5) &
        (np.abs(Q2 - Q3) <= 1e-5) &
        (np.abs(Q1) >= 1e-9) & (np.abs(Q2) >= 1e-9) &
        (np.abs(Q3) >= 1e-9))
    scattering_amplitude[row8] = det_T * (
        np.exp(1j * Q1[row8]) *
        (2j + 2 * Q1[row8] - 1j * Q1[row8]**2) /
        (2 * Q1[row8]**3) - 1j / Q1[row8]**3)

    # 4. q tending to 0 (first-order Taylor expansion)
    # The integral of u (or v, w) over the reference tetrahedron is 1/24
    row9 = (Q1**2 + Q2**2 + Q3**2 < 1e-6)
    scattering_amplitude[row9] = (
        det_T / 6 +
        1j * det_T * (Q1[row9] + Q2[row9] + Q3[row9]) / 24)

    # Regular case of Fq_ortho
    regular = ~(row2 | row3 | row4 | row5 | row6 | row7 | row8 | row9)
    scattering_amplitude[regular] = det_T * (
        1j * np.exp(1j * Q3[regular]) /
        (Q3[regular] * (Q3[regular] - Q2[regular]) * (Q3[regular] - Q1[regular])) +
        1j * np.exp(1j * Q2[regular]) /
        (Q2[regular] * (Q2[regular] - Q1[regular]) * (Q2[regular] - Q3[regular])) +
        1j * np.exp(1j * Q1[regular]) /
        (Q1[regular] * (Q1[regular] - Q2[regular]) * (Q1[regular] - Q3[regular])) -
        1j / (Q1[regular] * Q2[regular] * Q3[regular]))

    return scattering_amplitude

# To use Fq_ortho only for points near a singular case,
# we define a margin that is wider than all the thresholds of Fq_ortho.
NEAR_SINGULAR = 1e-3


def Fqabc(qa, qb, qc, n_sides, radius_polygon, half_height, truncation):
    """
    Computes the form factor amplitude of the sum of all the
    tetrahedrons that compose the prism from the three components of the
    scattering vector q, the number of sides of the regular polygon cross-section
    and the circumradius of the regular polygon cross-section.
    Parameters
    ----------
    n_sides : int
        Number of sides of the regular polygon cross-section.
    radius_polygon : float
        Circumradius of the regular polygon cross-section.
    half_height : float
        Half the height of the bifrustum.
    truncation : float
        The truncation of the bifrustum.
    qa, qb, qc : float
        Components of the scattering vector q
    Returns
    -------
    Fqabc :
        Form factor amplitude of the tetrahedron at the
        specific three dimensional q.
    """

    # Generate the vertices of the bifrustum and compute the tetrahedra volumes
    verts = bifrustum_verts(n_sides, radius_polygon, half_height, truncation)
    _ , filtered_tetrahedra, _ = compute_tetrahedra_volumes(verts)

    # Ensure that all tetrahedra are right-handed (positive volume)
    tetrahedra_fixed = fix_tetrahedron_orientations(verts, filtered_tetrahedra)

    # New optimization: shared vertices.
    # Each vertex of the bifrustum belongs to several tetrahedra. q·r_k and e^{i q·r_k} are therefore
    # computed once per vertex here, and reused by every tetrahedron that contains it
    # (instead of being recomputed for each tetrahedron).
    P = [qa * v[0] + qb * v[1] + qc * v[2] for v in verts]  # P[k] = q·r_k
    E = [np.exp(1j * p) for p in P]                         # E[k] = e^{i q·r_k}
    F = np.zeros_like(qa, dtype=complex)

    for tet in tetrahedra_fixed:
        i0, i1, i2, i3 = tet
        Q1, Q2, Q3 = P[i1] - P[i0], P[i2] - P[i0], P[i3] - P[i0]  # Q_k = q·(r_k - r_0)
        Q12, Q13, Q23 = Q1 - Q2, Q1 - Q3, Q2 - Q3

        # New optimization: only the points near a singular case are sent to Fq_ortho.

        # First, the regular case of Fq_ortho, written with the absolute phases of the vertices so that E is reused:
        # F = i·det_T·[ e^{iP1}/(Q1·Q12·Q13) - e^{iP2}/(Q2·Q12·Q23) + e^{iP3}/(Q3·Q13·Q23) - e^{iP0}/(Q1·Q2·Q3) ]
        # It is computed at every point; the points near a singular case will be replaced after.
        det_T = np.linalg.det(verts[[i1, i2, i3]] - verts[i0])  # positive after fix_tetrahedron_orientations
        with np.errstate(divide='ignore', invalid='ignore'): # if any of the denominators is zero, the result will be NaN or inf, which will be replaced by Fq_ortho
            f = 1j * det_T * (E[i1] / (Q1 * Q12 * Q13)
                              - E[i2] / (Q2 * Q12 * Q23)
                              + E[i3] / (Q3 * Q13 * Q23)
                              - E[i0] / (Q1 * Q2 * Q3))

        # Points near a singular case: one of the denominators Q1, Q2, Q3, Q12, Q13, Q23 is small.
        # This is deliberately a single, wide margin (NEAR_SINGULAR), not the exact conditions of Fq_ortho:
        # these points are sent to Fq_ortho, which handles all cases (regular included) and is the only
        # place where the singular cases and their thresholds are decided. The margin only has to cover
        # them (Fq_ortho uses 1e-9 for vertices/edges, 1e-5 for faces and |Q|² < 1e-6 for q -> 0).
        smallest = np.minimum.reduce([np.abs(Q1), np.abs(Q2), np.abs(Q3), # take the minimum of all the denominators
                                      np.abs(Q12), np.abs(Q13), np.abs(Q23)])
        near = smallest < NEAR_SINGULAR # larger than the thresholds of Fq_ortho
        if near.any():
            # Fq_ortho places r_0 at the origin, hence the translation phase e^{i q·r_0} = E[i0]
            r0 = verts[i0]
            V1, V2, V3 = verts[i1] - r0, verts[i2] - r0, verts[i3] - r0
            # (some cases of Fq_ortho overlap and are overwritten by a later case: ignore their warnings)
            with np.errstate(divide='ignore', invalid='ignore'):
                f[near] = Fq_ortho(qa[near], qb[near], qc[near], V1, V2, V3) * E[i0][near]

        F += f

    return F


def _Iqabc(qa, qb, qc, n_sides, radius_polygon, half_height, truncation): # proportional to the volume**2
    """
    Computes the scattered intensity of the tetrahedron at the specific
    three dimensional q components from the form factor amplitude Fqabc.
    Parameters
    ----------
    qa, qb, qc : float or array
        Components of the scattering vector q
    n_sides : int
        Number of sides of the regular polygon cross-section.
    radius : float
        Circumradius of the regular polygon cross-section.
    Returns
    -------
    Iqabc : float or array
        Scattered intensity of the tetrahedron at the specific three
        dimensional q components.
    """
    A = Fqabc(qa, qb, qc, n_sides, radius_polygon, half_height, truncation)
    intensity = (np.abs(A))**2  # intensity is proportional to the volume
    return intensity

def Iqabc(qa, qb, qc, sld, sld_solvent, n_sides, radius_polygon, half_height, truncation): # proportionnal to the volume**2
    """
    Computes the scattered intensity of the truncated tetrahedron at the specific
    three dimensional q components from the form factor amplitude Fqabc.
    Parameters
    ----------
    qa, qb, qc : float or array
        Components of the scattering vector q
    sld, sld_solvent :
        Contrast of scattering length density
    radius_polygon : float
        Circumradius of the full tetrahedron.
    half_height : float
        Half-height of the bipyramid.
    truncation : float
        Truncation parameter.
    Returns
    -------
    Iqabc : float or array
        Scattered intensity of the tetrahedron at the specific three
        dimensional q components.
    """
    intensity = _Iqabc(qa, qb, qc, n_sides, radius_polygon, half_height, truncation)  # shape (nq, npoints)
    return intensity * ((sld - sld_solvent)**2 * 10**-4)

def Iq(q, sld, sld_solvent, n_sides: int, radius_polygon: float, half_height: float, truncation: float, npoints_fibonacci: int = 500):
    """
    Computes the scattering intensity I(q) of tetrahedrons averaged over all
    orientations using the Fibonacci quadrature.
    The number of points on the sphere is set by npoints_fibonacci.
    Each point has an equal weight = 1/npoints_fibonacci.
    Parameters
    ----------
    q : float or array
        Norm of the scattering vector
    sld, sld_solvent :
        Contrast of scattering length density
    n_sides : int
        Number of sides of the regular polygon cross-section.
    radius_polygon : float
        Circumradius of the regular polygon cross-section.
    half_height : float
        Half-height of the bipyramid.
    truncation : float
        Truncation ratio of the bipyramid.
    npoints_fibonacci : int
        Number of Fibonacci points on the sphere, set to 500 by default
        (higher number increases accuracy but also computation time,
          500 is usually a good compromise)
    Returns
    -------
    Iq : ndarray
        Scattering intensity averaged over all orientations
    """

    q = np.atleast_1d(q)
    q_unit, w = fibonacci_sphere(npoints_fibonacci)   # shape (npoints, 3)

    # Projections
    qa = q[:, np.newaxis] * q_unit[:, 0][np.newaxis, :]
    qb = q[:, np.newaxis] * q_unit[:, 1][np.newaxis, :]
    qc = q[:, np.newaxis] * q_unit[:, 2][np.newaxis, :]
    # Compute intensity
    intensity = _Iqabc(qa, qb, qc, n_sides, radius_polygon, half_height, truncation)
    # Uniform average over the sphere
    integral = np.mean(intensity, axis=1)

    return (integral) * (sld - sld_solvent)**2 * 10**-4

Iq.vectorized = True

tests = [
    [{"background": 0, "scale": 1, "n_sides": 4, "radius_polygon": 50, "half_height": 500, "truncation": 0.5, "sld": 1., "sld_solvent": 0.},
     0.01, 145.11465],
    [{"background": 0, "scale": 1, "n_sides": 4, "radius_polygon": 50, "half_height": 500, "truncation": 0.5, "sld": 1., "sld_solvent": 0.},
     [0.01, 0.1], [145.11465, 0.37072]],
]

