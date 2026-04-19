r"""
This model provides the form factor for truncated tetrahedrons.
Orientation averaging is done by using the Fibonacci quadrature.
This quadrature provides a quasi-uniform distribution of points on the unit sphere
using the golden ratio. The number of points to generate on the unit sphere
is set to 500 points, it usually provides a good balance between accuracy
and computational efficiency.

Definition
----------

To define the truncated tetrahedron form factor, the regular tetrahedron form factor has to be defined first.
The truncated tetrahedron form factor will then be obtained by subtracting the four smaller tetrahedrons at its vertices.
 So first, let's consider a regular tetrahedron. The size of the tetrahedron is described by its circumradius :math:`R`,
which is the radius of the circumscribed sphere. The relationship between the circumradius
and the edge length is also implemented. The edge length :math:`L`
and volume :math:`V_T` are given by:

.. math::

    L = \frac{4}{\sqrt{6}} \, R

.. math::

    V_{T}  = \frac{\sqrt{2}}{12} \, L^3

.. figure:: img/tetrahedron_regular.png

    Regular tetrahedron in its reference orientation.

The four vertices of the tetrahedron, with :math:`\mathbf{v}_0` at the origin,
are defined as:
.. math::

    \mathbf{v}_0 = (0,\ 0,\ 0)

.. math::

    \mathbf{v}_1 =  (\frac{L}{\sqrt{2}},\ \frac{L}{\sqrt{2}},\ 0\right)

.. math::

    \mathbf{v}_2 = (0,\ \frac{L}{\sqrt{2}},\frac{L}{\sqrt{2}})

.. math::

    \mathbf{v}_3 = (\frac{L}{\sqrt{2}},\ 0,\frac{L}{\sqrt{2}})

The form factor amplitude of a regular tetrahedron is derived using the
projection method described by Yang *et al.* [#Yang2023]_. The tetrahedron
is defined by four vertices
:math:`(\mathbf{v}_0, \mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3)`,
where :math:`\mathbf{v}_0` is placed at the origin and serves
as the projection origin. The three remaining vertices are specified such
that the normal vector of each face points away from the tetrahedron.
The form factor amplitude is defined as the Fourier transform of the
particle shape:

.. math::

    F(\mathbf{q}) = \int_V \exp(i\mathbf{q} \cdot \mathbf{r}) \, d\mathbf{r}

The transformation matrix :math:`\mathbf{T}` is built from the three edge
vectors originating at :math:`\mathbf{v}_0`, with each vector forming one
column:

.. math::

    \mathbf{T} = \bigl[\, \mathbf{v}_1 \mid \mathbf{v}_2 \mid \mathbf{v}_3 \,\bigr]

where :math:`(x_i, y_i, z_i)` are the coordinates of vertex :math:`\mathbf{v}_i`.
Introducing the shorthand notation :math:`Q_i = \mathbf{q} \cdot \mathbf{v}_i`, the
general tetrahedral form factor is given by:

.. math::

    F_\mathrm{T}(\mathbf{q}) = |\det(\mathbf{T})| \left\{
        \frac{i}{Q_1 (Q_1 - Q_2)(Q_1 - Q_3)} \exp(i Q_1)
      + \frac{i}{Q_2 (Q_2 - Q_1)(Q_2 - Q_3)} \exp(i Q_2)
      + \frac{i}{Q_3 (Q_3 - Q_2)(Q_3 - Q_1)} \exp(i Q_3)
      - \frac{i}{Q_1 Q_2 Q_3}
    \right\}


The truncated tetrahedron is a polyhedron obtained by truncating the vertices of a regular tetrahedron.

.. figure:: img/truncated_tetrahedron.png
.. figure:: img/tetra-octa-center-V1V2V3.png

    Tetrahedron with different truncatures.

The volume of the truncated tetrahedron is given by:

.. math::

    V_{T_{truncated}}  = V_{T}  (1 - 4t^3)

The parameter :math:`t` represents the truncation level. At the maximum value :math:`t=1/2`,
the corresponding shape is a regular octahedron. At the minimum value :math:`t=0`, 
the shape is a regular tetrahedron. At :math:`t=1/3`, the corresponding shape is the
Friauf polyhedron, in which all edges are equal.
The truncated tetrahedron form factor is also obtained by subtracting the four
smaller tetrahedrons at its vertices. It is expressed as:    

.. math::

    F_{T_{truncated}} (\vec{q},t, \vec{v_1},\vec{v_2},\vec{v_3}) 
    = F_{T} (\vec{q},\vec{v_1},\vec{v_2},\vec{v_3})-F_{T} (\vec{q},t\vec{v_1},t\vec{v_2},t\vec{v_3})(1+e^{i(1-t)\vec{q} \cdot \vec{v_1}}+e^{i(1-t)\vec{q} \cdot \vec{v_2}}+\ e^{i(1-t)\vec{q} \cdot \vec{v_3}})\

where the phase terms correspond to four translations inside the shape. 

Singularities: the expression of the regular tetreahedron form factor presents numerical singularities in four distinct cases, each handled analytically.

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

Orientation average: The 1D form factor corresponds to the orientation
average with all the possible orientations having the same probability.
Instead of rotating the shape through all the possible orientations,
it is equivalent to integrate the 3D scattering vector over a sphere
of radius q with the shape in its reference orientation.

The sphere is sampled using Fibonacci quadrature to provide a quasi-uniform
distribution of points on the unit sphere.
The repartition of the N points is computed using the golden ratio
(see fibonacci.py). Each point of the quadrature on the unit sphere
correspond to a vector :math:`\mathbf{u}_{j}`. In the sum,
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

The model has been tested against experimental data obtained on gold
tetrahedrons at synchrotron Soleil.
Moreover, comparisons with Debye formula calculations were made using
DebyeCalculator library (https://github.com/FrederikLizakJohansen/DebyeCalculator).
Good agreement was found at q < 1 1/Angstrom.

References
----------

1. Patterson, A. L. (1939)
   "The Diffraction of X-Rays by Small Crystalline Particles". 
   In: Phys. Rev. 56 (1939), pp. 972–977. 
   https://link.aps.org/doi/10.1103/PhysRev.56.972.

2. Wei-Ren Chen et al. "Scattering functions of Platonic solids".
   In: Journal of Applied Crystallography - J APPL CRYST 44 (June 2011).
   DOI: 10.1107/S0021889811011691

3. Croset, B. (2017) 
   Form factor of any polyhedron: a general compact formula and its singularities. 
   Journal of Applied Crystallography, 50(5), 1245-1255. 
   https://doi.org/10.1107/s1600576717010147.

4. Yang, T., Chen, X., Zhang, J., Ma, J., & Liu, S. (2023).
   Form factor of any polyhedron and its singularities derived from a projection method. 
   Journal of Applied Crystallography, 56(1), 167–177.
   DOI:  https://doi.org/10.1107/s160057672201130x


Authorship and Verification
----------------------------

* **Authors:** Chi-Huan Tung (tungc@ornl.gov)
             Marianne Imperor-Clerc (marianne.imperor@cnrs.fr)
             Sara Mokhtari (smokhtari@insa-toulouse.fr)

* **Last Modified by:** SM **Date:** 17 April 2026

* **Last Reviewed by:** MIC **Date:** 15 April 2026

"""
import numpy as np
from numpy import inf

from sasmodels.special.fibonacci import fibonacci_sphere

name = "tetrahedron_truncated"
title = "Tetrahedron truncated"
description = """
        Model for a truncated tetrahedron with orientation average using the Fibonacci quadrature"""
category = "shape:polyhedron"
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld", "1e-6/Ang^2", 126., [-inf, inf], "sld",
               "tetrahedron scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 9.4, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["circumradius", "Ang", 100, [0., inf], "volume",
               "Circumradius of the full tetrahedron"],
              ["t", "", 0, [0, 0.5], "volume",
               "truncation"],
               ]

### Functions for geometrical calculations:
# volume, edge length and vertices of the tetrahedron.
def form_volume(circumradius, t):
    """
    Computes the volume of the truncated tetrahedron from its circumradius.
    Parameters
    ----------
    circumradius : float
        Circumradius of the full tetrahedron.
    t : float
        Truncation parameter.
    Returns
    -------
    volume : float
        Volume of the truncated tetrahedron.
    """
    edge = edge_from_circumradius(circumradius)
    volume = edge**3 * np.sqrt(2) / 12 * (1 - 4*t**3)
    return volume

def edge_from_circumradius(circumradius):
    """
    Computes the edge length of the full tetrahedron from its circumradius.
    Parameters
    ----------
    circumradius : float
        Circumradius of the full tetrahedron.
    Returns
    -------
    edge : float
        Edge length of the full tetrahedron.
    """
    return 4 / np.sqrt(6) * circumradius

def circumradius_from_edge(edge):
    """
    Computes the circumradius of the full tetrahedron from its edge length.
    Parameters
    ----------
    edge : float
        Edge length of the full tetrahedron.
    Returns
    -------
    circumradius : float
        Circumradius of the full tetrahedron.
    """
    return np.sqrt(6) / 4 * edge

def vertices(circumradius):
    """
    Computes the vertices of the full tetrahedron from its circumrdius.
    Parameters
    ----------
    circumradius : float
        Circumradius of the full tetrahedron.
    Returns
    -------
    vertices : list
        List of the coordinates of the vertices of the full tetrahedron.
    """
    edge = edge_from_circumradius(circumradius)
    # The origin is located at V0
    return np.array([[0, 0, 0], [edge/np.sqrt(2), edge/np.sqrt(2), 0],
                     [ 0, edge/np.sqrt(2), edge/np.sqrt(2)],
                     [ edge/np.sqrt(2), 0,  edge/np.sqrt(2)]])

# Form factor and intensity calculations
# Better to define first an intermediate function F_ortho(qa,qb,qc,V1,V2,V3)
# for any tetrahedron shape
# then Fqabc(qa,qb,qc, circumradius) will call directly F_ortho(qa,qb,qc,V1,V2,V3)
def Fq_ortho(qa, qb, qc, V1, V2, V3):
    """
    General function that computes the form factor amplitude
    of the tetrahedra.
    Takes in account the singularities:
     - q perpendicular to a vertex,
     - q perpendicular to a an edge
     - q perpendicular to a to a face
     - q tending to 0.
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
    row9 = (Q1**2 + Q2**2 + Q3**2 < 1e-6)
    scattering_amplitude[row9] = (
        det_T / 6 +
        1j * det_T * (Q1[row9] + Q2[row9] + Q3[row9]))

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


def Fqabc(qa, qb, qc, circumradius, t):
    """
    Computes the form factor amplitude of the truncated tetrahedron
    from its circumradius and the three components of the
    scattering vector q.
    Parameters
    ----------
    circumradius : float
        Circumradius of the tetrahedron.
    qa, qb, qc : float
        Components of the scattering vector q
    t : float
        Truncation parameter.
    Returns
    -------
    Fqabc :
        Form factor amplitude of the truncated tetrahedron at the
        specific three dimensional q.
    """
    V0, V1, V2, V3 = vertices(circumradius)
    scattering_amplitude =  Fq_ortho(qa, qb, qc, V1, V2, V3) - Fq_ortho(qa, qb, qc, V1*t, V2*t, V3*t) * (
                   (1 + np.exp(1j * (1-t) * (qa * V1[0] + qb * V1[1] + qc * V1[2])))+
                   np.exp(1j * (1-t) * (qa * V2[0] + qb * V2[1] + qc * V2[2])) +
                   np.exp(1j * (1-t) * (qa * V3[0] + qb * V3[1] + qc * V3[2])))

    return scattering_amplitude


def Iqabc(qa, qb, qc, circumradius, t): # proportionnal to the volume**2
    """
    Computes the scattered intensity of the truncated tetrahedron at the specific
    three dimensional q components from the form factor amplitude Fqabc.
    Parameters
    ----------
    qa, qb, qc : float or array
        Components of the scattering vector q
    circumradius : float
        Circumradius of the full tetrahedron.
    t : float
        Truncation parameter.
    Returns
    -------
    Iqabc : float or array
        Scattered intensity of the tetrahedron at the specific three
        dimensional q components.
    """
    A = Fqabc(qa, qb, qc, circumradius, t)
    intensity = (np.abs(A))**2  # intensity is proportional to the volume
    return intensity

def Iq(q, sld, sld_solvent, circumradius:float, t:float, npoints_fibonacci:int= 500):
    """
    Computes the scattering intensity I(q) of the truncated tetrahedrons averaged over all
    orientations using the Fibonacci quadrature.
    The number of points on the sphere is set by npoints_fibonacci.
    Each point has an equal weight = 1/npoints_fibonacci.
    Parameters
    ----------
    q : float or array
        Norm of the scattering vector
    sld, sld_solvent :
        Contrast of scattering length density
    circumradius : float
        Circumradius of the full tetrahedron.
    t : float
        Truncation parameter.
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
    intensity = Iqabc(qa, qb, qc, circumradius, t)  # shape (nq, npoints)
    # Uniform average over the sphere
    integral = np.mean(intensity, axis=1)
    return (integral) * (sld - sld_solvent)**2 * 10**-4

Iq.vectorized = True

tests = [
    [{"background": 0, "scale": 1, "circumradius": 100, "t": 0, "sld": 1., "sld_solvent": 0.},
     0.01, 48.00520],
    [{"background": 0, "scale": 1, "circumradius": 100, "t": 0, "sld": 1., "sld_solvent": 0.},
     [0.01, 0.1], [48.00520, 0.57463]],
    [{"background": 0, "scale": 1, "circumradius": 100, "t": 0.5, "sld": 1., "sld_solvent": 0.},
     0.0001, 25.65146],
    [{"background": 0, "scale": 1, "circumradius": 100, "t": 0.3, "sld": 1., "sld_solvent": 0.},
     0.0001, 45.75264],

]

