r"""
This model provides the form factor for a right prism whose cross-section is a regular polygon. 
Orientation averaging is done by using the Fibonacci quadrature.
This quadrature provides a quasi-uniform distribution of points on the unit sphere
using the golden ratio. The number of points to generate on the unit sphere is set to 500 points, it usually provides
a good balance between accuracy and computational efficiency.

Definition
----------

We consider particles having the shape of a right prism of length *L* (parameter called *length* in the model) and a cross section made of 
a regular polygon with a certain number of sides *n* (parameter called *n_sides* in the model) as illustrated in the figure below.

.. figure:: img/nanoprisms_geometry.jpg

The size of a regular polygon can be characterized by its edge length E or by R, the circumradius of the polygon shown in black in the figure, with :

.. math::

    R = \frac{E}{2\sin(\pi/n)}


The apothem is the radius of the dashed circle, it is defined as:

.. math::

    R\cos(\pi/n)

For comparison purposes, it is convenient to introduce an average radius :math:`R_ave` (parameter called *radius_average* in the model). It is shown in red in the figure.
The area of the n-sided regular polygon is given by :

.. math::

   A = \pi R^{2} \, \mathrm{sinc}\left(\frac{2\pi}{n}\right)
     = \pi R_{\mathrm{ave}}^{2}

and the volume of the nanoprism is therefore given by:

.. math::

    V = L \times {A}

where :math:`R_ave` is the radius of the equivalent disc having the same area as the n-sided polygon.
It is also the squared average of the distance from the center of the polygon to any point of its perimeter. It is related to the circumradius R by :

.. math::

    R_{\mathrm{ave}}^{2} = R^{2} \, \mathrm{sinc}\left(\frac{2\pi}{n}\right)

Form factor for a prism: Following Wuttke's expression, the form factor :math:`F(\mathbf{q})` for any right prism can be decomposed into the product of two factors.
One factor corresponds to the component :math:`\mathbf{q}_{\perp}` of the scattering vector that is perpendicular to the cross section and depends
only on the length :math:`L`. The other factor is coplanar with the cross section and involves the component :math:`\mathbf{q}_{\|}`;
it depends on the number of sides :math:`n` and the edge length :math:`E` of the polygon.
The perpendicular factor is:

.. math::

   f_{\perp}(\mathbf{q}_{\perp}, L)
   = L \mathrm{sinc}(
     \frac{(\mathbf{q}_{\perp} \cdot \hat{\mathbf{n}})\, L}{2}
     )

where :math:`\hat{\mathbf{n}}` is the direction normal to the cross section.
The length :math:`L` gives rise to a standard sinc function for the form factor.

On the other hand, the parallel factor for a regular :math:`n`-sided polygon
of circumradius :math:`R` can be expressed as

.. math::

   f_{\|}(\mathbf{q}_{\|}, n, R)
   =
   \frac{2}{i q_{\|}^{2}}
   \sum_{j=1}^{n}
   \hat{\mathbf{n}} \cdot
   \left( \mathbf{q}_{\|} \times \mathbf{E}_{j} \right)
   \,
   \mathrm{sinc}\left( \mathbf{q}_{\|} \cdot \mathbf{E}_{j} \right)
   \exp(
     i \mathbf{q}_{\|} \cdot \mathbf{M}_{j})

In the sum over all edges, :math:`\mathbf{M}_{j}` is the vector joining the center of the polygon to the middle of the jth edge and
:math:`\mathbf{E}_{j}` is the half-edge vector.

The scattered intensity for one prism is given by:

.. math::

   I(\mathbf{q}, n, R, L)
   =
   \left| F(\mathbf{q}, n, R, L) \right|^{2}
   =
   \left|
     f_{\perp}(\mathbf{q}_{\perp}, L)
     f_{\|}(\mathbf{q}_{\|}, n, R)
   \right|^{2}

Orientation average: The 1D form factor corresponds to the orientation average with all the possible orientations having the same probability.
Instead of rotating the shape through all the possible orientations,
it is equivalent to integrate the 3D scattering vector over a sphere of radius q with the shape in its reference orientation.

The sphere is sampled using Fibonacci quadrature to provide a quasi-uniform distribution of points on the unit sphere.
The distribution of the N points is computed using the golden ratio (see fibonacci.py). 
Each point of the quadrature on the unit sphere corresponds to a vector :math:`\mathbf{u}_{j}`.
In the sum, all weights :math:`w_j` are taken identical and equal to :math:`\frac{1}{N}`.

.. math::

    P(q) =  \sum_{j=1}^{N} w_j I(q\mathbf{u}_{j}, n, R, L)

.. figure:: img/fibonacci_sphere.png

    Fibonacci sphere using N=5810 points.

Validation
----------

The model has been tested against experimental data obtained on gold nanoprisms with pentagonal cross-section (see J. Marcone et al. JAC 2025).
Moreover, comparisons with Debye formula calculations were made using DebyeCalculator library (https://github.com/FrederikLizakJohansen/DebyeCalculator).
Good agreement was found at q < 0.1 1/Angstrom.

References
----------

1. Marcone, J., Trazo, J. G., Nag, R., Goldmann, C., Ratel-Ramond, N., Hamon, C., & Impéror-Clerc, M. (2025).
   Form factor of prismatic particles for small-angle scattering analysis.
   *Journal of Applied Crystallography*, 58(2), 543‑552. https://doi.org/10.1107/S1600576725000676

2. Li, X., Shew, C., He, L., Meilleur, F., Myles, D. A. A., Liu, E., Zhang, Y., Smith, G. S.,
   Herwig, K. W., Pynn, R., & Chen, W. (2011). Scattering functions of Platonic solids.
   *Journal Of Applied Crystallography*, 44(3), 545‑557. https://doi.org/10.1107/s0021889811011691

3. Croset, B. (2017). Form factor of any polyhedron : a general compact formula and its singularities.
   *Journal Of Applied Crystallography*, 50(5), 1245‑1255. https://doi.org/10.1107/s1600576717010147

4. Wuttke, J. (2021). Numerically stable form factor of any polygon and polyhedron.
   *Journal Of Applied Crystallography*, 54(2), 580‑587. https://doi.org/10.1107/s1600576721001710

Authorship and Verification
----------------------------

* **Authors:** Marianne Imperor-Clerc (marianne.imperor@cnrs.fr)
             Jules Marcone (julesmarcone@gmail.com)
             Sara Mokhtari (smokhtari@insa-toulouse.fr)

* **Last Modified by:** MIC **Date:** 11 December 2025

* **Last Reviewed by:** SM **Date:** 17 April 2026

"""
import numpy as np
from numpy import inf

from sasmodels.special import sas_sinx_x
from sasmodels.special.fibonacci import fibonacci_sphere

name = "prism"
title = "prism of different cross-sections"
description = """
        Model for prisms of different cross-sections with orientation average using the Fibonacci quadrature"""
category = "shape:polyhedron"
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld", "1e-6/Ang^2", 126., [-inf, inf], "sld",
               "Prism scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 9.4, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["n_sides", "", 5, [3, 50], "volume",
               "Number of sides"],
              ["radius_average", "Ang", 500, [0., inf], "volume",
               "Average radius"],
              ["length", "Ang", 5000, [0., inf], "volume",
               "Length"]
               ]

### Functions for geometrical calculations: volume, surface of the cross section, average radius, edge length etc.
def form_volume(n_sides, radius_average, length):
    """
    Computes the volume of a nanoprism given its number of sides, average radius and length.
    Parameters
    ----------
    n_sides : int
        Number of sides of the regular polygon cross-section.
    radius_average : float
        Average radius of the regular polygon cross-section.
    length : float
        Length of the nanoprism.
    Returns
    -------
    volume : float
        Volume of the nanoprism.
    """
    n_sides = int(n_sides)
    edge = edge_from_gyration_radius(n_sides, radius_average)
    radius = radius_from_edge(n_sides, edge)
    surface = surface_from_radius(n_sides, radius)
    return surface * length

def edge_from_gyration_radius(n_sides:int, gyr):
    """
    Computes the edge length of an n-sided regular polygon of average radius radius_average = gyr.
    Parameters
    ----------
    n_sides : int
        Number of sides of the regular polygon.
    gyr : float
        average radius of the regular polygon.
    Returns
    -------
    edge : float
        Edge length of the regular polygon.
    """
    n_sides = int(n_sides)
    return (gyr * 2 * np.sin(np.pi/n_sides)) / np.sqrt(np.sinc(2/n_sides))

def surface_from_radius(n_sides:int, radius):
    """
    Computes the area of an n-sided regular polygon of circumradius radius.
    Parameters
    ----------
    n_sides : int
        Number of sides of the regular polygon.
    radius : float
        Circumradius of the regular polygon.
    Returns
    -------
    area : float
        Area of the regular polygon.
    """
    return n_sides * (radius**2) * np.sin(np.pi/n_sides) * np.cos(np.pi/n_sides)

def radius_from_edge(n_sides:int, edge):
    """
    Computes the circumradius of an n-sided regular polygon of edge length edge.
    Parameters
    ----------
    n_sides : int
        Number of sides of the regular polygon.
    edge : float
        Edge length of the regular polygon.
    Returns
    -------
    radius : float
        Circumradius of the regular polygon.
    """
    return edge/(2*np.sin(np.pi/n_sides))

def shape_generator(n_sides, radius):
    """
    Computes the list of vertices of a n-sided regular polygon, of circumradius "radius".
    Parameters
    ----------
    number_of_sides : int
        Number of sides of the regular polygon.
    radius : float
        Circumradius of the regular polygon.
    Returns
    -------
    vertices : np.ndarray shape (2, n)
        List of the vertices of the regular polygon.
    """
    theta = np.linspace(0, 2 * np.pi, n_sides, endpoint=False)
    return radius*np.vstack((np.cos(theta), np.sin(theta)))  # (2 x n)

def edgecenters_generator(vertices:np.ndarray):
    """
    Computes the list of edge centers of a 2D shape defined by its vertices.
    Parameters
    ----------
    vertices : np.ndarray shape (2, n)
        list of the vertices of the 2D-shape each listed as a list of the 2D coordinates of the shape in the plane
    Returns
    -------
    edgecenter : np.ndarray shape (2, n)
        List of the edge centers of the 2D shape.
    Preconditions
    -------------
    a loop, all in the plane
    """
    return (vertices + np.roll(vertices, -1, axis=1))/2

def halfedges_generator(vertices:np.ndarray):
    """
    Computes the list of half-edges of a 2D shape defined by its vertices.
    Parameters
    ----------
    vertices : np.ndarray shape (2, n)
        list of the vertices of the 2D-shape each listed as a list of the 2D coordinates of the shape in the plane
    Returns
    -------
    halfedge : np.ndarray shape (2, n)
        List of the half-edges of the 2D shape.
    Preconditions
    -------------
    a loop, all in the plane
    """
    return (np.roll(vertices, -1, axis=1) - vertices)/2

# Form factor and intensity calculations
# Reminder: the form factor is defined as: parallel factor (= complex function) * perpendicular factor (= sinc function)
# see documentation and Reference [1] (Jules Marcone et al. "Form factor of prismatic particles for small-angle scattering analysis")

def parallel_factor(vertices:np.ndarray, q:np.ndarray, c:float): # gives the area of the polygon at q==[0,0]
    """
    Computes the parallel form factor of a 2D shape defined by its vertices at a specific in-plane scattering vector q.
    Parameters
    ----------
    vertices : np.ndarray shape (2, n) 
        list of the vertices of the 2D-shape each listed as a list of the
        2D coordinates of the shape in the plane
    q : np.ndarray shape (N, 2)
        listed as a list of coordinates in the plane (q is q//)
    c : float
        an arbitrary constant
    Returns
    -------
    parallel_factor : np.ndarray of complex, shape (N,)
        Parallel form factor of the 2D shape at the specific in-plane scattering vector q.
    Preconditions
    -------------
    a loop, all in the plane
    """
    qmodulus2 = q[:, 0]**2 + q[:, 1]**2
    cutoff = 10**-30
    qmodulus2_cutoff = np.where(qmodulus2==0, cutoff, qmodulus2) # replace by cutoff if equal to 0
    # This case starts occuring for computations of the prism formfactor, where during the
    # calculation of the orientational average, q// may be equal to 0.
    # if qmodulus2==0:
    #    return 0
    edgecenters = edgecenters_generator(vertices)
    halfedges = halfedges_generator(vertices)
    sum = 0
    # for i in range(len(vertices)):
    for i in range(vertices.shape[1]):
        qEj = np.dot(q, halfedges[:,i]) # scalar_product(q, halfedges[i])
        # triple_product = q[0]*halfedges[i][1]-q[1]*halfedges[i][0]
        triple_product = q[:, 0]*halfedges[1, i]  - q[:, 1]*halfedges[0, i]
        #The exp(iqRj) is rewritten as a sum of cos+isin because of the way math.exp() works (not allowing complex as input)
        qRj = np.dot(q, edgecenters[:,i])
        sum += triple_product * (sas_sinx_x(qEj)*(np.cos(qRj)+np.sin(qRj)*1J)-c)
    return (2/(1J*qmodulus2_cutoff)*sum)

def Fqabc(qa, qb, qc, n_sides, radius_average, length): # Form factor in 3D of the nanoprism for (qa, qb, qc) scattering vector
    """
    Computes the form factor amplitude of a prism with a n-sided regular polygon cross-section and length for (qa, qb, qc) scattering vector.
    Takes in account the parallel factor (complex function) and the perpendicular factor (sinc function).
    Parameters
    ----------
    qa, qb, qc : np.ndarray shape (N,)
        components of the scattering vector q
    n_sides : int
        Number of sides of the regular polygon cross-section.
    radius_average : float
        Average radius.
    length : float
        Length of the nanoprism.
    Returns
    -------
    Fqabc : np.ndarray of complex, shape (N,)
        Form factor amplitude of the nanoprism at the specific three dimensional q.
    """
    qab = np.vstack((qa, qb)).T
    edge = edge_from_gyration_radius(n_sides, radius_average)
    radius = radius_from_edge(n_sides,edge)
    vertices = shape_generator(n_sides,radius)
    perpendicular_factor_value = sas_sinx_x(qc*length/2) * length
    parallel_factor_value = parallel_factor(vertices,qab,0.)
    A = parallel_factor_value * perpendicular_factor_value
    return A

def Iqabc(qa,qb,qc,n_sides,radius_average, length): # proportionnal to the volume**2
    """
    Calls the function that computes the edge length and the scattered intensity.
    Parameters
    ----------
    qa, qb, qc :  np.ndarray shape (N,)
        components of the scattering vector q
    n_sides : int
        Number of sides of the regular polygon cross-section.
    radius_average : float
        Average radius of the regular polygon cross-section.
    length : float
        Length of the nanoprism.
    Returns
    -------
    Iqabc : np.ndarray of float, shape (N,)
        Scattered intensity of the nanoprism at the specific three dimensional q components.
    """
    n_sides = int(n_sides)
    A = Fqabc(qa, qb, qc, n_sides, radius_average, length)
    intensity = (np.abs(A))**2  # intensity is proportional to volume
    return intensity

def Iq(q, sld, sld_solvent, n_sides:int, radius_average:float, length:float, npoints_fibonacci:int= 500):
    """
    Computes the scattering intensity I(q) of nanoprisms averaged over all orientations using the Fibonacci quadrature.
    The number of points on the sphere is set by npoints_fibonacci. Each point has an equal weight = 1/npoints_fibonacci.
    Parameters
    ----------
    q : float ou array
        Norm of the scattering vector
    sld, sld_solvent :
        Contrast of scattering length density
    n_sides, radius_average, length :
        Geometrical parameters of the prism
    npoints_fibonacci : int
        Number of Fibonacci points on the sphere, set to 500 by default
        (higher number increases accuracy but also computation time, 500 is usually a good compromise)
    Returns
    -------
    Iq : ndarray
        Scattering intensity averaged over all orientations
    """
    n_sides = int(n_sides)
    q = np.atleast_1d(q)
    q_unit,w = fibonacci_sphere(npoints_fibonacci)   # shape (npoints, 3)
    # Projections
    qa = q[:, np.newaxis] * q_unit[:, 0][np.newaxis, :]
    qb = q[:, np.newaxis] * q_unit[:, 1][np.newaxis, :]
    qc = q[:, np.newaxis] * q_unit[:, 2][np.newaxis, :]
    # # Compute intensity
    intensity = Iqabc(qa.ravel(), qb.ravel(), qc.ravel(), n_sides, radius_average, length).reshape(qa.shape)
    # Uniform average over the sphere
    integral = np.mean(intensity, axis=1)
    return (integral) * (sld - sld_solvent)**2 * 10**-4

Iq.vectorized = True

tests = [
    [{"background": 0, "scale": 1, "n_sides": 4, "radius_average": 10, "length": 200, "sld": 1., "sld_solvent": 0.},
     0.01, 5.62789],
    [{"background": 0, "scale": 1, "n_sides": 4, "radius_average": 10, "length": 200, "sld": 1., "sld_solvent": 0.},
     [0.01, 0.1], [5.62789, 0.73696]],
]
