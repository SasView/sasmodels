r"""
This model provides the form factor for nanoprisms with different cross-sections.
Orientation averaging is done by using the Fibonacci quadrature.
This quadrature provides a quasi-uniform distribution of points on the unit sphere
using the golden ratio. The number of points to generate on the unit sphere is set to 500 points, it usually provides
a good balance between accuracy and computational efficiency.

Definition
----------

We consider particles having the shape of a right prism of length L and a cross section made of a regular polygon with n sides as illustrated in the figure below.

.. figure:: img/nanoprisms_geometry.jpg

The size of a regular polygon can be characterized by its edge length E or by R, the circumradius of the polygon shown in black in the figure, with :

.. math::

    R = \frac{E}{2\sin(\pi/n)}


The apothem is the radius of the dashed circle, it is defined as:

.. math::

    R\cos(\pi/n)

For comparison purposes, it is convenient to introduce an average radius R_ave (shown in red in the figure).
The area of the n-sided regular polygon is given by :

.. math::

   A = \pi R^{2} \, \mathrm{sinc}\left(\frac{2\pi}{n}\right)
     = \pi R_{\mathrm{ave}}^{2}

and the volume of the nanoprism is therefore given by:

.. math::

    V = L \times {A}

where R_ave is the radius of the equivalent disc having the same area as the n-sided polygon.
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
   = \mathrm{sinc}(
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

In the sum over all edges, :math:`\mathbf{M}_{j}` is the vector joining the center of the polygon to the middle of the jth edge andand
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

Orientation averaging: The 1D form factor corresponds to the orientation average with all the possible orientations having the same probability.
Instead of rotating the shape through all the possible orientations in the integral,
it is equivalent to integrate the 3D scattering vector over a sphere of radius q with the shape in its reference orientation.

.. math::

    P(q) =  \frac{2}{\pi} \int_0^{\frac{\pi}{2}} \,
    \int_0^{\frac{\pi}{2}} A_P^2(q) \, \sin\theta \, d\theta \, d\phi

The sphere is sampled using Fibonacci quadrature to provide a quasi-uniform distribution of points on the unit sphere.
The repartition of the points is computed using the golden ratio (see fibonacci.py).


.. figure:: img/fibonacci_sphere.png

    Fibonacci sphere using 5810 points.

Validation
----------

The model has been tested against experimental data obtained on gold nanoprisms with pentagonal cross-section (see J. Marcone et al. JAC 2025).
Moreover, comparisons with Debye formula calculations were made using DebyeCalculator library (https://github.com/FrederikLizakJohansen/DebyeCalculator).
Good agreement was found at q < 0.1 1/Angstrom.

References
----------

1. Jules Marcone et al. "Form factor of prismatic particles for small-angle scattering analysis".
   In:  J. Appl. Cryst. (2025) 58, 543-552
   DOI: https://doi.org/10.1107/S1600576725000676

2. Wei-Ren Chen et al. "Scattering functions of Platonic solids".
   In: Journal of Applied Crystallography - J. Appl. Cryst. (June 2011).
   DOI:  https://doi.org/10.1107/S0021889811011691

3. Croset, Bernard, "Form factor of any polyhedron: a general compact
   formula and its singularities" In: J. Appl. Cryst. (2017). 50, 1245–1255
   https://doi.org/10.1107/S1600576717010147

4. Wuttke, J. Numerically stable form factor of any polygon and polyhedron
   J Appl Cryst 54, 580-587 (2021)
   https://doi.org/10.1107/S160057672100171

Authorship and Verification
----------------------------

* **Authors:** Marianne Imperor-Clerc (marianne.imperor@cnrs.fr)
             Jules Marcone (jules.marcone@universite-paris-saclay.fr)
             Sara Mokhtari (smokhtari@insa-toulouse.fr)

* **Last Modified by:** MIC **Date:** 11 December 2025

* **Last Reviewed by:** SM **Date:** 16 January 2026

"""
import numpy as np
from numpy import inf

from sasmodels.quadratures.fibonacci import fibonacci_sphere

name = "nanoprisms"
title = "nanoprisms of different cross-sections"
description = """
        Model for nanoprisms of different cross-sections with orientation averaging using the Fibonacci quadrature"""
category = "plugin"
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld", "1e-6/Ang^2", 126., [-inf, inf], "sld",
               "nanoprism scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 9.4, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["nsides", "", 5, [3, 50], "volume",
               "nsides"],
              ["Rave", "Ang", 500, [0., inf], "volume",
               "Average radius"],
              ["L", "Ang", 5000, [0., inf], "volume",
               "length"]
               ]

### Functions for geometrical calculations: volume, surface of the cross section, average radius, edge length etc.
def form_volume(nsides, Rave, L):
    """
    Computes the volume of a nanoprism given its number of sides, average radius and length.
    Parameters
    ----------
    nsides : int
        Number of sides of the regular polygon cross-section.
    Rave : float
        Average radius of the regular polygon cross-section.
    L : float
        Length of the nanoprism.
    Returns
    -------
    volume : float
        Volume of the nanoprism.
    """
    nsides = int(nsides)
    edge = edge_from_gyration_radius(nsides,Rave)
    return volume_nanoprism(nsides,edge,L)

def edge_from_gyration_radius(nsides:int, gyr):
    """
    Computes the edge length of an n-sided regular polygon of gyration ratio gyr.
    Parameters
    ----------
    nsides : int
        Number of sides of the regular polygon.
    gyr : float
        Gyration radius of the regular polygon.
    Returns
    -------
    edge : float
        Edge length of the regular polygon.
    """
    nsides = int(nsides)
    return (gyr * 2 * np.sin(np.pi/nsides)) / np.sqrt(np.sinc(2/nsides))

def volume_nanoprism(nsides:int, edge, L):
    """ Computes the volume of a nanoprism with a cross-section
    in the shape of an n-sided regular polygon of edge length "edge", and length of the nanoprism "L".
    Parameters
    ----------
    nsides : int
        Number of sides of the regular polygon cross-section.
    edge : float
        Edge length of the regular polygon cross-section.
    L : float
        Length of the nanoprism.
    Returns
    -------
    volume : float
        Volume of the nanoprism.
    """
    radius = radius_from_edge(nsides, edge)
    surface = surface_from_radius(nsides, radius)
    return surface * L

def surface_from_radius(nsides:int, radius):
    """
    Computes the area of an n-sided regular polygon of circumradius radius.
    Parameters
    ----------
    nsides : int
        Number of sides of the regular polygon.
    radius : float
        Circumradius of the regular polygon.
    Returns
    -------
    area : float
        Area of the regular polygon.
    """
    return nsides * (radius**2) * np.sin(np.pi/nsides) * np.cos(np.pi/nsides)

def radius_from_edge(nsides:int, edge):
    """
    Computes the circumradius of an n-sided regular polygon of edge length edge.
    Parameters
    ----------
    nsides : int
        Number of sides of the regular polygon.
    edge : float
        Edge length of the regular polygon.
    Returns
    -------
    radius : float
        Circumradius of the regular polygon.
    """
    return edge/(2*np.sin(np.pi/nsides))

def shape_generator(number_of_sides, radius):
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
    vertices : list
        List of the vertices of the regular polygon.
    """
    vertices = [[radius*np.cos(2*step*np.pi/number_of_sides), radius*np.sin(2*step*np.pi/number_of_sides)] for step in range(0,number_of_sides)]
    return vertices

def edgecenters_generator(vertices:list):
    """
    Computes the list of edge centers of a 2D shape defined by its vertices.
    Parameters
    ----------
    vertices : list
        list of the vertices of the 2D-shape each listed as a list of the 2D coordinates of the shape in the plane
    Returns
    -------
    edgecenter : list
        List of the edge centers of the 2D shape.
    Preconditions
    -------------
    a loop, all in the plane
    """
    extended_vertices = [] + vertices
    extended_vertices.append(vertices[0])
    edgecenter = []
    for i in range(len(vertices)):
        coordinates=[]
        for j in range(2):
            coordinates.append((extended_vertices[i+1][j]+extended_vertices[i][j])/2)
        edgecenter.append(coordinates)
    return edgecenter

def halfedges_generator(vertices:list):
    """
    Computes the list of half-edges of a 2D shape defined by its vertices.
    Parameters
    ----------
    vertices : list
        list of the vertices of the 2D-shape each listed as a list of the 2D coordinates of the shape in the plane
    Returns
    -------
    halfedge : list
        List of the half-edges of the 2D shape.
    Preconditions
    -------------
    a loop, all in the plane
    """
    extended_vertices = [] + vertices
    extended_vertices.append(vertices[0])
    halfedge = []
    for i in range(len(vertices)):
        coordinates=[]
        for j in range(2):
            coordinates.append((extended_vertices[i+1][j]-extended_vertices[i][j])/2)
        halfedge.append(coordinates)
    return halfedge

# useful function sinc(x) = sin(x)/x
def sinc(x):
    """
    Computes the cardinal sinus of x since in numpy sinc is defined with pi*x for the argument.
    Parameters
    ----------
    x : float or array
        Input value(s).
    Returns
    -------
    sinc_x : float or array
        Cardinal sinus of x.
    """
    return np.sinc(x/np.pi)

# useful functions for scalar product
def scalar_product(u,v):
    """
    Computes the scalar product of two vectors u and v since numpy.dot does not work with different shapes of arrays.
    Parameters
    ----------
    u and v : list
        lists of the coordinates of the vectors.
    Returns
    -------
    scalar_product : float
        Scalar product of the two vectors.
    """
    somme = 0
    for i in range(len(u)):
        somme += u[i]*v[i]
    return somme

# Form factor and intensity calculation
# Reminder: the form factor is defined as: parallel factor (= complex function) * perpendicular factor (= sinc function)
# see documentation and Reference [1] (Jules Marcone et al. "Form factor of prismatic particles for small-angle scattering analysis")
def parallel_factor(vertices:list, q, c): # gives the area of the polygon at q==[0,0]
    """
    Computes the parallel form factor of a 2D shape defined by its vertices at a specific in-plane scattering vector q.
    Parameters
    ----------
    vertices : list
        list of the vertices of the 2D-shape each listed as a list of the
        2D coordinates of the shape in the plane
    q : list
        listed as a list of coordinates in the plane (q is q//)
    c : float
        an arbitrary constant
    Returns
    -------
    parallel_factor : complex
        Parallel form factor of the 2D shape at the specific in-plane scattering vector q.
    Preconditions
    -------------
    a loop, all in the plane
    """
    qmodulus2 = q[0]**2+q[1]**2
    qmodulus2 = np.asanyarray(qmodulus2) #conversion to array
    cutoff = 10**-30
    qmodulus2_cutoff = np.where(qmodulus2==0, cutoff, qmodulus2) # replace by cutoff if equal to 0
    #This case starts occuring for computations of the prism formfactor, where during the
    #calculation of the orientational average, q// may be equal to 0.
    #if qmodulus2==0:
    #    return 0
    edgecenters = edgecenters_generator(vertices)
    halfedges = halfedges_generator(vertices)
    sum = 0
    for i in range(len(vertices)):
        qEj = scalar_product(q, halfedges[i]) # scalar_product(q, halfedges[i])
        triple_product = q[0]*halfedges[i][1]-q[1]*halfedges[i][0]
        #The exp(iqRj) is rewritten as a sum of cos+isin because of the way math.exp() works (not allowing complex as input)
        qRj = scalar_product(q, edgecenters[i])
        sum += triple_product * (sinc(qEj)*(np.cos(qRj)+np.sin(qRj)*1J)-c)
    return (2/(1J*qmodulus2_cutoff)*sum)

def Amp_nanoprism(q,nsides,edge,L): # From factor in 3D of the nanoprism q with three components
    """
    Computes the form factor amplitude of a prism with a n-sided regular polygon cross-section and length L for a specific three dimensional q.
    Takes in account the parallel factor (complex function) and the perpendicular factor (sinc function).
    Parameters
    ----------
    q : list
        listed as a list of three components [qa,qb,qc]
    nsides : int
        Number of sides of the regular polygon cross-section.
    edge : float
        Edge length of the regular polygon cross-section.
    L : float
        Length of the nanoprism.
    Returns
    -------
    Amp_nanoprism : complex
        Form factor amplitude of the nanoprism at the specific three dimensional q.
    """
    qa, qb, qc = q[0], q[1], q[2]
    qab = [qa,qb]
    radius = radius_from_edge(nsides,edge)
    vertices = shape_generator(nsides,radius)
#   A=parallel_factor(vertices,qab,0.)*sinc(qc*L/2)/surface_from_radius(nsides,radius)
    perpendicular_factor = sinc(qc*L/2)
    A = parallel_factor(vertices,qab,0.)*perpendicular_factor # parallel form factor * perpendicular form factor
    return A

def Iqabc(qa,qb,qc,nsides,Rave,L): # proportionnal to the volume**2
    """
    Calls the function that computes the edge length and the scattered intensity.
    Parameters
    ----------
    qa, qb, qc : float or array
        components of the scattering vector q
    nsides : int
        Number of sides of the regular polygon cross-section.
    Rave : float
        Average radius of the regular polygon cross-section.
    L : float
        Length of the nanoprism.
    Returns
    -------
    Iqabc : float or array
        Scattered intensity of the nanoprism at the specific three dimensional q components.
    """
    nsides=int(nsides)
    edge = edge_from_gyration_radius(nsides, Rave)
    # intensity = I_nanoprism([qa,qb,qc], nsides, edge, L)
    A = Amp_nanoprism([qa,qb,qc], nsides, edge, L)
    intensity = (np.abs(A))**2
    intensity = intensity * (L)**2  # multiplication by the volume at the power 2 (? volume or L only ?)
    return intensity

def Iq(q, sld, sld_solvent, nsides:int, Rave, L, npoints_fibonacci:int= 500):
    """
    Computes the scattering intensity I(q) of nanoprisms averaged over all orientations using the Fibonacci quadrature.
    The number of points on the sphere is set by npoints_fibonacci. Each point has an equal weight = 1/npoints_fibonacci.
    Parameters
    ----------
    q : float ou array
        Norm of the scattering vector
    sld, sld_solvent :
        Contrast of scattering length density
    nsides, Rave, L :
        Geometrical parameters of the prism
    npoints_fibonacci : int
        Number of Fibonacci points on the sphere, set to 500 by default
        (higher number increases accuracy but also computation time, 500 is usually a good compromise)
    Returns
    -------
    Iq : ndarray
        Scattering intensity averaged over all orientations
    """
    nsides = int(nsides)
    q = np.atleast_1d(q)
    q_unit,w = fibonacci_sphere(npoints_fibonacci)   # shape (npoints, 3)
    # Projections
    qa = q[:, np.newaxis] * q_unit[:, 0][np.newaxis, :]
    qb = q[:, np.newaxis] * q_unit[:, 1][np.newaxis, :]
    qc = q[:, np.newaxis] * q_unit[:, 2][np.newaxis, :]
    # Compute intensity
    intensity = Iqabc(qa, qb, qc, nsides, Rave, L)  # shape (nq, npoints)
    # Uniform average over the sphere
    integral = np.sum(w[np.newaxis, :] * intensity, axis=1)
    integral = np.mean(intensity, axis=1)
    return integral * (sld - sld_solvent)**2  * 1e4 # Convert to [cm-1]

Iq.vectorized = True

tests = [
    [{"background": 0, "scale": 1, "nsides": 4, "Rave": 10, "L":200,"sld": 1., "sld_solvent": 0.},
     0.01, 5.62789],
    [{"background": 0, "scale": 1, "nsides": 4, "Rave": 10, "L":200,"sld": 1., "sld_solvent": 0.},
     [0.01, 0.1], [5.62789, 0.73696]],
]
