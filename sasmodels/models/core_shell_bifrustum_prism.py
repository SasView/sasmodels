r"""
This model provides the form factor for core-shell particles with a bifrustum core and a prism shell.
Orientation averaging is done by using the Fibonacci quadrature.
This quadrature provides a quasi-uniform distribution of points on the unit sphere
using the golden ratio. The number of points to generate on the unit sphere is set to 500 points, it usually provides
a good balance between accuracy and computational efficiency.

Definition
----------
**Geometry.** This model describes a core-shell particle with a bifrustum core and a prism shell. Both shapes are defined with their symmetry axis along the z-axis and centered at the origin. The bifrustum core can be rotated around the z-axis by an angle :math:`\alpha` (*core_rotation*)
and translated along the z-axis by a distance :math:`z_0` (*core_translation*) with respect to the prism shell.

.. figure:: img/cs5sides.png

    Core-shell bifrustum-prism geometry.

Both the bifrustum and prism are already existing models in SasView. For further details on the bifrustum and prism models,
please refer to their respective model documentations (see bifrustum.py and prism.py).

**Parameters notation.** The slds are called *sld_core*, *sld_shell* and *sld_solvent* for the bifrustum core, prism shell and solvent, respectively.
The geometric parameters of the bifrustum core and prism shell are prefixed with *core_* and *shell_*, respectively, to distinguish them.

**Core parameters (bifrustum).** The bifrustum is a truncated bipyramid with a polygonal cross-section. It is defined
by  :math:`n` the number of sides of the polygonal cross-section (core_n_sides),
:math:`R` the circumradius of the polygonal cross-section (core_radius_polygon),
:math:`H` the half-height of the bifrustum before truncation (core_half_height)
and :math:`t` the truncation ratio (core_truncation) with :math:`t = 0` giving
a sharp bipyramid and :math:`t = 1` a fully truncated bipyramid to 0 height.

**Shell parameters (prism).** The prism is a rod with a polygonal cross-section.
It is defined by :math:`n` the number of sides of the polygonal cross-section (shell_n_sides),
:math:`R` the average radius of the polygonal cross-section (shell_radius_average)
and :math:`L` the length of the prism (shell_length_prism).

.. figure:: img/core_shell_bifrustum_prism_parameters.png

    Prism and and bifrustum geometries and their parameters.

**Form factor.** The prism and bifrustum form factor expressions are described in their respective model documentations (see bifrustum.py and prism.py).
The total scattering amplitude is given by the sum of the prism and bifrustum form factors,
weighted by their respective scattering length densities (sld_shell, sld_core)
and the solvent scattering length density (sld_solvent).

The general expression for the total scattering amplitude :math:`F_{\text{cs}}(\mathbf{q})` is given by:

.. math::

    F_{\text{cs}}(\mathbf{q, \alpha, z_0}) = F_{\text{shell}}(\mathbf{q}) \cdot (\text{sld}_{\text{shell}} - \text{sld}_{\text{solvent}}) +
    F_{\text{core}}(\mathbf{\alpha, q}) \cdot e^{i \mathbf{q} \cdot \mathbf{z_0}} \cdot (\text{sld}_{\text{core}} - \text{sld}_{\text{shell}})

with :math:`F_{\text{shell}}(\mathbf{q})` and :math:`F_{\text{core}}(\mathbf{q})` being the form factor amplitudes of the prism and bifrustum,
respectively, :math:`\mathbf{z_0}` being the translation vector of the bifrustum core along the z-axis with respect to the prism shell and
:math:`\alpha` being the rotation matrix of the bifrustum core around the z-axis with respect to the prism shell.


**Orientation average**: The 1D form factor corresponds to the orientation average with all the possible orientations having the same probability.
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

F(q) checked against an independent method (Debye equation with Debye calculator: https://github.com/FrederikLizakJohansen/DebyeCalculator)
for different prism cross-sections. The Debye equation is based on atomic positions while SasView model is based on analytical expressions.
Good agreement was found at q < 0.1 1/Angstrom. The model has not been tested with experimental data yet.

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

* **Authors:** Chi-Huan Tung (tungc@ornl.gov)
             Marianne Imperor-Clerc (marianne.imperor@cnrs.fr)
             Sara Mokhtari (smokhtari@insa-toulouse.fr)

* **Last Modified by:** SM **Date:** 17 April 2026

* **Last Reviewed by:** MIC **Date:** 15 April 2026

"""
import numpy as np
from numpy import inf

from sasmodels.models import bifrustum as core_model
from sasmodels.models import prism as shell_model
from sasmodels.special.fibonacci import fibonacci_sphere

name = "core_shell_bifrustum_prism"
title = "Core-Shell Bifrustum Prism"
description = """
        Model for core-shell particles with a bifrustum core and a prism shell"""
category = "shape:polyhedron"


def _geometry(model):
    """
    Returns the list of geometry parameters for a model.
    No sdl and angles.
    """
    return  [p for p in model.parameters if p[4] == "volume"]

def _prefixed(prefix, pars):
    return [[prefix + p[0]] + list(p[1:]) for p in pars]

CORE, SHELL = _geometry(core_model), _geometry(shell_model)
N_CORE, N_SHELL = len(CORE), len(SHELL)


#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = (
    ###########################   SLD PARAMETERS    #########################
    [["sld_core", "1e-6/Ang^2", 126., [-inf, inf], "sld", "Core scattering length density"],
     ["sld_shell", "1e-6/Ang^2", 126., [-inf, inf], "sld", "Shell scattering length density"],
     ["sld_solvent", "1e-6/Ang^2", 9.4, [-inf, inf], "sld", "Solvent scattering length density"]]
     ####### Geometry parameters for the prism shell and bifrustum core #######
    + _prefixed("shell_", SHELL) + _prefixed("core_", CORE)
    ##################  Seed orientation parameters    ##################
    +[["core_rotation", "deg", 0.0, [-180., 180.], "",
     "Rotation angle of the bifrustum core around the z-axis (vertical axis) with respect to the prism shell"],
    ["core_translation", "Ang", 0.0, [-inf, inf], "",
     "Translation of the bifrustum core along the z-axis (vertical axis) with respect to the prism shell"],
    ################## orientation angles of the whole system ###########################
    ["theta", "degrees", 0, [-180, 180], "orientation", "c axis to beam angle"],
    ["phi", "degrees", 0, [-360, 360], "orientation", "rotation about beam"],
    ["psi", "degrees", 0, [-360, 360], "orientation", "rotation about c axis"]])


def _split(args):
    """Kernel arguments (in the order of `parameters`) -> shell geometry, core geometry, rotation, z0."""
    shell = args[:N_SHELL]
    core = args[N_SHELL:N_SHELL + N_CORE]
    core_rotation, core_translation = args[N_SHELL + N_CORE:]
    return shell, core, core_rotation, core_translation

def form_volume(*args):  # using the volume of the prism shell
    return shell_model.form_volume(*args[:N_SHELL])

def Fqabc(qa, qb, qc, sld_core, sld_shell, sld_solvent, *args):
    """
    Computes the form factor amplitude of the core-shell particle
    with a bifrustum core and a prism shell from their geomtric parameters,
    their slds and the three components of the scattering vector q.
    Parameters
    ----------
    qa, qb, qc : float or array
        Components of the scattering vector q
    sld_core, sld_shell, sld_solvent : float
        Scattering length densities of the core, shell and solvent
    args : list
        List of parameters in the order of `parameters` (see above)
    Returns
    -------
    float or array
        Form factor amplitude of the core-shell particle

    """
    # 0. Retreive geometric parameters of the prism shell and bifrustum core, r
    # and the rotation and translation of the bifrustum core
    shell, core, core_rotation, core_translation = _split(args)

    # 1. Shell form factor amplitude
    # Issue with Fqabc of prism, it expects flattened q values?
    # TO DO: if we aim to make a general core-shell model, we should modify the prism model to accept
    # 3D q values without flattening
    F_shell = shell_model.Fqabc(qa, qb, qc, *shell)

    # 2. Core form factor amplitude with rotation and translation
    # 2.1. Rotation of the bifrustum core: rotation around z-axis by core_rotation
    if core_rotation != 0.0:
        rot = np.radians(core_rotation)
        qa_core = qa * np.cos(rot) - qb * np.sin(rot)
        qb_core = qa * np.sin(rot) + qb * np.cos(rot)
        qc_core = qc # invariant
    else:
        qa_core, qb_core, qc_core = qa, qb, qc
    F_core = core_model.Fqabc(qa_core, qb_core, qc_core, *core)
    # 2.2. Phase if translation along z-axis is applied
    phase = np.exp(1j * qc * core_translation) # phase factor for translation along z-axis, if z = 0, then phase = 1

    # 3. Total scattering amplitude: prism + bifrustum with phase factor
    total_scattering_amplitude = F_shell * (sld_shell - sld_solvent) + (F_core * phase) * (sld_core - sld_shell)

    return total_scattering_amplitude

def _Iqabc(qa, qb, qc, sld_core, sld_shell, sld_solvent, *args):
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

    A = Fqabc(qa, qb, qc, sld_core, sld_shell, sld_solvent, *args)  # shape (nq, npoints)
    intensity = (np.abs(A))**2  # (1e-6/Ang^2)**2 = 1e-12/Ang^4
    return intensity

def Iqabc(qa, qb, qc, sld_core, sld_shell, sld_solvent, *args):
    """
    Computes the scattered intensity of the tetrahedron at the specific
    three dimensional q components from the form factor amplitude Fqabc.
    Parameters
    ----------
    qa, qb, qc : float or array
        Components of the scattering vector q
n.
    Returns
    -------
    Iqabc : float or array
        Scattered intensity of the tetrahedron at the specific three
        dimensional q components.
    """

    intensity = _Iqabc(qa, qb, qc, sld_core, sld_shell, sld_solvent, *args)

    return intensity * 10**-4

def Iq(q, sld_core, sld_shell, sld_solvent, *args, npoints_fibonacci: int = 500):
    """
    Computes the scattering intensity I(q) of tetrahedrons averaged over all
    orientations using the Fibonacci quadrature.
    The number of points on the sphere is set by npoints_fibonacci.
    Each point has an equal weight = 1/npoints_fibonacci.
    Parameters
    ----------
    q : float or array
        Norm of the scattering vector
    sld_core, sld_shell, sld_solvent :
        Scattering length densities of the core, shell and solvent
    args : list
        List of parameters in the order of `parameters` (see above)
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
    intensity = _Iqabc(qa, qb, qc, sld_core, sld_shell, sld_solvent, *args)  # shape (nq, npoints)
    # Uniform average over the sphere
    integral = np.mean(intensity, axis=1)

    return (integral) * 10**-4


Iq.vectorized = True

tests = [
    # core centred, core more contrasted than the shell
    [{"background": 0, "scale": 1, "sld_core": 2., "sld_shell": 1., "sld_solvent": 0.,
      "shell_n_sides": 4, "shell_radius_average": 10, "shell_length": 100,
      "core_n_sides": 4, "core_radius_polygon": 10, "core_half_height": 80, "core_truncation": 0.5},
     [0.01, 0.1], [5.15199, 1.32571]],
    # same particle with the core rotated by 45 deg and translated by 10 Ang along z
    [{"background": 0, "scale": 1, "sld_core": 2., "sld_shell": 1., "sld_solvent": 0.,
      "shell_n_sides": 4, "shell_radius_average": 10, "shell_length": 100,
      "core_n_sides": 4, "core_radius_polygon": 10, "core_half_height": 80, "core_truncation": 0.5,
      "core_rotation": 45, "core_translation": 10},
     [0.01, 0.1], [5.14899, 1.31752]],
]
