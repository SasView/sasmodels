r"""
Definition
----------

Calculates the form factor for a rectangular solid with a core-shell structure.
The thickness and the scattering length density of the shell or "rim" can be
different on each (pair) of faces. The three dimensions of the core of the
parallelepiped (strictly here a cuboid) may be given in *any* size order as
long as the particles are randomly oriented (i.e. take on all possible
orientations see notes on 2D below). To avoid multiple fit solutions,
especially with Monte-Carlo fit methods, it may be advisable to restrict their
ranges. There may be a number of closely similar "best fits", so some trial and
error, or fixing of some dimensions at expected values, may help.

The form factor is normalized by the particle volume $V$ such that

.. math::

    I(q) = \frac{\text{scale}}{V} \langle P(q,\alpha,\beta) \rangle
    + \text{background}

where $\langle \ldots \rangle$ is an average over all possible orientations
of the rectangular solid, and the usual $\Delta \rho^2 \ V^2$ term cannot be
pulled out of the form factor term due to the multiple slds in the model.

The core of the solid is defined by the dimensions $A$, $B$, $C$ here shown
such that $A < B < C$.

.. figure:: img/parallelepiped_geometry.jpg

   Core of the core shell parallelepiped with the corresponding definition
   of sides.


There are rectangular "slabs" of thickness $t_A$ that add to the $A$ dimension
(on the $BC$ faces). There are similar slabs on the $AC$ $(=t_B)$ and $AB$
$(=t_C)$ faces. The projection in the $AB$ plane is

.. figure:: img/core_shell_parallelepiped_projection.jpg

   AB cut through the core-shell parallelipiped showing the cross secion of
   four of the six shell slabs. As can be seen, this model leaves **"gaps"**
   at the corners of the solid.


The total volume of the solid is thus given as

.. math::

    V = ABC + 2t_ABC + 2t_BAC + 2t_CAB

The intensity calculated follows the :ref:`parallelepiped` model, with the
core-shell intensity being calculated as the square of the sum of the
amplitudes of the core and the slabs on the edges. The scattering amplitude is
computed for a particular orientation of the core-shell parallelepiped with
respect to the scattering vector and then averaged over all possible
orientations, where $\alpha$ is the angle between the $z$ axis and the $C$ axis
of the parallelepiped, and $\beta$ is the angle between the projection of the
particle in the $xy$ detector plane and the $y$ axis.

.. math::

    P(q)=\frac {\int_{0}^{\pi/2}\int_{0}^{\pi/2}F^2(q,\alpha,\beta) \ sin\alpha
    \ d\alpha \ d\beta} {\int_{0}^{\pi/2} \ sin\alpha \ d\alpha \ d\beta}

and

.. math::

    F(q,\alpha,\beta)
    &= (\rho_\text{core}-\rho_\text{solvent})
       S(Q_A, A) S(Q_B, B) S(Q_C, C) \\
    &+ (\rho_\text{A}-\rho_\text{solvent})
        \left[S(Q_A, A+2t_A) - S(Q_A, A)\right] S(Q_B, B) S(Q_C, C) \\
    &+ (\rho_\text{B}-\rho_\text{solvent})
        S(Q_A, A) \left[S(Q_B, B+2t_B) - S(Q_B, B)\right] S(Q_C, C) \\
    &+ (\rho_\text{C}-\rho_\text{solvent})
        S(Q_A, A) S(Q_B, B) \left[S(Q_C, C+2t_C) - S(Q_C, C)\right]

with

.. math::

    S(Q_X, L) = L \frac{\sin (\tfrac{1}{2} Q_X L)}{\tfrac{1}{2} Q_X L}

and

.. math::

    Q_A &= q \sin\alpha \sin\beta \\
    Q_B &= q \sin\alpha \cos\beta \\
    Q_C &= q \cos\alpha


where $\rho_\text{core}$, $\rho_\text{A}$, $\rho_\text{B}$ and $\rho_\text{C}$
are the scattering lengths of the parallelepiped core, and the rectangular
slabs of thickness $t_A$, $t_B$ and $t_C$, respectively. $\rho_\text{solvent}$
is the scattering length of the solvent.

.. note::

   the code actually implements two substitutions: $d(cos\alpha)$ is
   substituted for -$sin\alpha \ d\alpha$ (note that in the
   :ref:`parallelepiped` code this is explicitly implemented with
   $\sigma = cos\alpha$), and $\beta$ is set to $\beta = u \pi/2$ so that
   $du = \pi/2 \ d\beta$.  Thus both integrals go from 0 to 1 rather than 0
   to $\pi/2$.

FITTING NOTES
~~~~~~~~~~~~~

#. There are many parameters in this model. Hold as many fixed as possible with
   known values, or you will certainly end up at a solution that is unphysical.

#. The 2nd virial coefficient of the core_shell_parallelepiped is calculated
   based on the the averaged effective radius $(=\sqrt{(A+2t_A)(B+2t_B)/\pi})$
   and length $(C+2t_C)$ values, after appropriately sorting the three
   dimensions to give an oblate or prolate particle, to give an effective radius
   for $S(q)$ when $P(q) * S(q)$ is applied.

#. For 2d data the orientation of the particle is required, described using
   angles $\theta$, $\phi$ and $\Psi$ as in the diagrams below, where $\theta$
   and $\phi$ define the orientation of the director in the laboratory reference
   frame of the beam direction ($z$) and detector plane ($x-y$ plane), while
   the angle $\Psi$ is effectively the rotational angle around the particle
   $C$ axis. For $\theta = 0$ and $\phi = 0$, $\Psi = 0$ corresponds to the
   $B$ axis oriented parallel to the y-axis of the detector with $A$ along
   the x-axis. For other $\theta$, $\phi$ values, the order of rotations
   matters. In particular, the parallelepiped must first be rotated $\theta$
   degrees in the $x-z$ plane before rotating $\phi$ degrees around the $z$
   axis (in the $x-y$ plane). Applying orientational distribution to the
   particle orientation (i.e  `jitter` to one or more of these angles) can get
   more confusing as `jitter` is defined **NOT** with respect to the laboratory
   frame but the particle reference frame. It is thus highly recommended to
   read :ref:`orientation` for further details of the calculation and angular
   dispersions.

.. note:: For 2d, constraints must be applied during fitting to ensure that the
   order of sides chosen is not altered, and hence that the correct definition
   of angles is preserved. For the default choice shown here, that means
   ensuring that the inequality $A < B < C$ is not violated,  The calculation
   will not report an error, but the results may be not correct.

.. figure:: img/parallelepiped_angle_definition.png

    Definition of the angles for oriented core-shell parallelepipeds.
    Note that rotation $\theta$, initially in the $x-z$ plane, is carried
    out first, then rotation $\phi$ about the $z$ axis, finally rotation
    $\Psi$ is now around the $C$ axis of the particle. The neutron or X-ray
    beam is along the $z$ axis and the detecotr defines the $x-y$ plane.

.. figure:: img/parallelepiped_angle_projection.png

    Examples of the angles for oriented core-shell parallelepipeds against the
    detector plane.


Validation
----------

Cross-checked against hollow rectangular prism and rectangular prism for equal
thickness overlapping sides, and by Monte Carlo sampling of points within the
shape for non-uniform, non-overlapping sides.


References
----------

#. P Mittelbach and G Porod, *Acta Physica Austriaca*, 14 (1961) 185-211
   Equations (1), (13-14). (in German)
#. D Singh (2009). *Small angle scattering studies of self assembly in
   lipid mixtures*, Johns Hopkins University Thesis (2009) 223-225. `Available
   from Proquest <http://search.proquest.com/docview/304915826>`_
#.  L. Onsager, *Ann. New York Acad. Sci.*, 51 (1949) 627-659

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Converted to sasmodels by:** Miguel Gonzalez **Date:** February 26, 2016
* **Last Modified by:** Paul Kienzle **Date:** October 17, 2017
* **Last Reviewed by:** Paul Butler **Date:** May 24, 2018 - documentation
  updated
"""

import numpy as np
from numpy import inf

name = "core_shell_parallelepiped"
title = "Rectangular solid with a core-shell structure."
description = """
     P(q)=
"""
category = "shape:parallelepiped"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld_core", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Parallelepiped core scattering length density"],
              ["sld_a", "1e-6/Ang^2", 2, [-inf, inf], "sld",
               "Parallelepiped A rim scattering length density"],
              ["sld_b", "1e-6/Ang^2", 4, [-inf, inf], "sld",
               "Parallelepiped B rim scattering length density"],
              ["sld_c", "1e-6/Ang^2", 2, [-inf, inf], "sld",
               "Parallelepiped C rim scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 6, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["length_a", "Ang", 35, [0, inf], "volume",
               "Shorter side of the parallelepiped"],
              ["length_b", "Ang", 75, [0, inf], "volume",
               "Second side of the parallelepiped"],
              ["length_c", "Ang", 400, [0, inf], "volume",
               "Larger side of the parallelepiped"],
              ["thick_rim_a", "Ang", 10, [0, inf], "volume",
               "Thickness of A rim"],
              ["thick_rim_b", "Ang", 10, [0, inf], "volume",
               "Thickness of B rim"],
              ["thick_rim_c", "Ang", 10, [0, inf], "volume",
               "Thickness of C rim"],
              ["theta", "degrees", 0, [-360, 360], "orientation",
               "c axis to beam angle"],
              ["phi", "degrees", 0, [-360, 360], "orientation",
               "rotation about beam"],
              ["psi", "degrees", 0, [-360, 360], "orientation",
               "rotation about c axis"],
             ]

source = ["lib/gauss76.c", "core_shell_parallelepiped.c"]
has_shape_visualization = True

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for core-shell parallelepiped visualization."""
    import numpy as np
    length_a = params.get('length_a', 35)
    length_b = params.get('length_b', 75)
    length_c = params.get('length_c', 400)
    thick_a = params.get('thick_rim_a', 10)
    thick_b = params.get('thick_rim_b', 10)
    thick_c = params.get('thick_rim_c', 10)

    faces = {}

    # Outer dimensions (core + shells on each face)
    outer_a = length_a + 2 * thick_a
    outer_b = length_b + 2 * thick_b
    outer_c = length_c + 2 * thick_c

    # Create outer box faces
    y = np.linspace(-outer_b/2, outer_b/2, resolution//4)
    z = np.linspace(-outer_c/2, outer_c/2, resolution//2)
    y_mesh, z_mesh = np.meshgrid(y, z)
    faces['outer_front'] = (np.full_like(y_mesh, outer_a/2), y_mesh, z_mesh)
    faces['outer_back'] = (np.full_like(y_mesh, -outer_a/2), y_mesh, z_mesh)

    x = np.linspace(-outer_a/2, outer_a/2, resolution//4)
    z = np.linspace(-outer_c/2, outer_c/2, resolution//2)
    x_mesh, z_mesh = np.meshgrid(x, z)
    faces['outer_right'] = (x_mesh, np.full_like(x_mesh, outer_b/2), z_mesh)
    faces['outer_left'] = (x_mesh, np.full_like(x_mesh, -outer_b/2), z_mesh)

    x = np.linspace(-outer_a/2, outer_a/2, resolution//4)
    y = np.linspace(-outer_b/2, outer_b/2, resolution//4)
    x_mesh, y_mesh = np.meshgrid(x, y)
    faces['outer_top'] = (x_mesh, y_mesh, np.full_like(x_mesh, outer_c/2))
    faces['outer_bottom'] = (x_mesh, y_mesh, np.full_like(x_mesh, -outer_c/2))

    return faces

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of the core-shell parallelepiped."""
    import numpy as np
    length_a = params.get('length_a', 35)
    length_b = params.get('length_b', 75)
    length_c = params.get('length_c', 400)
    thick_a = params.get('thick_rim_a', 10)
    thick_b = params.get('thick_rim_b', 10)
    thick_c = params.get('thick_rim_c', 10)

    outer_a = length_a + 2 * thick_a
    outer_b = length_b + 2 * thick_b
    outer_c = length_c + 2 * thick_c

    # XY plane (top view) - shows A and B dimensions
    # Outer rectangle
    outer_x = [-outer_a/2, -outer_a/2, outer_a/2, outer_a/2, -outer_a/2]
    outer_y = [-outer_b/2, outer_b/2, outer_b/2, -outer_b/2, -outer_b/2]
    # Core rectangle
    core_x = [-length_a/2, -length_a/2, length_a/2, length_a/2, -length_a/2]
    core_y = [-length_b/2, length_b/2, length_b/2, -length_b/2, -length_b/2]

    ax_xy.fill(outer_x, outer_y, 'lightcoral', alpha=0.3, label='Shell')
    ax_xy.fill(core_x, core_y, 'lightblue', alpha=0.5, label='Core')
    ax_xy.plot(outer_x, outer_y, 'r-', linewidth=2)
    ax_xy.plot(core_x, core_y, 'b-', linewidth=2)
    max_xy = max(outer_a, outer_b) / 2 * 1.3
    ax_xy.set_xlim(-max_xy, max_xy)
    ax_xy.set_ylim(-max_xy, max_xy)
    ax_xy.set_xlabel('X (Å) - A')
    ax_xy.set_ylabel('Y (Å) - B')
    ax_xy.set_title(f'XY Cross-section (t_a={thick_a:.0f}, t_b={thick_b:.0f}Å)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)
    ax_xy.legend(fontsize=8)

    # XZ plane (side view) - shows A and C dimensions
    outer_xz_x = [-outer_a/2, -outer_a/2, outer_a/2, outer_a/2, -outer_a/2]
    outer_xz_z = [-outer_c/2, outer_c/2, outer_c/2, -outer_c/2, -outer_c/2]
    core_xz_x = [-length_a/2, -length_a/2, length_a/2, length_a/2, -length_a/2]
    core_xz_z = [-length_c/2, length_c/2, length_c/2, -length_c/2, -length_c/2]

    ax_xz.fill(outer_xz_x, outer_xz_z, 'lightcoral', alpha=0.3)
    ax_xz.fill(core_xz_x, core_xz_z, 'lightblue', alpha=0.5)
    ax_xz.plot(outer_xz_x, outer_xz_z, 'r-', linewidth=2)
    ax_xz.plot(core_xz_x, core_xz_z, 'b-', linewidth=2)
    max_xz = max(outer_a, outer_c) / 2 * 1.1
    ax_xz.set_xlim(-max_xz, max_xz)
    ax_xz.set_ylim(-max_xz, max_xz)
    ax_xz.set_xlabel('X (Å) - A')
    ax_xz.set_ylabel('Z (Å) - C')
    ax_xz.set_title(f'XZ Cross-section (t_c={thick_c:.0f}Å)')
    ax_xz.grid(True, alpha=0.3)

    # YZ plane (front view) - shows B and C dimensions
    outer_yz_y = [-outer_b/2, -outer_b/2, outer_b/2, outer_b/2, -outer_b/2]
    outer_yz_z = [-outer_c/2, outer_c/2, outer_c/2, -outer_c/2, -outer_c/2]
    core_yz_y = [-length_b/2, -length_b/2, length_b/2, length_b/2, -length_b/2]
    core_yz_z = [-length_c/2, length_c/2, length_c/2, -length_c/2, -length_c/2]

    ax_yz.fill(outer_yz_y, outer_yz_z, 'lightgreen', alpha=0.3)
    ax_yz.fill(core_yz_y, core_yz_z, 'moccasin', alpha=0.5)
    ax_yz.plot(outer_yz_y, outer_yz_z, 'g-', linewidth=2)
    ax_yz.plot(core_yz_y, core_yz_z, 'orange', linewidth=2)
    max_yz = max(outer_b, outer_c) / 2 * 1.1
    ax_yz.set_xlim(-max_yz, max_yz)
    ax_yz.set_ylim(-max_yz, max_yz)
    ax_yz.set_xlabel('Y (Å) - B')
    ax_yz.set_ylabel('Z (Å) - C')
    ax_yz.set_title('YZ Cross-section')
    ax_yz.grid(True, alpha=0.3)
have_Fq = True
radius_effective_modes = [
    "equivalent cylinder excluded volume",
    "equivalent volume sphere",
    "half outer length_a", "half outer length_b", "half outer length_c",
    "equivalent circular cross-section",
    "half outer ab diagonal", "half outer diagonal",
    ]

def random():
    """Return a random parameter set for the model."""
    outer = 10**np.random.uniform(1, 4.7, size=3)
    thick = np.random.beta(0.5, 0.5, size=3)*(outer-2) + 1
    length = outer - thick
    pars = dict(
        length_a=length[0],
        length_b=length[1],
        length_c=length[2],
        thick_rim_a=thick[0],
        thick_rim_b=thick[1],
        thick_rim_c=thick[2],
    )
    return pars

# rkh 7/4/17 add random unit test for 2d, note make all params different,
# 2d values not tested against other codes or models
if 0:  # pak: model rewrite; need to update tests
    qx, qy = 0.2 * np.cos(np.pi/6.), 0.2 * np.sin(np.pi/6.)
    tests = [[{}, 0.2, 0.533149288477],
             [{}, [0.2], [0.533149288477]],
             [{'theta':10.0, 'phi':20.0}, (qx, qy), 0.0853299803222],
             [{'theta':10.0, 'phi':20.0}, [(qx, qy)], [0.0853299803222]],
            ]
    del qx, qy  # not necessary to delete, but cleaner
