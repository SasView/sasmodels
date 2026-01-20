# rectangular_prism model
# Note: model title and parameter table are inserted automatically
r"""
Definition
----------

This model provides the form factor, $P(q)$, for a hollow rectangular
parallelepiped with a wall of thickness $\Delta$. The 1D scattering intensity
for this model is calculated by forming the difference of the amplitudes of two
massive parallelepipeds differing in their outermost dimensions in each
direction by the same length increment $2\Delta$ (\ [#Nayuk2012]_ Nayuk, 2012).

As in the case of the massive parallelepiped model (:ref:`rectangular-prism`),
the scattering amplitude is computed for a particular orientation of the
parallelepiped with respect to the scattering vector and then averaged over all
possible orientations, giving

.. math::
  P(q) =  \frac{1}{V^2} \frac{2}{\pi} \times \, \int_0^{\frac{\pi}{2}} \,
  \int_0^{\frac{\pi}{2}} A_{P\Delta}^2(q) \, \sin\theta \, d\theta \, d\phi

where $\theta$ is the angle between the $z$ axis and the $C$ axis
of the parallelepiped, $\phi$ is the angle between the scattering vector
(lying in the $xy$ plane) and the $y$ axis, and

.. math::
  :nowrap:

  \begin{align*}
  A_{P\Delta}(q) & =  A B C
    \left[\frac{\sin \bigl( q \frac{C}{2} \cos\theta \bigr)}
    {\left( q \frac{C}{2} \cos\theta \right)} \right]
    \left[\frac{\sin \bigl( q \frac{A}{2} \sin\theta \sin\phi \bigr)}
    {\left( q \frac{A}{2} \sin\theta \sin\phi \right)}\right]
    \left[\frac{\sin \bigl( q \frac{B}{2} \sin\theta \cos\phi \bigr)}
    {\left( q \frac{B}{2} \sin\theta \cos\phi \right)}\right] \\
    & {} - A' B' C'
    \left[ \frac{\sin \bigl[ q \frac{C'}{2} \cos\theta \bigr]}
    {q \frac{C'}{2} \cos\theta} \right]
    \left[ \frac{\sin \bigl[ q \frac{A'}{2} \sin\theta \sin\phi \bigr]}
    {q \frac{A'}{2} \sin\theta \sin\phi} \right]
    \left[ \frac{\sin \bigl[ q \frac{B'}{2} \sin\theta \cos\phi \bigr]}
    {q \frac{B'}{2} \sin\theta \cos\phi} \right]
  \end{align*}

$A'$, $B'$, $C'$ are the inner dimensions and $A = A' + 2\Delta$,
$B = B' + 2\Delta$, $C = C' + 2\Delta$ are the outer dimensions of the
parallelepiped shell, giving a shell volume of $V = ABC - A'B'C'$.

The 1D scattering intensity is then calculated as

.. math::
  I(q) = \text{scale} \times V \times (\rho_\text{p} -
  \rho_\text{solvent})^2 \times P(q) + \text{background}

where $\rho_\text{p}$ is the scattering length density of the parallelepiped,
$\rho_\text{solvent}$ is the scattering length density of the solvent,
and (if the data are in absolute units) *scale* represents the volume fraction
(which is unitless) of the rectangular shell of material (i.e. not including
the volume of the solvent filled core).

For 2d data the orientation of the particle is required, described using
angles $\theta$, $\phi$ and $\Psi$ as in the diagrams below, for further details
of the calculation and angular dispersions see :ref:`orientation`.
The angle $\Psi$ is the rotational angle around the long *C* axis. For example,
$\Psi = 0$ when the *B* axis is parallel to the *x*-axis of the detector.

For 2d, constraints must be applied during fitting to ensure that the inequality
$A < B < C$ is not violated, and hence the correct definition of angles is
preserved. The calculation will not report an error if the inequality is *not*
preserved, but the results may be not correct.

.. figure:: img/parallelepiped_angle_definition.png

    Definition of the angles for oriented hollow rectangular prism. Note that
    rotation $\theta$, initially in the $xz$ plane, is carried out first,
    then rotation $\phi$ about the $z$ axis, finally rotation $\Psi$ is now
    around the axis of the prism. The neutron or X-ray beam is along the $z$
    axis.

.. figure:: img/parallelepiped_angle_projection.png

    Examples of the angles for oriented hollow rectangular prisms against the
    detector plane.


Validation
----------

Validation of the code was conducted by qualitatively comparing the output
of the 1D model to the curves shown in (Nayuk, 2012).


References
----------

See also Onsager [#Onsager1949]_.

.. [#Nayuk2012] R Nayuk and K Huber, *Z. Phys. Chem.*, 226 (2012) 837-854
.. [#Onsager1949] L. Onsager, *Ann. New York Acad. Sci.*, 51 (1949) 627-659

Authorship and Verification
----------------------------

* **Author:** Miguel Gonzales **Date:** February 26, 2016
* **Last Modified by:** Paul Kienzle **Date:** December 14, 2017
* **Last Reviewed by:** Paul Butler **Date:** September 06, 2018
"""

import numpy as np
from numpy import inf

name = "hollow_rectangular_prism"
title = "Hollow rectangular parallelepiped with uniform scattering length density."
description = """
    I(q)= scale*V*(sld - sld_solvent)^2*P(q,theta,phi)+background
        P(q,theta,phi) = (2/pi/V^2) * double integral from 0 to pi/2 of ...
           (AP1-AP2)^2(q)*sin(theta)*dtheta*dphi
        AP1 = S(q*C*cos(theta)/2) * S(q*A*sin(theta)*sin(phi)/2) * S(q*B*sin(theta)*cos(phi)/2)
        AP2 = S(q*C'*cos(theta)) * S(q*A'*sin(theta)*sin(phi)) * S(q*B'*sin(theta)*cos(phi))
        C' = (C/2-thickness)
        B' = (B/2-thickness)
        A' = (A/2-thickness)
        S(x) = sin(x)/x
"""
category = "shape:parallelepiped"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 6.3, [-inf, inf], "sld",
               "Parallelepiped scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["length_a", "Ang", 35, [0, inf], "volume",
               "Shortest, external, size of the parallelepiped"],
              ["b2a_ratio", "Ang", 1, [0, inf], "volume",
               "Ratio sides b/a"],
              ["c2a_ratio", "Ang", 1, [0, inf], "volume",
               "Ratio sides c/a"],
              ["thickness", "Ang", 1, [0, inf], "volume",
               "Thickness of parallelepiped"],
              ["theta", "degrees", 0, [-360, 360], "orientation",
               "c axis to beam angle"],
              ["phi", "degrees", 0, [-360, 360], "orientation",
               "rotation about beam"],
              ["psi", "degrees", 0, [-360, 360], "orientation",
               "rotation about c axis"],
             ]

source = ["lib/gauss76.c", "hollow_rectangular_prism.c"]
has_shape_visualization = True

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for hollow rectangular prism visualization."""
    import numpy as np
    length_a = params.get('length_a', 35)
    b2a_ratio = params.get('b2a_ratio', 1)
    c2a_ratio = params.get('c2a_ratio', 1)
    thickness = params.get('thickness', 1)

    # Outer dimensions
    outer_a = length_a
    outer_b = length_a * b2a_ratio
    outer_c = length_a * c2a_ratio

    # Inner dimensions
    inner_a = outer_a - 2 * thickness
    inner_b = outer_b - 2 * thickness
    inner_c = outer_c - 2 * thickness

    # Ensure inner dimensions are positive
    inner_a = max(inner_a, 0.1)
    inner_b = max(inner_b, 0.1)
    inner_c = max(inner_c, 0.1)

    faces = {}

    # Outer faces
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
    """Plot 2D cross-sections of the hollow rectangular prism."""
    length_a = params.get('length_a', 35)
    b2a_ratio = params.get('b2a_ratio', 1)
    c2a_ratio = params.get('c2a_ratio', 1)
    thickness = params.get('thickness', 1)

    outer_a, outer_b, outer_c = length_a, length_a * b2a_ratio, length_a * c2a_ratio
    inner_a = max(outer_a - 2*thickness, 0.1)
    inner_b = max(outer_b - 2*thickness, 0.1)
    inner_c = max(outer_c - 2*thickness, 0.1)

    # XY plane
    outer_x = [-outer_a/2, -outer_a/2, outer_a/2, outer_a/2, -outer_a/2]
    outer_y = [-outer_b/2, outer_b/2, outer_b/2, -outer_b/2, -outer_b/2]
    inner_x = [-inner_a/2, -inner_a/2, inner_a/2, inner_a/2, -inner_a/2]
    inner_y = [-inner_b/2, inner_b/2, inner_b/2, -inner_b/2, -inner_b/2]

    ax_xy.fill(outer_x, outer_y, 'lightblue', alpha=0.3)
    ax_xy.fill(inner_x, inner_y, 'white', alpha=1.0)
    ax_xy.plot(outer_x, outer_y, 'b-', linewidth=2, label='Outer')
    ax_xy.plot(inner_x, inner_y, 'r--', linewidth=2, label='Inner (hollow)')
    max_xy = max(outer_a, outer_b) / 2 * 1.3
    ax_xy.set_xlim(-max_xy, max_xy)
    ax_xy.set_ylim(-max_xy, max_xy)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title('XY Cross-section (Top View)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)
    ax_xy.legend()

    # XZ plane
    outer_xz_x = [-outer_a/2, -outer_a/2, outer_a/2, outer_a/2, -outer_a/2]
    outer_xz_z = [-outer_c/2, outer_c/2, outer_c/2, -outer_c/2, -outer_c/2]
    inner_xz_x = [-inner_a/2, -inner_a/2, inner_a/2, inner_a/2, -inner_a/2]
    inner_xz_z = [-inner_c/2, inner_c/2, inner_c/2, -inner_c/2, -inner_c/2]

    ax_xz.fill(outer_xz_x, outer_xz_z, 'lightcoral', alpha=0.3)
    ax_xz.fill(inner_xz_x, inner_xz_z, 'white', alpha=1.0)
    ax_xz.plot(outer_xz_x, outer_xz_z, 'r-', linewidth=2)
    ax_xz.plot(inner_xz_x, inner_xz_z, 'b--', linewidth=2)
    max_xz = max(outer_a, outer_c) / 2 * 1.3
    ax_xz.set_xlim(-max_xz, max_xz)
    ax_xz.set_ylim(-max_xz, max_xz)
    ax_xz.set_xlabel('X (Å)')
    ax_xz.set_ylabel('Z (Å)')
    ax_xz.set_title('XZ Cross-section (Side View)')
    ax_xz.grid(True, alpha=0.3)

    # YZ plane
    outer_yz_y = [-outer_b/2, -outer_b/2, outer_b/2, outer_b/2, -outer_b/2]
    outer_yz_z = [-outer_c/2, outer_c/2, outer_c/2, -outer_c/2, -outer_c/2]
    inner_yz_y = [-inner_b/2, -inner_b/2, inner_b/2, inner_b/2, -inner_b/2]
    inner_yz_z = [-inner_c/2, inner_c/2, inner_c/2, -inner_c/2, -inner_c/2]

    ax_yz.fill(outer_yz_y, outer_yz_z, 'lightgreen', alpha=0.3)
    ax_yz.fill(inner_yz_y, inner_yz_z, 'white', alpha=1.0)
    ax_yz.plot(outer_yz_y, outer_yz_z, 'g-', linewidth=2)
    ax_yz.plot(inner_yz_y, inner_yz_z, 'orange', linewidth=2)
    max_yz = max(outer_b, outer_c) / 2 * 1.3
    ax_yz.set_xlim(-max_yz, max_yz)
    ax_yz.set_ylim(-max_yz, max_yz)
    ax_yz.set_xlabel('Y (Å)')
    ax_yz.set_ylabel('Z (Å)')
    ax_yz.set_title('YZ Cross-section (Front View)')
    ax_yz.grid(True, alpha=0.3)
have_Fq = True
radius_effective_modes = [
    "equivalent cylinder excluded volume", "equivalent outer volume sphere",
    "half length_a", "half length_b", "half length_c",
    "equivalent outer circular cross-section",
    "half ab diagonal", "half diagonal",
    ]

def random():
    """Return a random parameter set for the model."""
    a, b, c = 10**np.random.uniform(1, 4.7, size=3)
    # Thickness is limited to 1/2 the smallest dimension
    # Use a distribution with a preference for thin shell or thin core
    # Avoid core,shell radii < 1
    min_dim = 0.5*min(a, b, c)
    thickness = np.random.beta(0.5, 0.5)*(min_dim-2) + 1
    #print(a, b, c, thickness, thickness/min_dim)
    pars = dict(
        length_a=a,
        b2a_ratio=b/a,
        c2a_ratio=c/a,
        thickness=thickness,
    )
    return pars

tests = [[{}, 0.2, 0.76687283098],
         [{}, [0.2], [0.76687283098]],
        ]
