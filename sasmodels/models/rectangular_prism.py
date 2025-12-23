# rectangular_prism model
# Note: model title and parameter table are inserted automatically
r"""

This model provides the form factor, $P(q)$, for a rectangular prism.

Note that this model is almost totally equivalent to the existing
:ref:`parallelepiped` model.
The only difference is that the way the relevant
parameters are defined here ($a$, $b/a$, $c/a$ instead of $a$, $b$, $c$)
which allows use of polydispersity with this model while keeping the shape of
the prism (e.g. setting $b/a = 1$ and $c/a = 1$ and applying polydispersity
to *a* will generate a distribution of cubes of different sizes).


Definition
----------

The 1D scattering intensity for this model was calculated by (Mittelbach and
Porod, 1961 [#Mittelbach1961]_), but the implementation here is closer to
the equations given by (Nayuk and Huber, 2012 [#Nayuk2012]_).
Note also that the angle definitions used in the code and the present
documentation correspond to those used in Nayuk and Huber (see Fig. 1 of
that reference), with $\theta$ corresponding to $\alpha$ in that paper,
and not to the usual convention used for example in the
:ref:`parallelepiped` model.

In this model the scattering from a massive parallelepiped with an
orientation with respect to the scattering vector given by $\theta$
and $\phi$

.. math::

    A_P\,(q) =
      \frac{\sin \left( \tfrac{1}{2}qC \cos\theta \right) }{\tfrac{1}{2} qC \cos\theta}
      \,\times\,
      \frac{\sin \left( \tfrac{1}{2}qA \cos\theta \right) }{\tfrac{1}{2} qA \cos\theta}
      \,\times\ ,
      \frac{\sin \left( \tfrac{1}{2}qB \cos\theta \right) }{\tfrac{1}{2} qB \cos\theta}

where $A$, $B$ and $C$ are the sides of the parallelepiped and must fulfill
$A \le B \le C$, $\theta$ is the angle between the $z$ axis and the
$C$ axis of the parallelepiped, and $\phi$ is the angle between the
scattering vector (lying in the $xy$ plane) and the $y$ axis.

The normalized form factor in 1D is obtained averaging over all possible
orientations

.. math::
    P(q) =  \frac{2}{\pi} \int_0^{\frac{\pi}{2}} \,
    \int_0^{\frac{\pi}{2}} A_P^2(q) \, \sin\theta \, d\theta \, d\phi

And the 1D scattering intensity is calculated as

.. math::
    I(q) = \text{scale} \times V \times (\rho_\text{p} -
    \rho_\text{solvent})^2 \times P(q)

where $V$ is the volume of the rectangular prism, $\rho_\text{p}$
is the scattering length of the parallelepiped, $\rho_\text{solvent}$
is the scattering length of the solvent, and (if the data are in absolute
units) *scale* represents the volume fraction (which is unitless).

For 2d data the orientation of the particle is required, described using
angles $\theta$, $\phi$ and $\Psi$ as in the diagrams below, for further details
of the calculation and angular dispersions see :ref:`orientation`.
The angle $\Psi$ is the rotational angle around the long *C* axis. For example,
$\Psi = 0$ when the *B* axis is parallel to the *x*-axis of the detector.

For 2d, constraints must be applied during fitting to ensure that the
inequality $A < B < C$ is not violated, and hence the correct definition of
angles is preserved. The calculation will not report an error, but the
results may be not correct.

.. figure:: img/parallelepiped_angle_definition.png

    Definition of the angles for oriented core-shell parallelepipeds. Note
    that rotation $\theta$, initially in the $xz$ plane, is carried out
    first, then rotation $\phi$ about the $z$ axis, finally rotation $\Psi$
    is now around the axis of the cylinder. The neutron or X-ray beam is
    along the $z$ axis.

.. figure:: img/parallelepiped_angle_projection.png

    Examples of the angles for oriented rectangular prisms against the
    detector plane.


Validation
----------

Validation of the code was conducted by comparing the output of the 1D model
to the output of the existing :ref:`parallelepiped` model.


References
----------

See also Onsager [#Onsager1949]_.

.. [#Mittelbach1961] P Mittelbach and G Porod, *Acta Physica Austriaca*,
   14 (1961) 185-211
.. [#Nayuk2012] R Nayuk and K Huber, *Z. Phys. Chem.*, 226 (2012) 837-854
.. [#Onsager1949] L Onsager, *Ann. New York Acad. Sci.*, 51 (1949) 627-659

Authorship and Verification
----------------------------

* **Author:**
* **Last Modified by:**
* **Last Reviewed by:**
"""

import numpy as np
from numpy import inf

name = "rectangular_prism"
title = "Rectangular parallelepiped with uniform scattering length density."
description = """
    I(q)= scale*V*(sld - sld_solvent)^2*P(q,theta,phi)+background
        P(q,theta,phi) = (2/pi) * double integral from 0 to pi/2 of ...
           AP^2(q)*sin(theta)*dtheta*dphi
        AP = S(q*C*cos(theta)/2) * S(q*A*sin(theta)*sin(phi)/2) * S(q*B*sin(theta)*cos(phi)/2)
        S(x) = sin(x)/x
"""
category = "shape:parallelepiped"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 6.3, [-inf, inf], "sld",
               "Parallelepiped scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["length_a", "Ang", 35, [0, inf], "volume",
               "Shorter side of the parallelepiped"],
              ["b2a_ratio", "", 1, [0, inf], "volume",
               "Ratio sides b/a"],
              ["c2a_ratio", "", 1, [0, inf], "volume",
               "Ratio sides c/a"],
              ["theta", "degrees", 0, [-360, 360], "orientation",
               "c axis to beam angle"],
              ["phi", "degrees", 0, [-360, 360], "orientation",
               "rotation about beam"],
              ["psi", "degrees", 0, [-360, 360], "orientation",
               "rotation about c axis"],
             ]

source = ["lib/gauss76.c", "rectangular_prism.c"]
has_shape_visualization = True

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for rectangular prism visualization."""
    import numpy as np
    length_a = params.get('length_a', 35)
    b2a_ratio = params.get('b2a_ratio', 1)
    c2a_ratio = params.get('c2a_ratio', 1)

    length_b = length_a * b2a_ratio
    length_c = length_a * c2a_ratio

    faces = {}
    # Front and back faces (YZ planes at x = ±a/2)
    y = np.linspace(-length_b/2, length_b/2, resolution//4)
    z = np.linspace(-length_c/2, length_c/2, resolution//2)
    y_mesh, z_mesh = np.meshgrid(y, z)
    faces['front'] = (np.full_like(y_mesh, length_a/2), y_mesh, z_mesh)
    faces['back'] = (np.full_like(y_mesh, -length_a/2), y_mesh, z_mesh)

    # Left and right faces (XZ planes at y = ±b/2)
    x = np.linspace(-length_a/2, length_a/2, resolution//4)
    z = np.linspace(-length_c/2, length_c/2, resolution//2)
    x_mesh, z_mesh = np.meshgrid(x, z)
    faces['right'] = (x_mesh, np.full_like(x_mesh, length_b/2), z_mesh)
    faces['left'] = (x_mesh, np.full_like(x_mesh, -length_b/2), z_mesh)

    # Top and bottom faces (XY planes at z = ±c/2)
    x = np.linspace(-length_a/2, length_a/2, resolution//4)
    y = np.linspace(-length_b/2, length_b/2, resolution//4)
    x_mesh, y_mesh = np.meshgrid(x, y)
    faces['top'] = (x_mesh, y_mesh, np.full_like(x_mesh, length_c/2))
    faces['bottom'] = (x_mesh, y_mesh, np.full_like(x_mesh, -length_c/2))

    return faces

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of the rectangular prism."""
    import numpy as np
    length_a = params.get('length_a', 35)
    b2a_ratio = params.get('b2a_ratio', 1)
    c2a_ratio = params.get('c2a_ratio', 1)

    length_b = length_a * b2a_ratio
    length_c = length_a * c2a_ratio

    # XY plane (A x B)
    rect_x = [-length_a/2, -length_a/2, length_a/2, length_a/2, -length_a/2]
    rect_y = [-length_b/2, length_b/2, length_b/2, -length_b/2, -length_b/2]
    ax_xy.plot(rect_x, rect_y, 'b-', linewidth=2)
    ax_xy.fill(rect_x, rect_y, 'lightblue', alpha=0.3)
    max_xy = max(length_a, length_b) / 2 * 1.3
    ax_xy.set_xlim(-max_xy, max_xy)
    ax_xy.set_ylim(-max_xy, max_xy)
    ax_xy.set_xlabel('X (Å) - A')
    ax_xy.set_ylabel('Y (Å) - B')
    ax_xy.set_title('XY Cross-section (Top View)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)

    # Add dimensions
    ax_xy.text(0, -max_xy*0.9, f'A={length_a:.0f}Å, B={length_b:.0f}Å', ha='center', fontsize=9)

    # XZ plane (A x C)
    rect_x_xz = [-length_a/2, -length_a/2, length_a/2, length_a/2, -length_a/2]
    rect_z = [-length_c/2, length_c/2, length_c/2, -length_c/2, -length_c/2]
    ax_xz.plot(rect_x_xz, rect_z, 'r-', linewidth=2)
    ax_xz.fill(rect_x_xz, rect_z, 'lightcoral', alpha=0.3)
    max_xz = max(length_a, length_c) / 2 * 1.3
    ax_xz.set_xlim(-max_xz, max_xz)
    ax_xz.set_ylim(-max_xz, max_xz)
    ax_xz.set_xlabel('X (Å) - A')
    ax_xz.set_ylabel('Z (Å) - C')
    ax_xz.set_title('XZ Cross-section (Side View)')
    ax_xz.grid(True, alpha=0.3)

    # YZ plane (B x C)
    rect_y_yz = [-length_b/2, -length_b/2, length_b/2, length_b/2, -length_b/2]
    rect_z_yz = [-length_c/2, length_c/2, length_c/2, -length_c/2, -length_c/2]
    ax_yz.plot(rect_y_yz, rect_z_yz, 'g-', linewidth=2)
    ax_yz.fill(rect_y_yz, rect_z_yz, 'lightgreen', alpha=0.3)
    max_yz = max(length_b, length_c) / 2 * 1.3
    ax_yz.set_xlim(-max_yz, max_yz)
    ax_yz.set_ylim(-max_yz, max_yz)
    ax_yz.set_xlabel('Y (Å) - B')
    ax_yz.set_ylabel('Z (Å) - C')
    ax_yz.set_title('YZ Cross-section (Front View)')
    ax_yz.grid(True, alpha=0.3)
have_Fq = True
radius_effective_modes = [
    "equivalent cylinder excluded volume", "equivalent volume sphere",
    "half length_a", "half length_b", "half length_c",
    "equivalent circular cross-section", "half ab diagonal", "half diagonal",
    ]

def random():
    """Return a random parameter set for the model."""
    a, b, c = 10**np.random.uniform(1, 4.7, size=3)
    pars = dict(
        length_a=a,
        b2a_ratio=b/a,
        c2a_ratio=c/a,
    )
    return pars

tests = [[{}, 0.2, 0.375248406825],
         [{}, [0.2], [0.375248406825]],
        ]
