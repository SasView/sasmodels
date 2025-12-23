r"""
Definition
----------

This model provides the form factor for an elliptical cylinder with a
core-shell scattering length density profile [#Onsager1949]_.
Thus this is a variation
of the core-shell bicelle model, but with an elliptical cylinder for the core.
Outer shells on the rims and flat ends may be of different thicknesses and
scattering length densities. The form factor is normalized by the total
particle volume.

.. figure:: img/core_shell_bicelle_geometry.png

    Schematic cross-section of bicelle. Note however that the model here
    calculates for rectangular, not curved, rims as shown below.

.. figure:: img/core_shell_bicelle_parameters.png

   Cross section of model used here. Users will have
   to decide how to distribute "heads" and "tails" between the rim, face
   and core regions in order to estimate appropriate starting parameters.

Given the scattering length densities (sld) $\rho_c$, the core sld, $\rho_f$,
the face sld, $\rho_r$, the rim sld and $\rho_s$ the solvent sld, the
scattering length density variation along the bicelle axis is:

.. math::

    \rho(r) =
      \begin{cases}
      &\rho_c \text{ for } 0 \lt r \lt R; -L \lt z\lt L \\[1.5ex]
      &\rho_f \text{ for } 0 \lt r \lt R; -(L+2t) \lt z\lt -L;
      L \lt z\lt (L+2t) \\[1.5ex]
      &\rho_r\text{ for } 0 \lt r \lt R; -(L+2t) \lt z\lt -L; L \lt z\lt (L+2t)
      \end{cases}

The form factor for the bicelle is calculated in cylindrical coordinates, where
$\alpha$ is the angle between the $Q$ vector and the cylinder axis, and $\psi$
is the angle for the ellipsoidal cross section core, to give:

.. math::

    I(Q,\alpha,\psi) = \frac{\text{scale}}{V_t} \cdot
        F(Q,\alpha, \psi)^2 \cdot sin(\alpha) + \text{background}

where a numerical integration of $F(Q,\alpha, \psi)^2 \cdot sin(\alpha)$
is carried out over \alpha and \psi for:

.. math::
    :nowrap:

    \begin{align*}
    F(Q,\alpha,\psi) = &\bigg[
    (\rho_c - \rho_f) V_c
     \frac{2J_1(QR'sin \alpha)}{QR'sin\alpha}
     \frac{sin(QLcos\alpha/2)}{Q(L/2)cos\alpha} \\
    &+(\rho_f - \rho_r) V_{c+f}
     \frac{2J_1(QR'sin\alpha)}{QR'sin\alpha}
     \frac{sin(Q(L/2+t_f)cos\alpha)}{Q(L/2+t_f)cos\alpha} \\
    &+(\rho_r - \rho_s) V_t
     \frac{2J_1(Q(R'+t_r)sin\alpha)}{Q(R'+t_r)sin\alpha}
     \frac{sin(Q(L/2+t_f)cos\alpha)}{Q(L/2+t_f)cos\alpha}
    \bigg]
    \end{align*}

where

.. math::

    R'=\frac{R}{\sqrt{2}}\sqrt{(1+X_{core}^{2}) + (1-X_{core}^{2})cos(\psi)}


and $V_t = \pi.(R+t_r)(Xcore.R+t_r)^2.(L+2.t_f)$ is the total volume of
the bicelle, $V_c = \pi.Xcore.R^2.L$ the volume of the core,
$V_{c+f} = \pi.Xcore.R^2.(L+2.t_f)$ the volume of the core plus the volume
of the faces, $R$ is the radius of the core, $Xcore$ is the axial ratio of
the core, $L$ the length of the core, $t_f$ the thickness of the face, $t_r$
the thickness of the rim and $J_1$ the usual first order Bessel function.
The core has radii $R$ and $Xcore.R$ so is circular, as for the
core_shell_bicelle model, for $Xcore$ =1. Note that you may need to
limit the range of $Xcore$, especially if using the Monte-Carlo algorithm,
as setting radius to $R/Xcore$ and axial ratio to $1/Xcore$ gives an
equivalent solution!

The output of the 1D scattering intensity function for randomly oriented
bicelles is then given by integrating over all possible $\alpha$ and $\psi$.

For oriented bicelles the *theta*, *phi* and *psi* orientation parameters will
appear when fitting 2D data, see the :ref:`elliptical-cylinder` model
for further information.

.. figure:: img/elliptical_cylinder_angle_definition.png

    Definition of the angles for the oriented core_shell_bicelle_elliptical particles.

Model verified using Monte Carlo simulation for 1D and 2D scattering.

References
----------

.. [#Onsager1949] L. Onsager, *Ann. New York Acad. Sci.*, 51 (1949) 627-659

Authorship and Verification
----------------------------

* **Author:** Richard Heenan **Date:** December 14, 2016
* **Last Modified by:**  Richard Heenan **Date:** December 14, 2016
* **Last Reviewed by:**  Paul Kienzle **Date:** Feb 28, 2018
"""

import numpy as np
from numpy import cos, inf, pi, sin

name = "core_shell_bicelle_elliptical"
title = "Elliptical cylinder with a core-shell scattering length density profile.."
description = """
    core_shell_bicelle_elliptical
    Elliptical cylinder core, optional shell on the two flat faces, and shell of
    uniform thickness on its rim (extending around the end faces).
    Please see full documentation for equations and further details.
    Involves a double numerical integral around the ellipsoid diameter
    and the angle of the cylinder axis to Q.
    Compare also the core_shell_bicelle and elliptical_cylinder models.
      """
category = "shape:cylinder"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["radius",         "Ang",       30, [0, inf],    "volume",      "Cylinder core radius r_minor"],
    ["x_core",        "None",       3,  [0, inf],    "volume",      "Axial ratio of core, X = r_major/r_minor"],
    ["thick_rim",  "Ang",            8, [0, inf],    "volume",      "Rim shell thickness"],
    ["thick_face", "Ang",           14, [0, inf],    "volume",      "Cylinder face thickness"],
    ["length",         "Ang",       50, [0, inf],    "volume",      "Cylinder length"],
    ["sld_core",       "1e-6/Ang^2", 4, [-inf, inf], "sld",         "Cylinder core scattering length density"],
    ["sld_face",       "1e-6/Ang^2", 7, [-inf, inf], "sld",         "Cylinder face scattering length density"],
    ["sld_rim",        "1e-6/Ang^2", 1, [-inf, inf], "sld",         "Cylinder rim scattering length density"],
    ["sld_solvent",    "1e-6/Ang^2", 6, [-inf, inf], "sld",         "Solvent scattering length density"],
    ["theta",       "degrees",    90.0, [-360, 360], "orientation", "Cylinder axis to beam angle"],
    ["phi",         "degrees",    0,    [-360, 360], "orientation", "Rotation about beam"],
    ["psi",         "degrees",    0,    [-360, 360], "orientation", "Rotation about cylinder axis"]
    ]

# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sas_Si.c", "lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c",
          "core_shell_bicelle_elliptical.c"]
have_Fq = True
radius_effective_modes = [
    "equivalent cylinder excluded volume", "equivalent volume sphere",
    "outer rim average radius", "outer rim min radius",
    "outer max radius", "half outer thickness", "half diagonal",
    ]
has_shape_visualization = True

def create_shape_mesh(params, resolution=50):
    import numpy as np
    radius = params.get('radius', 30)  # r_minor
    x_core = params.get('x_core', 3)  # r_major/r_minor ratio
    thick_rim = params.get('thick_rim', 8)
    thick_face = params.get('thick_face', 14)
    length = params.get('length', 50)

    r_minor = radius
    r_major = radius * x_core

    # Outer dimensions
    outer_r_minor = r_minor + thick_rim
    outer_r_major = r_major + thick_rim
    outer_length = length + 2 * thick_face

    # Create core elliptical cylinder
    theta = np.linspace(0, 2*np.pi, resolution)
    z_core = np.linspace(-length/2, length/2, resolution//2)
    theta_core, z_core_mesh = np.meshgrid(theta, z_core)
    x_core_mesh = r_major * np.cos(theta_core)
    y_core_mesh = r_minor * np.sin(theta_core)

    # Create shell elliptical cylinder (outer surface)
    z_shell = np.linspace(-outer_length/2, outer_length/2, resolution//2)
    theta_shell, z_shell_mesh = np.meshgrid(theta, z_shell)
    x_shell = outer_r_major * np.cos(theta_shell)
    y_shell = outer_r_minor * np.sin(theta_shell)

    # Create end caps
    # Core end caps (elliptical)
    u = np.linspace(0, 1, resolution//4)
    theta_cap = np.linspace(0, 2*np.pi, resolution)
    u_mesh, theta_cap_mesh = np.meshgrid(u, theta_cap)

    x_cap_core = u_mesh * r_major * np.cos(theta_cap_mesh)
    y_cap_core = u_mesh * r_minor * np.sin(theta_cap_mesh)
    z_cap_core_top = np.full_like(x_cap_core, length/2)
    z_cap_core_bottom = np.full_like(x_cap_core, -length/2)

    # Shell end caps (elliptical)
    x_cap_shell = u_mesh * outer_r_major * np.cos(theta_cap_mesh)
    y_cap_shell = u_mesh * outer_r_minor * np.sin(theta_cap_mesh)
    z_cap_shell_top = np.full_like(x_cap_shell, outer_length/2)
    z_cap_shell_bottom = np.full_like(x_cap_shell, -outer_length/2)

    return {
        'core_cylinder': (x_core_mesh, y_core_mesh, z_core_mesh),
        'shell_cylinder': (x_shell, y_shell, z_shell_mesh),
        'core_cap_top': (x_cap_core, y_cap_core, z_cap_core_top),
        'core_cap_bottom': (x_cap_core, y_cap_core, z_cap_core_bottom),
        'shell_cap_top': (x_cap_shell, y_cap_shell, z_cap_shell_top),
        'shell_cap_bottom': (x_cap_shell, y_cap_shell, z_cap_shell_bottom),
    }

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    import numpy as np
    radius = params.get('radius', 30)
    x_core = params.get('x_core', 3)
    thick_rim = params.get('thick_rim', 8)
    thick_face = params.get('thick_face', 14)
    length = params.get('length', 50)

    r_minor = radius
    r_major = radius * x_core
    outer_r_minor = r_minor + thick_rim
    outer_r_major = r_major + thick_rim
    outer_length = length + 2 * thick_face

    # XY plane (top view) - nested ellipses
    theta = np.linspace(0, 2*np.pi, 100)

    # Core ellipse
    core_x = r_major * np.cos(theta)
    core_y = r_minor * np.sin(theta)

    # Shell ellipse
    shell_x = outer_r_major * np.cos(theta)
    shell_y = outer_r_minor * np.sin(theta)

    ax_xy.plot(shell_x, shell_y, 'r-', linewidth=2, label='Shell rim')
    ax_xy.fill(shell_x, shell_y, 'lightcoral', alpha=0.3)
    ax_xy.plot(core_x, core_y, 'b-', linewidth=2, label='Core')
    ax_xy.fill(core_x, core_y, 'lightblue', alpha=0.5)

    ax_xy.set_xlim(-outer_r_major*1.2, outer_r_major*1.2)
    ax_xy.set_ylim(-outer_r_minor*1.2, outer_r_minor*1.2)
    ax_xy.set_xlabel('X (Å) - Major axis')
    ax_xy.set_ylabel('Y (Å) - Minor axis')
    ax_xy.set_title('XY Cross-section (Top View)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)
    ax_xy.legend()

    # XZ plane (side view along major axis) - nested rectangles
    # Core rectangle (using major axis)
    core_rect_x = [-length/2, -length/2, length/2, length/2, -length/2]
    core_rect_z = [-r_major, r_major, r_major, -r_major, -r_major]

    # Full outer shell with face thickness
    shell_full_x = [-outer_length/2, -outer_length/2, outer_length/2, outer_length/2, -outer_length/2]
    shell_full_z = [-outer_r_major, outer_r_major, outer_r_major, -outer_r_major, -outer_r_major]

    ax_xz.plot(shell_full_x, shell_full_z, 'r-', linewidth=2, label='Shell (rim+face)')
    ax_xz.fill(shell_full_x, shell_full_z, 'lightcoral', alpha=0.3)
    ax_xz.plot(core_rect_x, core_rect_z, 'b-', linewidth=2, label='Core')
    ax_xz.fill(core_rect_x, core_rect_z, 'lightblue', alpha=0.5)

    ax_xz.set_xlim(-outer_length/2*1.2, outer_length/2*1.2)
    ax_xz.set_ylim(-outer_r_major*1.3, outer_r_major*1.3)
    ax_xz.set_xlabel('Z (Å)')
    ax_xz.set_ylabel('X (Å) - Major axis')
    ax_xz.set_title('XZ Cross-section (Major axis)')
    ax_xz.grid(True, alpha=0.3)
    ax_xz.legend()

    # Annotations
    ax_xz.text(0, -outer_r_major*1.15, f'L = {length:.0f} Å', ha='center', fontsize=9)
    ax_xz.text(0, outer_r_major*1.15, f'r_major = {r_major:.0f} Å', ha='center', fontsize=9)
    ax_xz.text(outer_length/2*0.7, outer_r_major*0.7, f't_rim = {thick_rim:.0f}', fontsize=8)
    ax_xz.text(length/2 + thick_face/2, 0, f't_face = {thick_face:.0f}', fontsize=8, rotation=90)

    # YZ plane (side view along minor axis)
    core_rect_y = [-r_minor, r_minor, r_minor, -r_minor, -r_minor]
    shell_body_y = [-outer_r_minor, outer_r_minor, outer_r_minor, -outer_r_minor, -outer_r_minor]

    ax_yz.plot(shell_full_x, shell_body_y, 'g-', linewidth=2, label='Shell (rim+face)')
    ax_yz.fill(shell_full_x, shell_body_y, 'lightgreen', alpha=0.3)
    ax_yz.plot(core_rect_x, core_rect_y, 'orange', linewidth=2, label='Core')
    ax_yz.fill(core_rect_x, core_rect_y, 'moccasin', alpha=0.5)

    ax_yz.set_xlim(-outer_length/2*1.2, outer_length/2*1.2)
    ax_yz.set_ylim(-outer_r_minor*1.3, outer_r_minor*1.3)
    ax_yz.set_xlabel('Z (Å)')
    ax_yz.set_ylabel('Y (Å) - Minor axis')
    ax_yz.set_title('YZ Cross-section (Minor axis)')
    ax_yz.grid(True, alpha=0.3)
    ax_yz.legend()

    ax_yz.text(0, outer_r_minor*1.15, f'r_minor = {r_minor:.0f} Å', ha='center', fontsize=9)
    ax_yz.text(0, -outer_r_minor*1.15, f'X_core = {x_core:.1f}', ha='center', fontsize=9)

def random():
    """Return a random parameter set for the model."""
    outer_major = 10**np.random.uniform(1, 4.7)
    outer_minor = 10**np.random.uniform(1, 4.7)
    # Use a distribution with a preference for thin shell or thin core,
    # limited by the minimum radius. Avoid core,shell radii < 1
    min_radius = min(outer_major, outer_minor)
    thick_rim = np.random.beta(0.5, 0.5)*(min_radius-2) + 1
    radius_major = outer_major - thick_rim
    radius_minor = outer_minor - thick_rim
    radius = radius_major
    x_core = radius_minor/radius_major
    outer_length = 10**np.random.uniform(1, 4.7)
    # Caps should be a small percentage of the total length, but at least one
    # angstrom long.  Since outer length >= 10, the following will suffice
    thick_face = 10**np.random.uniform(-np.log10(outer_length), -1)*outer_length
    length = outer_length - thick_face
    pars = dict(
        radius=radius,
        x_core=x_core,
        thick_rim=thick_rim,
        thick_face=thick_face,
        length=length
    )
    return pars


q = 0.1
# april 6 2017, rkh added a 2d unit test, NOT READY YET pull #890 branch assume correct!
qx = q*cos(pi/6.0)
qy = q*sin(pi/6.0)

tests = [
    #[{'radius': 30.0, 'x_core': 3.0,
    #  'thick_rim': 8.0, 'thick_face': 14.0, 'length': 50.0}, 'ER', 1],
    #[{'radius': 30.0, 'x_core': 3.0,
    #  'thick_rim': 8.0, 'thick_face': 14.0, 'length': 50.0}, 'VR', 1],

    [{'radius': 30.0, 'x_core': 3.0,
      'thick_rim': 8.0, 'thick_face': 14.0, 'length': 50.0,
      'sld_core': 4.0, 'sld_face': 7.0, 'sld_rim': 1.0,
      'sld_solvent': 6.0, 'background': 0.0},
     0.015, 286.540286],
    #[{'theta':80., 'phi':10.}, (qx, qy), 7.88866563001],
]

del qx, qy  # not necessary to delete, but cleaner
