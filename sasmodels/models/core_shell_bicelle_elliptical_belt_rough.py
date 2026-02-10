r"""
Definition
----------

This model provides the form factor for an elliptical cylinder with a
core-shell scattering length density profile [#Onsager1949]_.
Thus this is a variation
of the core-shell bicelle model, but with an elliptical cylinder for the core.
In this version the "rim" or "belt" does NOT extend the full length of
the particle, but has the same length as the core. Outer shells on the
rims and flat ends may be of different thicknesses and scattering length
densities. The form factor is normalized by the total particle volume.
This version includes an approximate "interfacial roughness".


.. figure:: img/core_shell_bicelle_belt_geometry.png

    Schematic cross-section of bicelle with belt. Note however that the model
    here calculates for rectangular, not curved, rims as shown below.

.. figure:: img/core_shell_bicelle_belt_parameters.png

   Cross section of model used here. Users will have
   to decide how to distribute "heads" and "tails" between the rim, face
   and core regions in order to estimate appropriate starting parameters.

Given the scattering length densities (sld) $\rho_c$, the core sld, $\rho_f$,
the face sld, $\rho_r$, the rim sld and $\rho_s$ the solvent sld, the
scattering length density variation along the bicelle axis is:

.. math::

    \rho(r) =
      \begin{cases}
      &\rho_c \text{ for } 0 \lt r \lt R;   -L/2 \lt z\lt L/2 \\[1.5ex]
      &\rho_f \text{ for } 0 \lt r \lt R;   -(L/2 +t_\text{face}) \lt z\lt -L/2;
      L/2 \lt z\lt (L/2+t_\text{face}) \\[1.5ex]
      &\rho_r\text{ for } R \lt r \lt R+t_\text{rim}; -L/2 \lt z\lt L/2
      \end{cases}

The form factor for the bicelle is calculated in cylindrical coordinates, where
$\alpha$ is the angle between the $Q$ vector and the cylinder axis, and $\psi$
is the angle for the ellipsoidal cross section core, to give:

.. math::

    I(Q,\alpha,\psi) = \frac{\text{scale}}{V_t}
        \cdot F(Q,\alpha, \psi)^2 \cdot \sin(\alpha)
        \cdot\exp\left\{ -\frac{1}{2}Q^2\sigma^2 \right\} + \text{background}

where a numerical integration of $F(Q,\alpha, \psi)^2\sin(\alpha)$ is
carried out over $\alpha$ and $\psi$ for:

.. math::

    F(Q,\alpha,\psi) = &\bigg[
      (\rho_c -\rho_r - \rho_f + \rho_s) V_c
      \frac{2J_1(QR'\sin \alpha)}{QR'\sin\alpha}
      \frac{\sin(QL\cos\alpha/2)}{Q(L/2)\cos\alpha} \\
    &+(\rho_f - \rho_s) V_{c+f}
      \frac{2J_1(QR'\sin\alpha)}{QR'\sin\alpha}
      \frac{\sin(Q(L/2+t_f)\cos\alpha)}{Q(L/2+t_f)\cos\alpha} \\
    &+(\rho_r - \rho_s) V_{c+r}
      \frac{2J_1(Q(R'+t_r)\sin\alpha)}{Q(R'+t_r)\sin\alpha}
      \frac{\sin(Q(L/2)\cos\alpha)}{Q(L/2)\cos\alpha}
    \bigg]

where

.. math::

    R' = \frac{R}{\sqrt{2}}
        \sqrt{(1+X_\text{core}^{2}) + (1-X_\text{core}^{2})\cos(\psi)}


and $V_t = \pi (R+t_r)(X_\text{core} R+t_r) L + 2 \pi X_\text{core} R^2 t_f$ is
the total volume of the bicelle, $V_c = \pi X_\text{core} R^2 L$ the volume of
the core, $V_{c+f} = \pi X_\text{core} R^2 (L+2 t_f)$ the volume of the core
plus the volume of the faces, $V_{c+r} = \pi (R+t_r)(X_\text{core} R+t_r) L$
the volume of the core plus the rim, $R$ is the radius of the core,
$X_\text{core}$ is the axial ratio of the core, $L$ the length of the core,
$t_f$ the thickness of the face, $t_r$ the thickness of the rim and $J_1$ the
usual first order Bessel function. The core has radii $R$ and $X_\text{core} R$
so is circular, as for the core_shell_bicelle model, for $X_\text{core}=1$.
Note that you may need to limit the range of $X_\text{core}$, especially if
using the Monte-Carlo algorithm, as setting radius to $R/X_\text{core}$ and
axial ratio to $1/X_\text{core}$ gives an equivalent solution!

An approximation for the effects of "Gaussian interfacial roughness" $\sigma$
is included, by multiplying $I(Q)$ by
$\exp\left \{ -\frac{1}{2}Q^2\sigma^2 \right \}$. This applies, in some way, to
all interfaces in the model not just the external ones. (Note that for a one
dimensional system convolution of the scattering length density profile with
a Gaussian of standard deviation $\sigma$ does exactly this multiplication.)
Leave $\sigma$ set to zero for the usual sharp interfaces.

The output of the 1D scattering intensity function for randomly oriented
bicelles is then given by integrating over all possible $\alpha$ and $\psi$.

For oriented bicelles the *theta*, *phi* and *psi* orientation parameters
will appear when fitting 2D data, for further details of the calculation
and angular dispersions  see :ref:`orientation` .

.. figure:: img/elliptical_cylinder_angle_definition.png

    Definition of the angles for the oriented core_shell_bicelle_elliptical
    particles.



References
----------

.. [#Onsager1949] L. Onsager, *Ann. New York Acad. Sci.*, 51 (1949) 627-659

Authorship and Verification
----------------------------

* **Author:** Richard Heenan **Date:** October 5, 2017
* **Last Modified by:**  Richard Heenan new 2d orientation **Date:** October 5, 2017
* **Last Reviewed by:**  Richard Heenan 2d calc seems agree with 1d **Date:** Nov 2, 2017
"""

from numpy import cos, inf, pi, sin

name = "core_shell_bicelle_elliptical_belt_rough"
title = "Elliptical cylinder with a core-shell scattering length density profile.."
description = """
    core_shell_bicelle_elliptical_belt_rough
    Elliptical cylinder core, optional shell on the two flat faces, and "belt" shell of
    uniform thickness on its rim (in this case NOT extending around the end faces).
    with approximate interfacial roughness.
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
    ["thick_rim",  "Ang",            8, [0, inf],    "volume",      "Rim or belt shell thickness"],
    ["thick_face", "Ang",           14, [0, inf],    "volume",      "Cylinder face thickness"],
    ["length",         "Ang",       50, [0, inf],    "volume",      "Cylinder length"],
    ["sld_core",       "1e-6/Ang^2", 4, [-inf, inf], "sld",         "Cylinder core scattering length density"],
    ["sld_face",       "1e-6/Ang^2", 7, [-inf, inf], "sld",         "Cylinder face scattering length density"],
    ["sld_rim",        "1e-6/Ang^2", 1, [-inf, inf], "sld",         "Cylinder rim scattering length density"],
    ["sld_solvent",    "1e-6/Ang^2", 6, [-inf, inf], "sld",         "Solvent scattering length density"],
    ["sigma",       "Ang",        0,    [0, inf],    "",            "Interfacial roughness"],
    ["theta",       "degrees",    90.0, [-360, 360], "orientation", "Cylinder axis to beam angle"],
    ["phi",         "degrees",    0,    [-360, 360], "orientation", "Rotation about beam"],
    ["psi",         "degrees",    0,    [-360, 360], "orientation", "Rotation about cylinder axis"],
    ]

# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sas_Si.c", "lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c",
          "core_shell_bicelle_elliptical_belt_rough.c"]
has_shape_visualization = True
have_Fq = True
radius_effective_modes = [
    "equivalent cylinder excluded volume", "equivalent volume sphere",
    "outer rim average radius", "outer rim min radius",
    "outer max radius", "half outer thickness", "half diagonal",
    ]

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for elliptical bicelle with belt visualization."""
    import numpy as np
    radius = params.get('radius', 30)
    x_core = params.get('x_core', 3)
    thick_rim = params.get('thick_rim', 8)
    thick_face = params.get('thick_face', 14)
    length = params.get('length', 50)

    r_minor = radius
    r_major = radius * x_core
    rim_r_minor = r_minor + thick_rim
    rim_r_major = r_major + thick_rim

    theta = np.linspace(0, 2*np.pi, resolution)

    # Core elliptical cylinder
    z_core = np.linspace(-length/2, length/2, resolution//2)
    theta_core, z_core_mesh = np.meshgrid(theta, z_core)
    x_core_mesh = r_major * np.cos(theta_core)
    y_core_mesh = r_minor * np.sin(theta_core)

    # Rim/belt cylinder (same height as core, wider cross-section)
    theta_rim, z_rim_mesh = np.meshgrid(theta, z_core)
    x_rim = rim_r_major * np.cos(theta_rim)
    y_rim = rim_r_minor * np.sin(theta_rim)

    # Face shell top surface (core cross-section, extends above core)
    z_face_top = np.linspace(length/2, length/2 + thick_face, max(resolution//8, 2))
    theta_ft, z_ft_mesh = np.meshgrid(theta, z_face_top)
    x_ft = r_major * np.cos(theta_ft)
    y_ft = r_minor * np.sin(theta_ft)

    # Face shell bottom surface
    z_face_bot = np.linspace(-length/2 - thick_face, -length/2, max(resolution//8, 2))
    theta_fb, z_fb_mesh = np.meshgrid(theta, z_face_bot)
    x_fb = r_major * np.cos(theta_fb)
    y_fb = r_minor * np.sin(theta_fb)

    # End caps (elliptical disks)
    u = np.linspace(0, 1, resolution//4)
    u_mesh, theta_cap = np.meshgrid(u, theta)

    # Face cap top
    x_cap_top = u_mesh * r_major * np.cos(theta_cap)
    y_cap_top = u_mesh * r_minor * np.sin(theta_cap)
    z_cap_top = np.full_like(x_cap_top, length/2 + thick_face)

    # Face cap bottom
    z_cap_bot = np.full_like(x_cap_top, -length/2 - thick_face)

    # Rim annular caps at z = +/- length/2
    r_frac = np.linspace(0, 1, resolution//4)
    rf_mesh, theta_ann = np.meshgrid(r_frac, theta)
    x_ann = (r_major + rf_mesh * thick_rim) * np.cos(theta_ann)
    y_ann = (r_minor + rf_mesh * thick_rim) * np.sin(theta_ann)
    z_ann_top = np.full_like(x_ann, length/2)
    z_ann_bot = np.full_like(x_ann, -length/2)

    return {
        'core_cylinder': (x_core_mesh, y_core_mesh, z_core_mesh),
        'rim_cylinder': (x_rim, y_rim, z_rim_mesh),
        'face_top': (x_ft, y_ft, z_ft_mesh),
        'face_bottom': (x_fb, y_fb, z_fb_mesh),
        'face_cap_top': (x_cap_top, y_cap_top, z_cap_top),
        'face_cap_bottom': (x_cap_top, y_cap_top, z_cap_bot),
        'rim_cap_top': (x_ann, y_ann, z_ann_top),
        'rim_cap_bottom': (x_ann, y_ann, z_ann_bot),
    }

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of the elliptical bicelle with belt."""
    import numpy as np
    radius = params.get('radius', 30)
    x_core = params.get('x_core', 3)
    thick_rim = params.get('thick_rim', 8)
    thick_face = params.get('thick_face', 14)
    length = params.get('length', 50)

    r_minor = radius
    r_major = radius * x_core
    rim_r_minor = r_minor + thick_rim
    rim_r_major = r_major + thick_rim
    half_L = length / 2

    # --- XY plane (top view, z=0): core ellipse + rim annulus ---
    theta = np.linspace(0, 2*np.pi, 100)
    core_x = r_major * np.cos(theta)
    core_y = r_minor * np.sin(theta)
    rim_x = rim_r_major * np.cos(theta)
    rim_y = rim_r_minor * np.sin(theta)

    ax_xy.plot(rim_x, rim_y, 'r-', linewidth=2, label='Rim/belt')
    ax_xy.fill(rim_x, rim_y, 'lightcoral', alpha=0.3)
    ax_xy.plot(core_x, core_y, 'b-', linewidth=2, label='Core')
    ax_xy.fill(core_x, core_y, 'lightblue', alpha=0.5)
    ax_xy.set_xlim(-rim_r_major*1.2, rim_r_major*1.2)
    ax_xy.set_ylim(-rim_r_minor*1.2, rim_r_minor*1.2)
    ax_xy.set_xlabel('X (Å) - Major axis')
    ax_xy.set_ylabel('Y (Å) - Minor axis')
    ax_xy.set_title('XY Cross-section (Top View)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)
    ax_xy.legend(fontsize=8)

    # --- XZ plane (side view, major axis): cross / plus shape ---
    # Rim rectangle (wider, core height only)
    rim_rect_x = [-half_L, -half_L, half_L, half_L, -half_L]
    rim_rect_z = [-rim_r_major, rim_r_major, rim_r_major, -rim_r_major, -rim_r_major]
    # Face rectangle (core width, full height)
    face_rect_x = [-(half_L + thick_face), -(half_L + thick_face),
                   half_L + thick_face, half_L + thick_face, -(half_L + thick_face)]
    face_rect_z = [-r_major, r_major, r_major, -r_major, -r_major]
    # Core rectangle
    core_rect_x = [-half_L, -half_L, half_L, half_L, -half_L]
    core_rect_z = [-r_major, r_major, r_major, -r_major, -r_major]

    ax_xz.fill(rim_rect_x, rim_rect_z, 'lightcoral', alpha=0.3, label='Rim/belt')
    ax_xz.plot(rim_rect_x, rim_rect_z, 'r-', linewidth=2)
    ax_xz.fill(face_rect_x, face_rect_z, 'lightyellow', alpha=0.5, label='Face')
    ax_xz.plot(face_rect_x, face_rect_z, color='goldenrod', linewidth=2)
    ax_xz.fill(core_rect_x, core_rect_z, 'lightblue', alpha=0.5, label='Core')
    ax_xz.plot(core_rect_x, core_rect_z, 'b-', linewidth=2)
    max_x = (half_L + thick_face) * 1.2
    max_z = rim_r_major * 1.2
    ax_xz.set_xlim(-max_x, max_x)
    ax_xz.set_ylim(-max_z, max_z)
    ax_xz.set_xlabel('Z (Å)')
    ax_xz.set_ylabel('X (Å) - Major axis')
    ax_xz.set_title('XZ Cross-section (Major axis)')
    ax_xz.grid(True, alpha=0.3)
    ax_xz.legend(fontsize=8)

    # --- YZ plane (side view, minor axis): cross / plus shape ---
    rim_rect_y = [-rim_r_minor, rim_r_minor, rim_r_minor, -rim_r_minor, -rim_r_minor]
    face_rect_y = [-r_minor, r_minor, r_minor, -r_minor, -r_minor]
    core_rect_y = [-r_minor, r_minor, r_minor, -r_minor, -r_minor]

    ax_yz.fill(rim_rect_x, rim_rect_y, 'lightgreen', alpha=0.3, label='Rim/belt')
    ax_yz.plot(rim_rect_x, rim_rect_y, 'g-', linewidth=2)
    ax_yz.fill(face_rect_x, face_rect_y, 'lightyellow', alpha=0.5, label='Face')
    ax_yz.plot(face_rect_x, face_rect_y, color='goldenrod', linewidth=2)
    ax_yz.fill(core_rect_x, core_rect_y, 'lightblue', alpha=0.5, label='Core')
    ax_yz.plot(core_rect_x, core_rect_y, 'b-', linewidth=2)
    max_y = rim_r_minor * 1.2
    ax_yz.set_xlim(-max_x, max_x)
    ax_yz.set_ylim(-max_y, max_y)
    ax_yz.set_xlabel('Z (Å)')
    ax_yz.set_ylabel('Y (Å) - Minor axis')
    ax_yz.set_title('YZ Cross-section (Minor axis)')
    ax_yz.grid(True, alpha=0.3)
    ax_yz.legend(fontsize=8)

# TODO: No random() for core-shell bicelle elliptical belt rough

q = 0.1
# april 6 2017, rkh added a 2d unit test, NOT READY YET pull #890 branch assume correct!
qx = q*cos(pi/6.0)
qy = q*sin(pi/6.0)

tests = [
    #[{'radius': 30.0, 'x_core': 3.0, 'thick_rim':8.0, 'thick_face':14.0, 'length':50.0}, 'ER', 1],
    #[{'radius': 30.0, 'x_core': 3.0, 'thick_rim':8.0, 'thick_face':14.0, 'length':50.0}, 'VR', 1],

    [{'radius': 30.0, 'x_core': 3.0, 'thick_rim': 8.0, 'thick_face': 14.0,
      'length': 50.0, 'sld_core': 4.0, 'sld_face': 7.0, 'sld_rim': 1.0,
      'sld_solvent': 6.0, 'background': 0.0},
     0.015, 189.328],
    #[{'theta':80., 'phi':10.}, (qx, qy), 7.88866563001 ],
]

del qx, qy  # not necessary to delete, but cleaner
