r"""
Definitions
-----------

Calculates the scattering from a cylinder with spherical section end-caps.
Like :ref:`barbell`, this is a sphereocylinder with end caps that have a
radius larger than that of the cylinder, but with the center of the end cap
radius lying within the cylinder. This model simply becomes a convex
lens when the length of the cylinder $L=0$. See the diagram for the details
of the geometry and restrictions on parameter values.

.. figure:: img/capped_cylinder_geometry.jpg

    Capped cylinder geometry, where $r$ is *radius*, $R$ is *radius_cap* and
    $L$ is *length*. Since the end cap radius $R \geq r$ and by definition
    for this geometry $h \le 0$, $h$ is then defined by $r$ and $R$ as
    $h = -\sqrt{R^2 - r^2}$

The scattered intensity $I(q)$ is calculated as

.. math::

    I(q) = \frac{\Delta \rho^2}{V} \left<A^2(q,\alpha).sin(\alpha)\right>

where the amplitude $A(q,\alpha)$ with the rod axis at angle $\alpha$ to $q$
is given as

.. math::

    A(q) =&\ \pi r^2L
        \frac{\sin\left(\tfrac12 qL\cos\alpha\right)}
            {\tfrac12 qL\cos\alpha}
        \frac{2 J_1(qr\sin\alpha)}{qr\sin\alpha} \\
        &\ + 4 \pi R^3 \int_{-h/R}^1 dt
        \cos\left[ q\cos\alpha
            \left(Rt + h + {\tfrac12} L\right)\right]
        \times (1-t^2)
        \frac{J_1\left[qR\sin\alpha \left(1-t^2\right)^{1/2}\right]}
             {qR\sin\alpha \left(1-t^2\right)^{1/2}}

The $\left<\ldots\right>$ brackets denote an average of the structure over
all orientations. $\left< A^2(q)\right>$ is then the form factor, $P(q)$.
The scale factor is equivalent to the volume fraction of cylinders, each of
volume, $V$. Contrast $\Delta\rho$ is the difference of scattering length
densities of the cylinder and the surrounding solvent.

The volume of the capped cylinder is (with $h$ as a positive value here)

.. math::

    V = \pi r_c^2 L + 2\pi\left(\tfrac23R^3 + R^2h - \tfrac13h^3\right)

and its radius of gyration is

.. math::

    R_g^2 =&\ \left[ \tfrac{12}{5}R^5
        + R^4\left(6h+\tfrac32 L\right)
        + R^3\left(4h^2 + L^2 + 4Lh\right)
        + R^2\left(3Lh^2 + \tfrac32 L^2h\right) \right. \\
        &\ \left. + \tfrac25 h^5 - \tfrac12 Lh^4 - \tfrac12 L^2h^3
        + \tfrac14 L^3r^2 + \tfrac32 Lr^4 \right]
        \left( 4R^3 + 6R^2h - 2h^3 + 3r^2L \right)^{-1}


.. note::

    The requirement that $R \geq r$ is not enforced in the model!
    It is up to you to restrict this during analysis.

The 2D scattering intensity is calculated similar to the 2D cylinder model.

.. figure:: img/cylinder_angle_definition.png

    Definition of the angles for oriented 2D cylinders.


References
----------

#. H Kaya, *J. Appl. Cryst.*, 37 (2004) 223-230

#. H Kaya and N R deSouza, *J. Appl. Cryst.*, 37 (2004) 508-509
   (addenda and errata)

#. L. Onsager, *Ann. New York Acad. Sci.*, 51 (1949) 627-659

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Paul Butler **Date:** September 30, 2016
* **Last Reviewed by:** Richard Heenan **Date:** January 4, 2017
"""

import numpy as np
from numpy import cos, inf, pi, sin

name = "capped_cylinder"
title = "Right circular cylinder with spherical end caps and uniform SLD"
description = """That is, a sphereocylinder
    with end caps that have a radius larger than
    that of the cylinder and the center of the
    end cap radius lies within the cylinder.
    Note: As the length of cylinder -->0,
    it becomes a Convex Lens.
    It must be that radius <(=) radius_cap.
    [Parameters];
    scale: volume fraction of spheres,
    background:incoherent background,
    radius: radius of the cylinder,
    length: length of the cylinder,
    radius_cap: radius of the semi-spherical cap,
    sld: SLD of the capped cylinder,
    sld_solvent: SLD of the solvent.
"""
category = "shape:cylinder"
# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld",         "1e-6/Ang^2", 4, [-inf, inf], "sld",    "Cylinder scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld",    "Solvent scattering length density"],
              ["radius",      "Ang",       20, [0, inf],    "volume", "Cylinder radius"],

              # TODO: use an expression for cap radius with fixed bounds.
              # The current form requires cap radius R bigger than cylinder radius r.
              # Could instead use R/r in [1,inf], r/R in [0,1], or the angle between
              # cylinder and cap in [0,90].  The problem is similar for the barbell
              # model.  Propose r/R in [0,1] in both cases, with the model specifying
              # cylinder radius in the capped cylinder model and sphere radius in the
              # barbell model.  This leads to the natural value of zero for no cap
              # in the capped cylinder, and zero for no bar in the barbell model.  In
              # both models, one would be a pill.
              ["radius_cap", "Ang",     20, [0, inf],    "volume", "Cap radius"],
              ["length",     "Ang",    400, [0, inf],    "volume", "Cylinder length"],
              ["theta",      "degrees", 60, [-360, 360], "orientation", "cylinder axis to beam angle"],
              ["phi",        "degrees", 60, [-360, 360], "orientation", "rotation about beam"],
             ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c", "capped_cylinder.c"]
valid = "radius_cap >= radius"
have_Fq = True
radius_effective_modes = [
    "equivalent cylinder excluded volume", "equivalent volume sphere",
    "radius", "half length", "half total length",
    ]
has_shape_visualization = True

def create_shape_mesh(params, resolution=50):
    import numpy as np
    radius = params.get('radius', 20)
    radius_cap = params.get('radius_cap', 25)
    length = params.get('length', 400)

    if radius_cap < radius:
        raise ValueError(f"Cap radius ({radius_cap}) must be >= cylinder radius ({radius})")

    # Calculate cap geometry
    h = np.sqrt(radius_cap**2 - radius**2)

    # Create cylinder body
    theta = np.linspace(0, 2*np.pi, resolution)
    z_cyl = np.linspace(-length/2, length/2, resolution//2)
    theta_cyl, z_cyl_mesh = np.meshgrid(theta, z_cyl)
    x_cyl = radius * np.cos(theta_cyl)
    y_cyl = radius * np.sin(theta_cyl)

    # Create spherical caps
    phi_max = np.arccos(h / radius_cap)
    phi = np.linspace(0, phi_max, resolution//4)
    phi_mesh, theta_mesh = np.meshgrid(phi, theta)

    # Top cap
    x_cap_top = radius_cap * np.sin(phi_mesh) * np.cos(theta_mesh)
    y_cap_top = radius_cap * np.sin(phi_mesh) * np.sin(theta_mesh)
    z_cap_top = length/2 - h + radius_cap * np.cos(phi_mesh)

    # Bottom cap
    x_cap_bottom = radius_cap * np.sin(phi_mesh) * np.cos(theta_mesh)
    y_cap_bottom = radius_cap * np.sin(phi_mesh) * np.sin(theta_mesh)
    z_cap_bottom = -length/2 + h - radius_cap * np.cos(phi_mesh)

    return {
        'cylinder': (x_cyl, y_cyl, z_cyl_mesh),
        'cap_top': (x_cap_top, y_cap_top, z_cap_top),
        'cap_bottom': (x_cap_bottom, y_cap_bottom, z_cap_bottom)
    }

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot cross-sections matching SasView documentation figure."""
    import numpy as np

    radius = params.get('radius', 20)  # r in docs
    radius_cap = params.get('radius_cap', 25)  # R in docs
    length = params.get('length', 400)  # L in docs

    if radius_cap < radius:
        return  # Skip if invalid parameters

    # h = sqrt(R^2 - r^2), cap centers are INSIDE cylinder at z = ±(L/2 - h)
    h = np.sqrt(radius_cap**2 - radius**2)

    theta = np.linspace(0, 2*np.pi, 100)

    # XY plane (top view) - circle (cylinder cross-section)
    circle_x = radius * np.cos(theta)
    circle_y = radius * np.sin(theta)
    ax_xy.fill(circle_x, circle_y, 'lightblue', alpha=0.5)
    ax_xy.plot(circle_x, circle_y, 'b-', linewidth=2, label=f'r={radius:.0f}Å')

    ax_xy.set_xlim(-radius*1.4, radius*1.4)
    ax_xy.set_ylim(-radius*1.4, radius*1.4)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title('XY Cross-section')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)
    ax_xy.legend(fontsize=9)

    # XZ plane (side view) - matching documentation figure
    # Cylinder body
    cyl_z = np.array([-length/2, length/2, length/2, -length/2, -length/2])
    cyl_r = np.array([-radius, -radius, radius, radius, -radius])
    ax_xz.fill(cyl_z, cyl_r, 'lightblue', alpha=0.5)
    ax_xz.plot(cyl_z, cyl_r, 'b-', linewidth=2)

    # Right cap (positive z) - center at z = L/2 - h (inside cylinder)
    cap_center_right = length/2 - h
    angle_junction = np.arcsin(radius / radius_cap)
    cap_angles = np.linspace(-angle_junction, angle_junction, 50)
    # Cap extends from cylinder end outward
    cap_z_right = cap_center_right + radius_cap * np.cos(cap_angles)
    cap_r_right = radius_cap * np.sin(cap_angles)

    # Create closed polygon for right cap
    right_cap_z = np.concatenate([[length/2], cap_z_right, [length/2]])
    right_cap_r = np.concatenate([[radius], cap_r_right, [-radius]])
    ax_xz.fill(right_cap_z, right_cap_r, 'lightcoral', alpha=0.5)
    ax_xz.plot(cap_z_right, cap_r_right, 'r-', linewidth=2, label='Caps')

    # Left cap (negative z) - center at z = -L/2 + h (inside cylinder)
    cap_center_left = -length/2 + h
    cap_z_left = cap_center_left - radius_cap * np.cos(cap_angles)
    cap_r_left = radius_cap * np.sin(cap_angles)

    left_cap_z = np.concatenate([[-length/2], cap_z_left, [-length/2]])
    left_cap_r = np.concatenate([[radius], cap_r_left, [-radius]])
    ax_xz.fill(left_cap_z, left_cap_r, 'lightcoral', alpha=0.5)
    ax_xz.plot(cap_z_left, cap_r_left, 'r-', linewidth=2)

    # Show cap centers (inside cylinder per documentation)
    ax_xz.plot(cap_center_right, 0, 'ko', markersize=4)
    ax_xz.plot(cap_center_left, 0, 'ko', markersize=4)

    # Add dimension annotations
    ax_xz.annotate('', xy=(length/2, 0), xytext=(-length/2, 0),
                  arrowprops=dict(arrowstyle='<->', color='black'))
    ax_xz.text(0, -radius*1.3, f'L={length:.0f}Å', ha='center', fontsize=9)

    ax_xz.set_xlim(-length/2 - radius_cap*0.3, length/2 + radius_cap*0.3)
    ax_xz.set_ylim(-radius*1.5, radius*1.5)
    ax_xz.set_xlabel('Z (Å) - along axis')
    ax_xz.set_ylabel('Radial (Å)')
    ax_xz.set_title('Side View (caps inside)')
    ax_xz.grid(True, alpha=0.3)
    ax_xz.legend(fontsize=9)

    ax_yz.fill(left_cap_z, left_cap_r, 'moccasin', alpha=0.5)
    ax_yz.plot(cap_z_left, cap_r_left, 'orange', linewidth=2)

    ax_yz.set_xlim(-length/2 - radius_cap*0.3, length/2 + radius_cap*0.3)
    ax_yz.set_ylim(-radius*1.5, radius*1.5)
    ax_yz.set_xlabel('Z (Å)')
    ax_yz.set_ylabel('Radial (Å)')
    ax_yz.set_title('YZ Cross-section')
    ax_yz.grid(True, alpha=0.3)

def random():
    """Return a random parameter set for the model."""
    # TODO: increase volume range once problem with bell radius is fixed
    # The issue is that bell radii of more than about 200 fail at high q
    volume = 10**np.random.uniform(7, 9)
    bar_volume = 10**np.random.uniform(-4, -1)*volume
    bell_volume = volume - bar_volume
    bell_radius = (bell_volume/6)**0.3333  # approximate
    min_bar = bar_volume/np.pi/bell_radius**2
    bar_length = 10**np.random.uniform(0, 3)*min_bar
    bar_radius = np.sqrt(bar_volume/bar_length/np.pi)
    if bar_radius > bell_radius:
        bell_radius, bar_radius = bar_radius, bell_radius
    pars = dict(
        #background=0,
        radius_cap=bell_radius,
        radius=bar_radius,
        length=bar_length,
    )
    return pars


q = 0.1
# 2017-04-06: rkh add unit tests, NOT compared with any other calc method, assume correct!
# 2019-05-17: pak added barbell/capped cylinder to realspace sampling tests
qx = q*cos(pi/6.0)
qy = q*sin(pi/6.0)
tests = [
    [{}, 0.075, 26.0698570695],
    [{'theta':80., 'phi':10.}, (qx, qy), 0.561811990502],
]
del qx, qy  # not necessary to delete, but cleaner
