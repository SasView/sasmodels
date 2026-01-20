r"""
Definition
----------

Calculates the scattering from a barbell-shaped cylinder.  Like
:ref:`capped-cylinder`, this is a spherocylinder with spherical end
caps that have a radius larger than that of the cylinder, but with the center
of the end cap radius lying outside of the cylinder. See the diagram for
the details of the geometry and restrictions on parameter values.

.. figure:: img/barbell_geometry.jpg

    Barbell geometry, where $r$ is *radius*, $R$ is *radius_bell* and
    $L$ is *length*. Since the end cap radius $R \geq r$ and by definition
    for this geometry $h \ge 0$, $h$ is then defined by $r$ and $R$ as
    $h = \sqrt{R^2 - r^2}$

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
all orientations. $\left<A^2(q,\alpha)\right>$ is then the form factor, $P(q)$.
The scale factor is equivalent to the volume fraction of cylinders, each of
volume, $V$. Contrast $\Delta\rho$ is the difference of scattering length
densities of the cylinder and the surrounding solvent.

The volume of the barbell is

.. math::

    V = \pi r_c^2 L + 2\pi\left(\tfrac23R^3 + R^2h-\tfrac13h^3\right)

and its radius of gyration is

.. math::

    R_g^2 =&\ \left[ \tfrac{12}{5}R^4
        + R^3\left(3L + \tfrac{18}{5} h\right)
        + R^2\left(L^2 + Lh + \tfrac25 h^2\right)
        + R\left(\tfrac14 L^3 + \tfrac12 L^2h - Lh^2\right) \right. \\
        &\ \left. + Lh^4 - \tfrac12 L^2h^3 - \tfrac14 L^3h + \tfrac25 h^4\right]
        \left( 4R^2 + 3LR + 2Rh - 3Lh - 2h^2\right)^{-1}

.. note::
    The requirement that $R \geq r$ is not enforced in the model! It is
    up to you to restrict this during analysis.

The 2D scattering intensity is calculated similar to the 2D cylinder model.

.. figure:: img/cylinder_angle_definition.png

    Definition of the angles for oriented 2D barbells.


References
----------

#. H Kaya, *J. Appl. Cryst.*, 37 (2004) 223-230

#. H Kaya and N R deSouza, *J. Appl. Cryst.*, 37 (2004) 508-509
   (addenda and errata)

#. L. Onsager, *Ann. New York Acad. Sci.*, 51 (1949) 627-659

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Paul Butler **Date:** March 20, 2016
* **Last Reviewed by:** Richard Heenan **Date:** January 4, 2017
"""

import numpy as np
from numpy import cos, inf, pi, sin

name = "barbell"
title = "Cylinder with spherical end caps"
description = """
    Calculates the scattering from a barbell-shaped cylinder.
    That is a sphereocylinder with spherical end caps that have a radius larger
    than that of the cylinder and the center of the end cap radius lies outside
    of the cylinder.
    Note: As the length of cylinder(bar) -->0,it becomes a dumbbell. And when
    rad_bar = rad_bell, it is a spherocylinder.
    It must be that rad_bar <(=) rad_bell.
"""
category = "shape:cylinder"
# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["sld",         "1e-6/Ang^2",   4, [-inf, inf], "sld",         "Barbell scattering length density"],
    ["sld_solvent", "1e-6/Ang^2",   1, [-inf, inf], "sld",         "Solvent scattering length density"],
    ["radius_bell", "Ang",         40, [0, inf],    "volume",      "Spherical bell radius"],
    ["radius",      "Ang",         20, [0, inf],    "volume",      "Cylindrical bar radius"],
    ["length",      "Ang",        400, [0, inf],    "volume",      "Cylinder bar length"],
    ["theta",       "degrees",     60, [-360, 360], "orientation", "Barbell axis to beam angle"],
    ["phi",         "degrees",     60, [-360, 360], "orientation", "Rotation about beam"],
    ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c", "barbell.c"]
has_shape_visualization = True

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for barbell visualization."""
    import numpy as np
    radius = params.get('radius', 20)
    radius_bell = params.get('radius_bell', 40)
    length = params.get('length', 400)

    if radius_bell < radius:
        radius_bell = radius  # Ensure valid geometry

    # Height where bell meets cylinder
    h = np.sqrt(radius_bell**2 - radius**2)

    # Create cylinder body
    theta = np.linspace(0, 2*np.pi, resolution)
    z_cyl = np.linspace(-length/2, length/2, resolution//2)
    theta_cyl, z_cyl_mesh = np.meshgrid(theta, z_cyl)
    x_cyl = radius * np.cos(theta_cyl)
    y_cyl = radius * np.sin(theta_cyl)

    # Create spherical bells (larger than cylinder caps)
    # The bell center is at distance h inside the cylinder end
    phi_max = np.arcsin(radius / radius_bell)
    phi = np.linspace(0, phi_max, resolution//4)
    phi_mesh, theta_mesh = np.meshgrid(phi, theta)

    # Top bell
    x_bell_top = radius_bell * np.sin(phi_mesh) * np.cos(theta_mesh)
    y_bell_top = radius_bell * np.sin(phi_mesh) * np.sin(theta_mesh)
    z_bell_top = length/2 + h + radius_bell * np.cos(phi_mesh) - radius_bell

    # Bottom bell
    x_bell_bottom = radius_bell * np.sin(phi_mesh) * np.cos(theta_mesh)
    y_bell_bottom = radius_bell * np.sin(phi_mesh) * np.sin(theta_mesh)
    z_bell_bottom = -length/2 - h - radius_bell * np.cos(phi_mesh) + radius_bell

    return {
        'cylinder': (x_cyl, y_cyl, z_cyl_mesh),
        'bell_top': (x_bell_top, y_bell_top, z_bell_top),
        'bell_bottom': (x_bell_bottom, y_bell_bottom, z_bell_bottom)
    }

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of the barbell matching SasView documentation."""
    import numpy as np

    radius = params.get('radius', 20)  # r in docs
    radius_bell = params.get('radius_bell', 40)  # R in docs
    length = params.get('length', 400)  # L in docs

    if radius_bell < radius:
        radius_bell = radius

    # h = sqrt(R^2 - r^2) per documentation
    h = np.sqrt(radius_bell**2 - radius**2)

    theta = np.linspace(0, 2*np.pi, 100)

    # XY plane - circle (cylinder cross-section at center)
    circle_x = radius * np.cos(theta)
    circle_y = radius * np.sin(theta)
    ax_xy.fill(circle_x, circle_y, 'lightblue', alpha=0.5)
    ax_xy.plot(circle_x, circle_y, 'b-', linewidth=2, label=f'r={radius:.0f}Å')

    # Show bell outline for reference
    bell_x = radius_bell * np.cos(theta)
    bell_y = radius_bell * np.sin(theta)
    ax_xy.plot(bell_x, bell_y, 'r--', linewidth=1.5, alpha=0.6, label=f'R={radius_bell:.0f}Å')

    ax_xy.set_xlim(-radius_bell*1.4, radius_bell*1.4)
    ax_xy.set_ylim(-radius_bell*1.4, radius_bell*1.4)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title('XY Cross-section')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)
    ax_xy.legend(fontsize=9)

    # XZ plane (side view) - matching documentation figure
    # Draw cylinder body
    cyl_z = np.array([-length/2, length/2, length/2, -length/2, -length/2])
    cyl_r = np.array([-radius, -radius, radius, radius, -radius])
    ax_xz.fill(cyl_z, cyl_r, 'lightblue', alpha=0.5)
    ax_xz.plot(cyl_z, cyl_r, 'b-', linewidth=2)

    # Right bell (positive z) - sphere center at z = L/2 + h
    bell_center_right = length/2 + h
    # Arc from angle where it meets cylinder to the tip
    # At cylinder junction: z = L/2, r = ±radius
    # The arc angle at junction: sin(angle) = radius/radius_bell
    angle_junction = np.arcsin(radius / radius_bell)
    bell_angles = np.linspace(-angle_junction, angle_junction, 50)
    # Parametric: z = center + R*cos(angle), r = R*sin(angle)
    bell_z_right = bell_center_right - radius_bell * np.cos(bell_angles)
    bell_r_right = radius_bell * np.sin(bell_angles)

    # Create closed polygon for right bell
    right_bell_z = np.concatenate([[length/2], bell_z_right, [length/2]])
    right_bell_r = np.concatenate([[radius], bell_r_right, [-radius]])
    ax_xz.fill(right_bell_z, right_bell_r, 'lightcoral', alpha=0.5)
    ax_xz.plot(bell_z_right, bell_r_right, 'r-', linewidth=2, label='Bell')

    # Left bell (negative z) - sphere center at z = -L/2 - h
    bell_center_left = -length/2 - h
    bell_z_left = bell_center_left + radius_bell * np.cos(bell_angles)
    bell_r_left = radius_bell * np.sin(bell_angles)

    left_bell_z = np.concatenate([[-length/2], bell_z_left, [-length/2]])
    left_bell_r = np.concatenate([[radius], bell_r_left, [-radius]])
    ax_xz.fill(left_bell_z, left_bell_r, 'lightcoral', alpha=0.5)
    ax_xz.plot(bell_z_left, bell_r_left, 'r-', linewidth=2)

    # Add dimension annotations like in documentation
    total_length = length + 2*h + 2*(radius_bell - h)
    ax_xz.annotate('', xy=(length/2, 0), xytext=(-length/2, 0),
                  arrowprops=dict(arrowstyle='<->', color='black'))
    ax_xz.text(0, -radius*1.5, f'L={length:.0f}Å', ha='center', fontsize=9)

    ax_xz.set_xlim(-length/2 - radius_bell*1.3, length/2 + radius_bell*1.3)
    ax_xz.set_ylim(-radius_bell*1.4, radius_bell*1.4)
    ax_xz.set_xlabel('Z (Å) - along axis')
    ax_xz.set_ylabel('Radial (Å)')
    ax_xz.set_title('Side View (like documentation)')
    ax_xz.grid(True, alpha=0.3)
    ax_xz.legend(fontsize=9)

    # YZ plane - same profile
    ax_yz.fill(cyl_z, cyl_r, 'lightgreen', alpha=0.5)
    ax_yz.plot(cyl_z, cyl_r, 'g-', linewidth=2)
    ax_yz.fill(right_bell_z, right_bell_r, 'moccasin', alpha=0.5)
    ax_yz.plot(bell_z_right, bell_r_right, 'orange', linewidth=2)
    ax_yz.fill(left_bell_z, left_bell_r, 'moccasin', alpha=0.5)
    ax_yz.plot(bell_z_left, bell_r_left, 'orange', linewidth=2)

    ax_yz.set_xlim(-length/2 - radius_bell*1.3, length/2 + radius_bell*1.3)
    ax_yz.set_ylim(-radius_bell*1.4, radius_bell*1.4)
    ax_yz.set_xlabel('Z (Å)')
    ax_yz.set_ylabel('Radial (Å)')
    ax_yz.set_title('YZ Cross-section')
    ax_yz.grid(True, alpha=0.3)
valid = "radius_bell >= radius"
have_Fq = True
radius_effective_modes = [
    "equivalent cylinder excluded volume", "equivalent volume sphere",
    "radius", "half length", "half total length",
    ]

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
        radius_bell=bell_radius,
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
    [{}, 0.075, 25.5691260532],
    [{'theta':80., 'phi':10.}, (qx, qy), 3.04233067789],
]
del qx, qy  # not necessary to delete, but cleaner
