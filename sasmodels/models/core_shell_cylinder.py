r"""
Definition
----------

The output of the 2D scattering intensity function for oriented core-shell
cylinders is given by Kline [#Kline2006]_. The form factor is normalized
by the particle volume. Note that in this model the shell envelops the entire
core so that besides a "sleeve" around the core, the shell also provides two
flat end caps of thickness = shell thickness. In other words the length of the
total cylinder is the length of the core cylinder plus twice the thickness of
the shell. If no end caps are desired one should use the
:ref:`core-shell-bicelle` and set the thickness of the end caps (in this case
the "thick_face") to zero.

.. math::

    I(q,\alpha) = \frac{\text{scale}}{V_s} F^2(q,\alpha).sin(\alpha) + \text{background}

where

.. math::

    F(q,\alpha) = &\ (\rho_c - \rho_s) V_c
           \frac{\sin \left( q \tfrac12 L\cos\alpha \right)}
                {q \tfrac12 L\cos\alpha}
           \frac{2 J_1 \left( qR\sin\alpha \right)}
                {qR\sin\alpha} \\
         &\ + (\rho_s - \rho_\text{solv}) V_s
           \frac{\sin \left( q \left(\tfrac12 L+T\right) \cos\alpha \right)}
                {q \left(\tfrac12 L +T \right) \cos\alpha}
           \frac{ 2 J_1 \left( q(R+T)\sin\alpha \right)}
                {q(R+T)\sin\alpha}

and

.. math::

    V_s = \pi (R + T)^2 (L + 2T)

and $\alpha$ is the angle between the axis of the cylinder and $\vec q$,
$V_s$ is the total volume (i.e. including both the core and the outer shell),
$V_c$ is the volume of the core, $L$ is the length of the core,
$R$ is the radius of the core, $T$ is the thickness of the shell, $\rho_c$
is the scattering length density of the core, $\rho_s$ is the scattering
length density of the shell, $\rho_\text{solv}$ is the scattering length
density of the solvent, and *background* is the background level.  The outer
radius of the shell is given by $R+T$ and the total length of the outer
shell is given by $L+2T$. $J_1$ is the first order Bessel function.

.. _core-shell-cylinder-geometry:

.. figure:: img/core_shell_cylinder_geometry.jpg

    Core shell cylinder schematic.

To provide easy access to the orientation of the core-shell cylinder, we
define the axis of the cylinder using two angles $\theta$ and $\phi$.
(see :ref:`cylinder model <cylinder-angle-definition>`)

NB: The 2nd virial coefficient of the cylinder is calculated based on
the radius and 2 length values, and used as the effective radius for
$S(q)$ when $P(q) \cdot S(q)$ is applied.

The $\theta$ and $\phi$ parameters are not used for the 1D output.

Reference
---------

See also Livsey [#Livsey1987]_ and Onsager [#Onsager1949]_.

.. [#Livsey1987] I Livsey, *J. Chem. Soc., Faraday Trans. 2*, 83 (1987) 1445-1452

.. [#Kline2006] S R Kline, *J Appl. Cryst.*, 39 (2006) 895

.. [#Onsager1949] L. Onsager, *Ann. New York Acad. Sci.*, 51 (1949) 627-659

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Paul Kienzle **Date:** Aug 8, 2016
* **Last Reviewed by:** Richard Heenan **Date:** March 18, 2016
"""

import numpy as np
from numpy import cos, inf, pi, sin

name = "core_shell_cylinder"
title = "Right circular cylinder with a core-shell scattering length density profile."
description = """
P(q,alpha)= scale/Vs*f(q)^(2) + background,
      where: f(q)= 2(sld_core - solvant_sld)
        * Vc*sin[qLcos(alpha/2)]
        /[qLcos(alpha/2)]*J1(qRsin(alpha))
        /[qRsin(alpha)]+2(sld_shell-sld_solvent)
        *Vs*sin[q(L+T)cos(alpha/2)][[q(L+T)
        *cos(alpha/2)]*J1(q(R+T)sin(alpha))
        /q(R+T)sin(alpha)]

    alpha:is the angle between the axis of
        the cylinder and the q-vector
    Vs: the volume of the outer shell
    Vc: the volume of the core
    L: the length of the core
        sld_shell: the scattering length density of the shell
    sld_solvent: the scattering length density of the solvent
    background: the background
    T: the thickness
        R+T: is the outer radius
     L+2T: The total length of the outershell
    J1: the first order Bessel function
     theta: axis_theta of the cylinder
     phi: the axis_phi of the cylinder
"""
category = "shape:cylinder"

#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld_core", "1e-6/Ang^2", 4, [-inf, inf], "sld",
               "Cylinder core scattering length density"],
              ["sld_shell", "1e-6/Ang^2", 4, [-inf, inf], "sld",
               "Cylinder shell scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["radius", "Ang", 20, [0, inf], "volume",
               "Cylinder core radius"],
              ["thickness", "Ang", 20, [0, inf], "volume",
               "Cylinder shell thickness"],
              ["length", "Ang", 400, [0, inf], "volume",
               "Cylinder length"],
              ["theta", "degrees", 60, [-360, 360], "orientation",
               "cylinder axis to beam angle"],
              ["phi", "degrees", 60, [-360, 360], "orientation",
               "rotation about beam"],
             ]

source = ["lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c", "core_shell_cylinder.c"]
have_Fq = True
radius_effective_modes = [
    "excluded volume", "equivalent volume sphere", "outer radius", "half outer length",
    "half min outer dimension", "half max outer dimension", "half outer diagonal",
    ]
has_shape_visualization = True

def create_shape_mesh(params, resolution=50):
    import numpy as np
    radius = params.get('radius', 20)
    thickness = params.get('thickness', 20)
    length = params.get('length', 400)

    # Outer dimensions
    outer_radius = radius + thickness
    outer_length = length + 2 * thickness

    # Create core cylinder
    theta = np.linspace(0, 2*np.pi, resolution)
    z_core = np.linspace(-length/2, length/2, resolution//2)
    theta_core, z_core_mesh = np.meshgrid(theta, z_core)
    x_core = radius * np.cos(theta_core)
    y_core = radius * np.sin(theta_core)

    # Create shell cylinder (outer surface)
    z_shell = np.linspace(-outer_length/2, outer_length/2, resolution//2)
    theta_shell, z_shell_mesh = np.meshgrid(theta, z_shell)
    x_shell = outer_radius * np.cos(theta_shell)
    y_shell = outer_radius * np.sin(theta_shell)

    # Create end caps (shell only, as annular disks)
    r_cap_inner = np.linspace(0, radius, resolution//8)
    r_cap_outer = np.linspace(radius, outer_radius, resolution//8)
    theta_cap = np.linspace(0, 2*np.pi, resolution)

    # Inner caps (on core)
    r_inner_mesh, theta_inner_mesh = np.meshgrid(r_cap_inner, theta_cap)
    x_cap_core = r_inner_mesh * np.cos(theta_inner_mesh)
    y_cap_core = r_inner_mesh * np.sin(theta_inner_mesh)
    z_cap_core_top = np.full_like(x_cap_core, length/2)
    z_cap_core_bottom = np.full_like(x_cap_core, -length/2)

    # Outer shell caps (annular rings on ends)
    r_outer_mesh, theta_outer_mesh = np.meshgrid(r_cap_outer, theta_cap)
    x_cap_shell = r_outer_mesh * np.cos(theta_outer_mesh)
    y_cap_shell = r_outer_mesh * np.sin(theta_outer_mesh)
    z_cap_shell_top = np.full_like(x_cap_shell, outer_length/2)
    z_cap_shell_bottom = np.full_like(x_cap_shell, -outer_length/2)

    # Middle shell caps (between core and outer shell)
    r_full = np.linspace(0, outer_radius, resolution//4)
    r_full_mesh, theta_full_mesh = np.meshgrid(r_full, theta_cap)
    x_cap_middle = r_full_mesh * np.cos(theta_full_mesh)
    y_cap_middle = r_full_mesh * np.sin(theta_full_mesh)
    z_cap_middle_top = np.full_like(x_cap_middle, length/2)
    z_cap_middle_bottom = np.full_like(x_cap_middle, -length/2)

    return {
        'core_cylinder': (x_core, y_core, z_core_mesh),
        'shell_cylinder': (x_shell, y_shell, z_shell_mesh),
        'shell_cap_top': (x_cap_middle, y_cap_middle, z_cap_middle_top),
        'shell_cap_bottom': (x_cap_middle, y_cap_middle, z_cap_middle_bottom),
        'end_cap_top': (x_cap_shell, y_cap_shell, z_cap_shell_top),
        'end_cap_bottom': (x_cap_shell, y_cap_shell, z_cap_shell_bottom),
    }

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    import numpy as np
    radius = params.get('radius', 20)
    thickness = params.get('thickness', 20)
    length = params.get('length', 400)

    outer_radius = radius + thickness
    outer_length = length + 2 * thickness

    # XY plane (top view) - concentric circles
    theta = np.linspace(0, 2*np.pi, 100)

    # Core circle
    core_x = radius * np.cos(theta)
    core_y = radius * np.sin(theta)

    # Shell circle
    shell_x = outer_radius * np.cos(theta)
    shell_y = outer_radius * np.sin(theta)

    ax_xy.plot(shell_x, shell_y, 'r-', linewidth=2, label='Shell')
    ax_xy.fill(shell_x, shell_y, 'lightcoral', alpha=0.3)
    ax_xy.plot(core_x, core_y, 'b-', linewidth=2, label='Core')
    ax_xy.fill(core_x, core_y, 'lightblue', alpha=0.5)

    ax_xy.set_xlim(-outer_radius*1.3, outer_radius*1.3)
    ax_xy.set_ylim(-outer_radius*1.3, outer_radius*1.3)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title('XY Cross-section (Top View)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)
    ax_xy.legend()

    # XZ plane (side view) - nested rectangles
    # Core rectangle
    core_rect_x = [-length/2, -length/2, length/2, length/2, -length/2]
    core_rect_z = [-radius, radius, radius, -radius, -radius]

    # Shell rectangle
    shell_rect_x = [-outer_length/2, -outer_length/2, outer_length/2, outer_length/2, -outer_length/2]
    shell_rect_z = [-outer_radius, outer_radius, outer_radius, -outer_radius, -outer_radius]

    ax_xz.plot(shell_rect_x, shell_rect_z, 'r-', linewidth=2, label='Shell')
    ax_xz.fill(shell_rect_x, shell_rect_z, 'lightcoral', alpha=0.3)
    ax_xz.plot(core_rect_x, core_rect_z, 'b-', linewidth=2, label='Core')
    ax_xz.fill(core_rect_x, core_rect_z, 'lightblue', alpha=0.5)

    ax_xz.set_xlim(-outer_length/2*1.2, outer_length/2*1.2)
    ax_xz.set_ylim(-outer_radius*1.3, outer_radius*1.3)
    ax_xz.set_xlabel('Z (Å)')
    ax_xz.set_ylabel('X (Å)')
    ax_xz.set_title('XZ Cross-section (Side View)')
    ax_xz.grid(True, alpha=0.3)
    ax_xz.legend()

    # Add dimension annotations
    ax_xz.annotate('', xy=(-length/2, -outer_radius*1.5), xytext=(length/2, -outer_radius*1.5),
                  arrowprops=dict(arrowstyle='<->', color='blue'))
    ax_xz.text(0, -outer_radius*1.6, f'L = {length:.0f} Å (core)', ha='center', fontsize=9)

    ax_xz.annotate('', xy=(-outer_length/2, -outer_radius*1.8), xytext=(outer_length/2, -outer_radius*1.8),
                  arrowprops=dict(arrowstyle='<->', color='red'))
    ax_xz.text(0, -outer_radius*1.9, f'L+2T = {outer_length:.0f} Å (total)', ha='center', fontsize=9, color='red')

    ax_xz.annotate('', xy=(outer_length/2*0.7, 0), xytext=(outer_length/2*0.7, radius),
                  arrowprops=dict(arrowstyle='<->', color='blue'))
    ax_xz.text(outer_length/2*0.75, radius/2, f'r={radius:.0f}', fontsize=9, color='blue')

    ax_xz.annotate('', xy=(outer_length/2*0.85, 0), xytext=(outer_length/2*0.85, outer_radius),
                  arrowprops=dict(arrowstyle='<->', color='red'))
    ax_xz.text(outer_length/2*0.9, outer_radius/2, f'R+T={outer_radius:.0f}', fontsize=9, color='red')

    # YZ plane (front view) - same as XZ
    ax_yz.plot(shell_rect_x, shell_rect_z, 'g-', linewidth=2, label='Shell')
    ax_yz.fill(shell_rect_x, shell_rect_z, 'lightgreen', alpha=0.3)
    ax_yz.plot(core_rect_x, core_rect_z, 'orange', linewidth=2, label='Core')
    ax_yz.fill(core_rect_x, core_rect_z, 'moccasin', alpha=0.5)

    ax_yz.set_xlim(-outer_length/2*1.2, outer_length/2*1.2)
    ax_yz.set_ylim(-outer_radius*1.3, outer_radius*1.3)
    ax_yz.set_xlabel('Z (Å)')
    ax_yz.set_ylabel('Y (Å)')
    ax_yz.set_title('YZ Cross-section (Front View)')
    ax_yz.grid(True, alpha=0.3)
    ax_yz.legend()

def random():
    """Return a random parameter set for the model."""
    outer_radius = 10**np.random.uniform(1, 4.7)
    # Use a distribution with a preference for thin shell or thin core
    # Avoid core,shell radii < 1
    radius = np.random.beta(0.5, 0.5)*(outer_radius-2) + 1
    thickness = outer_radius - radius
    length = np.random.uniform(1, 4.7)
    pars = dict(
        radius=radius,
        thickness=thickness,
        length=length,
    )
    return pars

q = 0.1
# april 6 2017, rkh add unit tests, NOT compared with any other calc method, assume correct!
qx = q*cos(pi/6.0)
qy = q*sin(pi/6.0)
tests = [
    [{}, 0.075, 10.8552692237],
    [{}, (qx, qy), 0.444618752741],
]
del qx, qy  # not necessary to delete, but cleaner
