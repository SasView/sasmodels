r"""
.. _core_shell_sphere:

This model provides the form factor, $P(q)$, for a spherical particle with
a core-shell structure. The form factor is normalized by the particle volume.

For information about polarised and magnetic scattering, see
the :ref:`magnetism` documentation.

Definition
----------

The 1D scattering intensity is calculated in the following way (Guinier, 1955)

.. math::

    P(q) = \frac{\text{scale}}{V} F^2(q) + \text{background}

where

.. math::

    F(q) = 3\left[
       V_c(\rho_c-\rho_s)\frac{\sin(qr_c)-qr_c\cos(qr_c)}{(qr_c)^3} +
       V(\rho_s-\rho_\text{solv})\frac{\sin(qr_s)-qr_s\cos(qr_s)}{(qr_s)^3}
       \right]

$V$ is the volume of the whole particle, $V_c$ is the volume of the
core, $r_s$ = $radius$ + $thickness$ is the radius of the particle, $r_c$
is the radius of the core, $\rho_c$ is the scattering length density of the
core, $\rho_s$ is the scattering length density of the shell,
$\rho_\text{solv}$, is the scattering length density of the solvent.

The 2D scattering intensity is the same as $P(q)$ above, regardless of the
orientation of the $q$ vector.

NB: The outer most radius (ie, = radius + thickness) is used as the
effective radius for $S(Q)$ when $P(Q) \cdot S(Q)$ is applied.

Validation
----------

Validation of our code was done by comparing the output of the 1D model to
the output of the software provided by NIST (Kline, 2006). Figure 1 shows a
comparison of the output of our model and the output of the NIST software.

References
----------

#. A Guinier and G Fournet, *Small-Angle Scattering of X-Rays*,
   John Wiley and Sons, New York, (1955)

Authorship and Verification
----------------------------

* **Author:**
* **Last Modified by:**
* **Last Reviewed by:**
"""

import numpy as np
from numpy import inf, pi

name = "core_shell_sphere"
title = "Form factor for a monodisperse spherical particle with particle with a core-shell structure."
description = """
    F(q) = [V_c (sld_core-sld_shell) 3 (sin(q*radius)-q*radius*cos(q*radius))/(q*radius)^3
            + V_s (sld_shell-sld_solvent) 3 (sin(q*r_s)-q*r_s*cos(q*r_s))/(q*r_s)^3]

            V_s: Volume of the sphere shell
            V_c: Volume of the sphere core
            r_s: Shell radius = radius + thickness
"""
category = "shape:sphere"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["radius",      "Ang",        60.0, [0, inf],    "volume", "Sphere core radius"],
              ["thickness",   "Ang",        10.0, [0, inf],    "volume", "Sphere shell thickness"],
              ["sld_core",    "1e-6/Ang^2", 1.0,  [-inf, inf], "sld",    "core scattering length density"],
              ["sld_shell",   "1e-6/Ang^2", 2.0,  [-inf, inf], "sld",    "shell scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 3.0,  [-inf, inf], "sld",    "Solvent scattering length density"]]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sas_3j1x_x.c", "lib/core_shell.c", "core_shell_sphere.c"]
have_Fq = True
radius_effective_modes = ["outer radius", "core radius"]
has_shape_visualization = True

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for core_shell_sphere visualization."""
    import numpy as np
    radius = params.get('radius', 60)
    thickness = params.get('thickness', 10)
    
    phi = np.linspace(0, np.pi, resolution//2)
    theta = np.linspace(0, 2*np.pi, resolution)
    phi_mesh, theta_mesh = np.meshgrid(phi, theta)
    
    # Core
    x_core = radius * np.sin(phi_mesh) * np.cos(theta_mesh)
    y_core = radius * np.sin(phi_mesh) * np.sin(theta_mesh)
    z_core = radius * np.cos(phi_mesh)
    
    # Shell
    shell_radius = radius + thickness
    x_shell = shell_radius * np.sin(phi_mesh) * np.cos(theta_mesh)
    y_shell = shell_radius * np.sin(phi_mesh) * np.sin(theta_mesh)
    z_shell = shell_radius * np.cos(phi_mesh)
    
    return {
        'core': (x_core, y_core, z_core),
        'shell': (x_shell, y_shell, z_shell)
    }

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of the core_shell_sphere."""
    import numpy as np
    radius = params.get('radius', 60)
    thickness = params.get('thickness', 10)
    shell_radius = radius + thickness
    
    # Create circles for core and shell
    theta = np.linspace(0, 2*np.pi, 100)
    
    # Core circles
    core_x = radius * np.cos(theta)
    core_y = radius * np.sin(theta)
    
    # Shell circles
    shell_x = shell_radius * np.cos(theta)
    shell_y = shell_radius * np.sin(theta)
    
    # XY plane (top view)
    ax_xy.plot(shell_x, shell_y, 'r-', linewidth=2, label='Shell')
    ax_xy.fill(shell_x, shell_y, 'lightcoral', alpha=0.3)
    ax_xy.plot(core_x, core_y, 'b-', linewidth=2, label='Core')
    ax_xy.fill(core_x, core_y, 'lightblue', alpha=0.5)
    
    ax_xy.set_xlim(-shell_radius*1.2, shell_radius*1.2)
    ax_xy.set_ylim(-shell_radius*1.2, shell_radius*1.2)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title('XY Cross-section (Top View)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)
    ax_xy.legend()
    
    # XZ plane (side view)
    ax_xz.plot(shell_x, shell_y, 'r-', linewidth=2, label='Shell')
    ax_xz.fill(shell_x, shell_y, 'lightcoral', alpha=0.3)
    ax_xz.plot(core_x, core_y, 'b-', linewidth=2, label='Core')
    ax_xz.fill(core_x, core_y, 'lightblue', alpha=0.5)
    
    ax_xz.set_xlim(-shell_radius*1.2, shell_radius*1.2)
    ax_xz.set_ylim(-shell_radius*1.2, shell_radius*1.2)
    ax_xz.set_xlabel('X (Å)')
    ax_xz.set_ylabel('Z (Å)')
    ax_xz.set_title('XZ Cross-section (Side View)')
    ax_xz.set_aspect('equal')
    ax_xz.grid(True, alpha=0.3)
    ax_xz.legend()
    
    # YZ plane (front view)
    ax_yz.plot(shell_x, shell_y, 'g-', linewidth=2, label='Shell')
    ax_yz.fill(shell_x, shell_y, 'lightgreen', alpha=0.3)
    ax_yz.plot(core_x, core_y, 'orange', linewidth=2, label='Core')
    ax_yz.fill(core_x, core_y, 'moccasin', alpha=0.5)
    
    ax_yz.set_xlim(-shell_radius*1.2, shell_radius*1.2)
    ax_yz.set_ylim(-shell_radius*1.2, shell_radius*1.2)
    ax_yz.set_xlabel('Y (Å)')
    ax_yz.set_ylabel('Z (Å)')
    ax_yz.set_title('YZ Cross-section (Front View)')
    ax_yz.set_aspect('equal')
    ax_yz.grid(True, alpha=0.3)
    ax_yz.legend()
    
    # Add dimension annotations
    ax_xz.annotate('', xy=(-radius, 0), xytext=(radius, 0),
                  arrowprops=dict(arrowstyle='<->', color='blue'))
    ax_xz.text(0, -radius*0.3, f'Core R = {radius:.0f} Å', ha='center', fontsize=10, color='blue')
    
    ax_xz.annotate('', xy=(-shell_radius, -radius*0.7), xytext=(shell_radius, -radius*0.7),
                  arrowprops=dict(arrowstyle='<->', color='red'))
    ax_xz.text(0, -radius*0.9, f'Shell R = {shell_radius:.0f} Å', ha='center', fontsize=10, color='red')
    
    ax_xz.annotate('', xy=(radius*0.7, 0), xytext=(shell_radius*0.7, 0),
                  arrowprops=dict(arrowstyle='<->', color='black'))
    ax_xz.text(radius*0.85, radius*0.2, f't = {thickness:.0f} Å', ha='center', fontsize=10, rotation=90)

def random():
    """Return a random parameter set for the model."""
    outer_radius = 10**np.random.uniform(1.3, 4.3)
    # Use a distribution with a preference for thin shell or thin core
    # Avoid core,shell radii < 1
    radius = np.random.beta(0.5, 0.5)*(outer_radius-2) + 1
    thickness = outer_radius - radius
    pars = dict(
        radius=radius,
        thickness=thickness,
    )
    return pars

tests = [
    [{'radius': 20.0, 'thickness': 10.0}, 0.1, None, None, 30.0, 4.*pi/3*30**3, 1.0],

    # The SasView test result was 0.00169, with a background of 0.001
    [{'radius': 60.0, 'thickness': 10.0, 'sld_core': 1.0, 'sld_shell': 2.0,
      'sld_solvent': 3.0, 'background': 0.0}, 0.4, 0.000698838],
]
