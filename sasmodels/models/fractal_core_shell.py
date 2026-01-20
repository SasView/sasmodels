r"""
Definition
----------
Calculates the scattering from a fractal structure with a primary building
block of core-shell spheres, as opposed to just homogeneous spheres in the
fractal model. It is an extension of the well known Teixeira\ [#Teixeira1988]_
fractal model replacing the $P(q)$ of a solid sphere with that of a core-shell
sphere. This model could find use for aggregates of coated particles, or
aggregates of vesicles for example.

.. math::

    I(q) = P(q)S(q) + \text{background}

Where $P(q)$ is the core-shell form factor and $S(q)$ is the
Teixeira\ [#Teixeira1988]_ fractal structure factor both of which are given
again below:

.. math::

    P(q) &= \frac{\phi}{V_s}\left[3V_c(\rho_c-\rho_s)
    \frac{\sin(qr_c)-qr_c\cos(qr_c)}{(qr_c)^3}+
    3V_s(\rho_s-\rho_{solv})
    \frac{\sin(qr_s)-qr_s\cos(qr_s)}{(qr_s)^3}\right]^2 \\
    S(q) &= 1 + \frac{D_f\ \Gamma\!(D_f-1)}{[1+1/(q\xi)^2]^{(D_f-1)/2}}
    \frac{\sin[(D_f-1)\tan^{-1}(q\xi)]}{(qr_s)^{D_f}}

where $\phi$ is the volume fraction of particles, $V_s$ is the volume of the
whole particle, $V_c$ is the volume of the core, $\rho_c$, $\rho_s$, and
$\rho_{solv}$ are the scattering length densities of the core, shell, and
solvent respectively, $r_c$ and $r_s$ are the radius of the core and the
radius of the whole particle respectively, $D_f$ is the fractal dimension,
and $\xi$ the correlation length.

Polydispersity of radius and thickness are also provided for.

This model does not allow for anisotropy and thus the 2D scattering intensity
is calculated in the same way as 1D, where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}

Our model is derived from the form factor calculations implemented in IGOR
macros by the NIST Center for Neutron Research\ [#Kline2006]_

References
----------

.. [#Teixeira1988] J Teixeira, *J. Appl. Cryst.*, 21 (1988) 781-785
.. [#Kline2006]  S R Kline, *J Appl. Cryst.*, 39 (2006) 895

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Paul Butler and Paul Kienzle **Date:** November 27, 2016
* **Last Reviewed by:** Paul Butler and Paul Kienzle **Date:** November 27, 2016
"""

import numpy as np
from numpy import inf

name = "fractal_core_shell"
title = "Scattering from a fractal structure formed from core shell spheres"
description = """\
    Model for fractal aggregates of core-shell primary particles. It is based on
    the Teixeira model for the S(q) of a fractal * P(q) for a core-shell sphere

    radius =  the radius of the core
    thickness = thickness of the shell
    thick_layer = thickness of a layer
    sld_core = the SLD of the core
    sld_shell = the SLD of the shell
    sld_solvent = the SLD of the solvent
    volfraction = volume fraction of core-shell particles
    fractal_dim = fractal dimension
    cor_length = correlation length of the fractal like aggretates
    """
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["radius",      "Ang",        60.0, [0.0, inf],  "volume", "Sphere core radius"],
    ["thickness",   "Ang",        10.0, [0.0, inf],  "volume", "Sphere shell thickness"],
    ["sld_core",    "1e-6/Ang^2", 1.0,  [-inf, inf], "sld",    "Sphere core scattering length density"],
    ["sld_shell",   "1e-6/Ang^2", 2.0,  [-inf, inf], "sld",    "Sphere shell scattering length density"],
    ["sld_solvent", "1e-6/Ang^2", 3.0,  [-inf, inf], "sld",    "Solvent scattering length density"],
    ["volfraction", "",           0.05,  [0.0, inf],  "",       "Volume fraction of building block spheres"],
    ["fractal_dim",    "",        2.0,  [0.0, 6.0],  "",       "Fractal dimension"],
    ["cor_length",  "Ang",      100.0,  [0.0, inf],  "",       "Correlation length of fractal-like aggregates"],
]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sas_3j1x_x.c", "lib/sas_gamma.c", "lib/core_shell.c",
          "lib/fractal_sq.c", "fractal_core_shell.c"]

has_shape_visualization = True

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for fractal core-shell aggregate visualization."""
    import numpy as np
    radius = params.get('radius', 60)
    thickness = params.get('thickness', 10)
    fractal_dim = params.get('fractal_dim', 2)

    outer_radius = radius + thickness

    phi = np.linspace(0, np.pi, resolution//3)
    theta = np.linspace(0, 2*np.pi, resolution//2)
    phi_mesh, theta_mesh = np.meshgrid(phi, theta)

    mesh_data = {}

    # Generate fractal-like distribution
    np.random.seed(42)
    n_spheres = min(20, int(8 * fractal_dim))

    positions = [(0, 0, 0)]
    for i in range(1, n_spheres):
        parent = positions[np.random.randint(len(positions))]
        angle_theta = np.random.uniform(0, 2*np.pi)
        angle_phi = np.random.uniform(0, np.pi)
        dist = 2 * outer_radius * (1 + 0.1 * np.random.randn())

        new_x = parent[0] + dist * np.sin(angle_phi) * np.cos(angle_theta)
        new_y = parent[1] + dist * np.sin(angle_phi) * np.sin(angle_theta)
        new_z = parent[2] + dist * np.cos(angle_phi)
        positions.append((new_x, new_y, new_z))

    for i, (px, py, pz) in enumerate(positions):
        # Outer shell
        x_out = outer_radius * np.sin(phi_mesh) * np.cos(theta_mesh) + px
        y_out = outer_radius * np.sin(phi_mesh) * np.sin(theta_mesh) + py
        z_out = outer_radius * np.cos(phi_mesh) + pz
        mesh_data[f'shell_{i}'] = (x_out, y_out, z_out)

        # Inner core
        x_core = radius * np.sin(phi_mesh) * np.cos(theta_mesh) + px
        y_core = radius * np.sin(phi_mesh) * np.sin(theta_mesh) + py
        z_core = radius * np.cos(phi_mesh) + pz
        mesh_data[f'core_{i}'] = (x_core, y_core, z_core)

    return mesh_data

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of fractal core-shell aggregate."""
    import numpy as np
    radius = params.get('radius', 60)
    thickness = params.get('thickness', 10)
    fractal_dim = params.get('fractal_dim', 2)

    outer_radius = radius + thickness
    theta = np.linspace(0, 2*np.pi, 100)

    # Generate same positions as mesh
    np.random.seed(42)
    n_spheres = min(20, int(8 * fractal_dim))

    positions = [(0, 0, 0)]
    for i in range(1, n_spheres):
        parent = positions[np.random.randint(len(positions))]
        angle_theta = np.random.uniform(0, 2*np.pi)
        angle_phi = np.random.uniform(0, np.pi)
        dist = 2 * outer_radius * (1 + 0.1 * np.random.randn())

        new_x = parent[0] + dist * np.sin(angle_phi) * np.cos(angle_theta)
        new_y = parent[1] + dist * np.sin(angle_phi) * np.sin(angle_theta)
        new_z = parent[2] + dist * np.cos(angle_phi)
        positions.append((new_x, new_y, new_z))

    positions = np.array(positions)
    max_extent = np.max(np.abs(positions)) + outer_radius * 1.5

    # XY plane
    for px, py, pz in positions:
        # Shell
        shell_x = outer_radius * np.cos(theta) + px
        shell_y = outer_radius * np.sin(theta) + py
        ax_xy.plot(shell_x, shell_y, 'b-', linewidth=0.8)
        ax_xy.fill(shell_x, shell_y, 'lightblue', alpha=0.2)
        # Core
        core_x = radius * np.cos(theta) + px
        core_y = radius * np.sin(theta) + py
        ax_xy.fill(core_x, core_y, 'coral', alpha=0.4)

    ax_xy.set_xlim(-max_extent, max_extent)
    ax_xy.set_ylim(-max_extent, max_extent)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title(f'XY Projection (Df={fractal_dim:.1f})')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)

    # XZ plane
    for px, py, pz in positions:
        shell_x = outer_radius * np.cos(theta) + px
        shell_z = outer_radius * np.sin(theta) + pz
        ax_xz.plot(shell_x, shell_z, 'b-', linewidth=0.8)
        ax_xz.fill(shell_x, shell_z, 'lightblue', alpha=0.2)
        core_x = radius * np.cos(theta) + px
        core_z = radius * np.sin(theta) + pz
        ax_xz.fill(core_x, core_z, 'coral', alpha=0.4)

    ax_xz.set_xlim(-max_extent, max_extent)
    ax_xz.set_ylim(-max_extent, max_extent)
    ax_xz.set_xlabel('X (Å)')
    ax_xz.set_ylabel('Z (Å)')
    ax_xz.set_title('XZ Projection')
    ax_xz.set_aspect('equal')
    ax_xz.grid(True, alpha=0.3)

    # YZ plane
    for px, py, pz in positions:
        shell_y = outer_radius * np.cos(theta) + py
        shell_z = outer_radius * np.sin(theta) + pz
        ax_yz.plot(shell_y, shell_z, 'b-', linewidth=0.8)
        ax_yz.fill(shell_y, shell_z, 'lightblue', alpha=0.2)
        core_y = radius * np.cos(theta) + py
        core_z = radius * np.sin(theta) + pz
        ax_yz.fill(core_y, core_z, 'coral', alpha=0.4)

    ax_yz.set_xlim(-max_extent, max_extent)
    ax_yz.set_ylim(-max_extent, max_extent)
    ax_yz.set_xlabel('Y (Å)')
    ax_yz.set_ylabel('Z (Å)')
    ax_yz.set_title('YZ Projection')
    ax_yz.set_aspect('equal')
    ax_yz.grid(True, alpha=0.3)

def random():
    """Return a random parameter set for the model."""
    outer_radius = 10**np.random.uniform(0.7, 4)
    # Use a distribution with a preference for thin shell or thin core
    # Avoid core,shell radii < 1
    thickness = np.random.beta(0.5, 0.5)*(outer_radius-2) + 1
    radius = outer_radius - thickness
    cor_length = 10**np.random.uniform(0.7, 2)*outer_radius
    volfraction = 10**np.random.uniform(-3, -1)
    fractal_dim = 2*np.random.beta(3, 4) + 1
    pars = dict(
        #background=0, sld_block=1, sld_solvent=0,
        volfraction=volfraction,
        radius=radius,
        cor_length=cor_length,
        fractal_dim=fractal_dim,
    )
    return pars

#tests = [[{'radius': 20.0, 'thickness': 10.0}, 'ER', 30.0],
tests = [
    # At some point the SasView 3.x test result was deemed incorrect.  The
    # following tests were verified against NIST IGOR macros ver 7.850.
    # NOTE: NIST macros do only provide for a polydisperse core (no option
    # for a poly shell or for a monodisperse core.  The results seemed
    # extremely sensitive to the core PD, varying non monotonically all
    # the way to a PD of 1e-6. From 1e-6 to 1e-9 no changes in the
    # results were observed and the values below were taken using PD=1e-9.
    # Non-monotonically = I(0.001)=188 to 140 to 177 back to 160 etc.
    [{'radius': 20.0,
      'thickness': 5.0,
      'sld_core': 3.5,
      'sld_shell': 1.0,
      'sld_solvent': 6.35,
      'volfraction': 0.05,
      'background': 0.0},
     [0.001, 0.00291, 0.0107944, 0.029923, 0.100726, 0.476304],
     [177.146, 165.151, 84.1596, 20.1466, 1.40906, 0.00622666]
    ]
]
