r"""
Definition
----------
This model calculates the scattering from fractal-like aggregates of spherical
building blocks according the following equation:

.. math::

    I(q) = \phi\ V_\text{block} (\rho_\text{block}
          - \rho_\text{solvent})^2 P(q)S(q) + \text{background}

where $\phi$ is The volume fraction of the spherical "building block" particles
of radius $R_0$, $V_{block}$ is the volume of a single building block,
$\rho_{solvent}$ is the scattering length density of the solvent, and
$\rho_{block}$ is the scattering length density of the building blocks, and
P(q), S(q) are the scattering from randomly distributed spherical particles
(the building blocks) and the interference from such building blocks organized
in a fractal-like clusters.  P(q) and S(q) are calculated as:

.. math::

    P(q)&= F(qR_0)^2 \\
    F(q)&= \frac{3 (\sin x - x \cos x)}{x^3} \\
    V_\text{particle} &= \frac{4}{3}\ \pi R_0 \\
    S(q) &= 1 + \frac{D_f\  \Gamma\!(D_f-1)}{[1+1/(q \xi)^2\  ]^{(D_f -1)/2}}
    \frac{\sin[(D_f-1) \tan^{-1}(q \xi) ]}{(q R_0)^{D_f}}

where $\xi$ is the correlation length representing the cluster size and $D_f$
is the fractal dimension, representing the self similarity of the structure.
Note that S(q) here goes negative if $D_f$ is too large, and the Gamma function
diverges at $D_f=0$ and $D_f=1$.

**Polydispersity on the radius is provided for.**

For 2D data: The 2D scattering intensity is calculated in the same way as
1D, where the *q* vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


References
----------

#.  J Teixeira, *J. Appl. Cryst.*, 21 (1988) 781-785

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Converted to sasmodels by:** Paul Butler **Date:** March 19, 2016
* **Last Modified by:** Paul Butler **Date:** March 12, 2017
* **Last Reviewed by:** Paul Butler **Date:** March 12, 2017
"""

import numpy as np
from numpy import inf

name = "fractal"
title = "Calculates the scattering from fractal-like aggregates of spheres \
following theTexiera reference."
description = """
        The scattering intensity is given by
        I(q) = scale * V * delta^(2) * P(q) * S(q) + background, where
        p(q)= F(q*radius)^(2)
        F(x) = 3*[sin(x)-x cos(x)]/x**3
        delta = sld_block -sld_solv
        scale        =  scale * volfraction
        radius       =  Block radius
        sld_block    =  SDL block
        sld_solv  =  SDL solvent
        background   =  background
        and S(q) is the interference term between building blocks given
        in the full documentation and depending on the parameters
        fractal_dim  =  Fractal dimension
        cor_length  =  Correlation Length    """

category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["volfraction", "", 0.05, [0.0, 1], "",
               "volume fraction of blocks"],
              ["radius",    "Ang",  5.0, [0.0, inf], "volume",
               "radius of particles"],
              ["fractal_dim",      "",  2.0, [0.0, 6.0], "",
               "fractal dimension"],
              ["cor_length", "Ang", 100.0, [0.0, inf], "",
               "cluster correlation length"],
              ["sld_block", "1e-6/Ang^2", 2.0, [-inf, inf], "sld",
               "scattering length density of particles"],
              ["sld_solvent", "1e-6/Ang^2", 6.4, [-inf, inf], "sld",
               "scattering length density of solvent"],
             ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sas_3j1x_x.c", "lib/sas_gamma.c", "lib/fractal_sq.c", "fractal.c"]
has_shape_visualization = True
valid = "fractal_dim >= 0.0"

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for fractal aggregate visualization."""
    import numpy as np
    radius = params.get('radius', 5)
    cor_length = params.get('cor_length', 100)
    fractal_dim = params.get('fractal_dim', 2)
    
    phi = np.linspace(0, np.pi, resolution//3)
    theta = np.linspace(0, 2*np.pi, resolution//2)
    phi_mesh, theta_mesh = np.meshgrid(phi, theta)
    
    mesh_data = {}
    
    # Generate fractal-like distribution of spheres
    # Use a simplified DLA-like arrangement
    np.random.seed(42)  # For reproducibility
    n_spheres = min(30, int(10 * fractal_dim))
    
    # Central sphere
    x = radius * np.sin(phi_mesh) * np.cos(theta_mesh)
    y = radius * np.sin(phi_mesh) * np.sin(theta_mesh)
    z = radius * np.cos(phi_mesh)
    mesh_data['sphere_0'] = (x, y, z)
    
    # Add surrounding spheres in a fractal-like pattern
    positions = [(0, 0, 0)]
    for i in range(1, n_spheres):
        # Random walk from existing positions
        parent = positions[np.random.randint(len(positions))]
        angle_theta = np.random.uniform(0, 2*np.pi)
        angle_phi = np.random.uniform(0, np.pi)
        dist = 2 * radius * (1 + 0.1 * np.random.randn())
        
        new_x = parent[0] + dist * np.sin(angle_phi) * np.cos(angle_theta)
        new_y = parent[1] + dist * np.sin(angle_phi) * np.sin(angle_theta)
        new_z = parent[2] + dist * np.cos(angle_phi)
        
        positions.append((new_x, new_y, new_z))
        
        xs = radius * np.sin(phi_mesh) * np.cos(theta_mesh) + new_x
        ys = radius * np.sin(phi_mesh) * np.sin(theta_mesh) + new_y
        zs = radius * np.cos(phi_mesh) + new_z
        mesh_data[f'sphere_{i}'] = (xs, ys, zs)
    
    return mesh_data

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of fractal aggregate."""
    import numpy as np
    radius = params.get('radius', 5)
    fractal_dim = params.get('fractal_dim', 2)
    
    theta = np.linspace(0, 2*np.pi, 100)
    
    # Generate same positions as mesh
    np.random.seed(42)
    n_spheres = min(30, int(10 * fractal_dim))
    
    positions = [(0, 0, 0)]
    for i in range(1, n_spheres):
        parent = positions[np.random.randint(len(positions))]
        angle_theta = np.random.uniform(0, 2*np.pi)
        angle_phi = np.random.uniform(0, np.pi)
        dist = 2 * radius * (1 + 0.1 * np.random.randn())
        
        new_x = parent[0] + dist * np.sin(angle_phi) * np.cos(angle_theta)
        new_y = parent[1] + dist * np.sin(angle_phi) * np.sin(angle_theta)
        new_z = parent[2] + dist * np.cos(angle_phi)
        positions.append((new_x, new_y, new_z))
    
    positions = np.array(positions)
    max_extent = np.max(np.abs(positions)) + radius * 2
    
    # XY plane
    for px, py, pz in positions:
        circle_x = radius * np.cos(theta) + px
        circle_y = radius * np.sin(theta) + py
        ax_xy.plot(circle_x, circle_y, 'b-', linewidth=0.8)
        ax_xy.fill(circle_x, circle_y, 'lightblue', alpha=0.3)
    
    ax_xy.set_xlim(-max_extent, max_extent)
    ax_xy.set_ylim(-max_extent, max_extent)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title(f'XY Projection (Df={fractal_dim:.1f})')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)
    
    # XZ plane
    for px, py, pz in positions:
        circle_x = radius * np.cos(theta) + px
        circle_z = radius * np.sin(theta) + pz
        ax_xz.plot(circle_x, circle_z, 'r-', linewidth=0.8)
        ax_xz.fill(circle_x, circle_z, 'lightcoral', alpha=0.3)
    
    ax_xz.set_xlim(-max_extent, max_extent)
    ax_xz.set_ylim(-max_extent, max_extent)
    ax_xz.set_xlabel('X (Å)')
    ax_xz.set_ylabel('Z (Å)')
    ax_xz.set_title('XZ Projection')
    ax_xz.set_aspect('equal')
    ax_xz.grid(True, alpha=0.3)
    
    # YZ plane
    for px, py, pz in positions:
        circle_y = radius * np.cos(theta) + py
        circle_z = radius * np.sin(theta) + pz
        ax_yz.plot(circle_y, circle_z, 'g-', linewidth=0.8)
        ax_yz.fill(circle_y, circle_z, 'lightgreen', alpha=0.3)
    
    ax_yz.set_xlim(-max_extent, max_extent)
    ax_yz.set_ylim(-max_extent, max_extent)
    ax_yz.set_xlabel('Y (Å)')
    ax_yz.set_ylabel('Z (Å)')
    ax_yz.set_title('YZ Projection')
    ax_yz.set_aspect('equal')
    ax_yz.grid(True, alpha=0.3)

def random():
    """Return a random parameter set for the model."""
    radius = 10**np.random.uniform(0.7, 4)
    #radius = 5
    cor_length = 10**np.random.uniform(0.7, 2)*radius
    #cor_length = 20*radius
    volfraction = 10**np.random.uniform(-3, -1)
    #volfraction = 0.05
    fractal_dim = 2*np.random.beta(3, 4) + 1
    #fractal_dim = 2
    pars = dict(
        #background=0, sld_block=1, sld_solvent=0,
        volfraction=volfraction,
        radius=radius,
        cor_length=cor_length,
        fractal_dim=fractal_dim,
    )
    return pars

# NOTE: test results taken from values returned by SasView 3.1.2
tests = [
    [{}, 0.0005, 40.4980069872],
    [{}, 0.234734468938, 0.0947143166058],
    [{}, 0.5, 0.0176878183458],
    ]
