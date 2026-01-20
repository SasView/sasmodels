r"""
Calculates the scattering from fractal-like aggregates based on
the Mildner reference.

Definition
----------

The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = scale \times P(q)S(q) + background

.. math::

    P(q) = F(qR)^2

.. math::

    F(x) = \frac{3\left[sin(x)-xcos(x)\right]}{x^3}

.. math::

    S(q) = \frac{\Gamma(D_m-1)\zeta^{D_m-1}}{\left[1+(q\zeta)^2
    \right]^{(D_m-1)/2}}
    \frac{sin\left[(D_m - 1) tan^{-1}(q\zeta) \right]}{q}

.. math::

    scale = scale\_factor \times NV^2(\rho_\text{particle} - \rho_\text{solvent})^2

.. math::

    V = \frac{4}{3}\pi R^3

where $R$ is the radius of the building block, $D_m$ is the **mass** fractal
dimension, $\zeta$ is the cut-off length, $\rho_\text{solvent}$ is the
scattering length density of the solvent, and $\rho_\text{particle}$ is the
scattering length density of particles.

.. note::

    The mass fractal dimension ( $D_m$ ) is only
    valid if $1 < mass\_dim < 6$. It is also only valid over a limited
    $q$ range (see the reference for details).


References
----------

#. D Mildner and P Hall,
   *J. Phys. D: Appl. Phys.*, 19 (1986) 1535-1545 Equation(9)

Authorship and Verification
----------------------------

* **Author:**
* **Last Modified by:**
* **Last Reviewed by:**
"""

import numpy as np
from numpy import inf

name = "mass_fractal"
title = "Mass Fractal model"
description = """
        The scattering intensity  I(x) = scale*P(x)*S(x) + background, where
        scale = scale_factor  * V * delta^(2)
        p(x)=  F(x*radius)^(2)
        F(x) = 3*[sin(x)-x cos(x)]/x**3
        S(x) = [(gamma(Dm-1)*colength^(Dm-1)*[1+(x^2*colength^2)]^((1-Dm)/2)
        * sin[(Dm-1)*arctan(x*colength)])/x]
        where delta = sldParticle -sldSolv.
        radius       =  Particle radius
        fractal_dim_mass  =  Mass fractal dimension
        cutoff_length  =  Cut-off length
        background   =  background
        Ref.:Mildner, Hall,J Phys D Appl Phys(1986), 9, 1535-1545
        Note I: This model is valid for 1<fractal_dim_mass<6.
        Note II: This model is not in absolute scale.
        """
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["radius",           "Ang",  10.0, [0.0, inf], "", "Particle radius"],
    ["fractal_dim_mass", "",      1.9, [1.0, 6.0], "", "Mass fractal dimension"],
    ["cutoff_length",    "Ang", 100.0, [0.0, inf], "", "Cut-off length"],
]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sas_3j1x_x.c", "lib/sas_gamma.c", "mass_fractal.c"]
has_shape_visualization = True
valid = "fractal_dim_mass >= 1.0"

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for mass fractal aggregate visualization."""
    import numpy as np
    radius = params.get('radius', 10)
    fractal_dim_mass = params.get('fractal_dim_mass', 1.9)
    
    phi = np.linspace(0, np.pi, resolution//3)
    theta = np.linspace(0, 2*np.pi, resolution//2)
    phi_mesh, theta_mesh = np.meshgrid(phi, theta)
    
    mesh_data = {}
    
    np.random.seed(42)
    n_spheres = min(25, int(8 * fractal_dim_mass))
    
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
        
        xs = radius * np.sin(phi_mesh) * np.cos(theta_mesh) + new_x
        ys = radius * np.sin(phi_mesh) * np.sin(theta_mesh) + new_y
        zs = radius * np.cos(phi_mesh) + new_z
        mesh_data[f'sphere_{i}'] = (xs, ys, zs)
    
    # Central sphere
    x = radius * np.sin(phi_mesh) * np.cos(theta_mesh)
    y = radius * np.sin(phi_mesh) * np.sin(theta_mesh)
    z = radius * np.cos(phi_mesh)
    mesh_data['sphere_0'] = (x, y, z)
    
    return mesh_data

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of mass fractal aggregate."""
    import numpy as np
    radius = params.get('radius', 10)
    fractal_dim_mass = params.get('fractal_dim_mass', 1.9)
    
    theta = np.linspace(0, 2*np.pi, 100)
    
    np.random.seed(42)
    n_spheres = min(25, int(8 * fractal_dim_mass))
    
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
    
    for px, py, pz in positions:
        circle_x = radius * np.cos(theta) + px
        circle_y = radius * np.sin(theta) + py
        ax_xy.plot(circle_x, circle_y, 'b-', linewidth=0.8)
        ax_xy.fill(circle_x, circle_y, 'lightblue', alpha=0.3)
    
    ax_xy.set_xlim(-max_extent, max_extent)
    ax_xy.set_ylim(-max_extent, max_extent)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title(f'XY Projection (Dm={fractal_dim_mass:.1f})')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)
    
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
    cutoff_length = 10**np.random.uniform(0.7, 2)*radius
    # TODO: fractal dimension should range from 1 to 6
    fractal_dim_mass = 2*np.random.beta(3, 4) + 1
    #volfrac = 10**np.random.uniform(-4, -1)
    pars = dict(
        #background=0,
        scale=1, #1e5*volfrac/radius**(fractal_dim_mass),
        radius=radius,
        cutoff_length=cutoff_length,
        fractal_dim_mass=fractal_dim_mass,
    )
    return pars

tests = [

    # Accuracy tests based on content in test/utest_other_models.py
    [{'radius':         10.0,
      'fractal_dim_mass':        1.9,
      'cutoff_length': 100.0,
     }, 0.05, 279.59422],

    # Additional tests with larger range of parameters
    [{'radius':        2.0,
      'fractal_dim_mass':      3.3,
      'cutoff_length': 1.0,
     }, 0.5, 1.29116774904],

    [{'radius':        1.0,
      'fractal_dim_mass':      1.3,
      'cutoff_length': 1.0,
      'background':    0.8,
     }, 0.001, 1.69747015932],

    [{'radius':        1.0,
      'fractal_dim_mass':      2.3,
      'cutoff_length': 1.0,
      'scale':        10.0,
     }, 0.051, 11.6237966145],
    ]
