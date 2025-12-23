r"""
This model provides the form factor for $N$ spherical pearls of radius $R$
linearly joined by short strings (or segment length or edge separation)
$l$ $(= A - 2R)$. $A$ is the center-to-center pearl separation distance.
The thickness of each string is assumed to be negligible.

.. figure:: img/linear_pearls_geometry.jpg


Definition
----------

The output of the scattering intensity function for the linear_pearls model
is given by (Dobrynin, 1996)

.. math::

    P(Q) = \frac{\text{scale}}{V}\left[ m_{p}^2
    \left(N+2\sum_{n-1}^{N-1}(N-n)\frac{\sin(qnl)}{qnl}\right)
    \left( 3\frac{\sin(qR)-qR\cos(qR)}{(qr)^3}\right)^2\right]

where the mass $m_p$ is $(SLD_{pearl}-SLD_{solvent})*(volume\ of\ N\ pearls)$.
V is the total volume.

The 2D scattering intensity is the same as P(q) above,
regardless of the orientation of the q vector.

References
----------

#.  A V Dobrynin, M Rubinstein and S P Obukhov, *Macromol.*, 29 (1996) 2974-2979

Authorship and Verification
----------------------------

* **Author:**
* **Last Modified by:**
* **Last Reviewed by:**
"""

import numpy as np
from numpy import inf

name = "linear_pearls"
title = "Linear pearls model of scattering from spherical pearls."
description = """
    Calculate form factor for Pearl Necklace Model
    [Macromol. 1996, 29, 2974-2979]
    Parameters:

    sld_pearl: the SLD of the pearl spheres
    sld_solv: the SLD of the solvent
    num_pearls: number of the pearls
    radius: the radius of a pearl
    edge_separation: the length of string segment; surface to surface
    """
category = "shape:sphere"

# pylint: disable=bad-whitespace, line-too-long
#            ["name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["radius",      "Ang",       80.0, [0, inf],     "", "Radius of the pearls"],
    ["edge_sep",    "Ang",      350.0, [0, inf],     "", "Length of the string segment - surface to surface"],
    ["num_pearls",  "",           3.0, [1, inf],     "", "Number of the pearls"],
    ["sld",   "1e-6/Ang^2", 1.0, [-inf, inf],  "sld", "SLD of the pearl spheres"],
    ["sld_solvent", "1e-6/Ang^2", 6.3, [-inf, inf],  "sld", "SLD of the solvent"],
    ]
# pylint: enable=bad-whitespace, line-too-long
has_shape_visualization = True

def create_shape_mesh(params, resolution=40):
    """Create 3D mesh for linear pearls visualization."""
    import numpy as np
    radius = params.get('radius', 80.0)
    edge_sep = params.get('edge_sep', 350.0)
    num_pearls = int(round(params.get('num_pearls', 3)))
    num_pearls = max(num_pearls, 1)

    # Spacing between pearl centers
    center_step = 2 * radius + edge_sep
    z_positions = [
        (i - (num_pearls - 1) / 2.0) * center_step
        for i in range(num_pearls)
    ]

    # Sphere mesh
    phi = np.linspace(0, np.pi, resolution // 2)
    theta = np.linspace(0, 2 * np.pi, resolution)
    phi_mesh, theta_mesh = np.meshgrid(phi, theta)

    pearls = {}
    for i, z0 in enumerate(z_positions):
        x = radius * np.sin(phi_mesh) * np.cos(theta_mesh)
        y = radius * np.sin(phi_mesh) * np.sin(theta_mesh)
        z = radius * np.cos(phi_mesh) + z0
        pearls[f'pearl_{i}'] = (x, y, z)

    return pearls

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of the linear pearls."""
    import numpy as np
    radius = params.get('radius', 80.0)
    edge_sep = params.get('edge_sep', 350.0)
    num_pearls = int(round(params.get('num_pearls', 3)))
    num_pearls = max(num_pearls, 1)

    center_step = 2 * radius + edge_sep
    z_positions = [
        (i - (num_pearls - 1) / 2.0) * center_step
        for i in range(num_pearls)
    ]

    theta = np.linspace(0, 2 * np.pi, 100)
    circle_x = radius * np.cos(theta)
    circle_y = radius * np.sin(theta)

    # XY plane - single pearl cross-section
    ax_xy.plot(circle_x, circle_y, 'b-', linewidth=2)
    ax_xy.fill(circle_x, circle_y, 'lightblue', alpha=0.3)
    ax_xy.set_xlim(-radius*1.5, radius*1.5)
    ax_xy.set_ylim(-radius*1.5, radius*1.5)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title('XY Cross-section (Single Pearl)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)

    # XZ plane - all pearls in a line
    colors = ['lightblue', 'lightcoral', 'lightgreen', 'lightyellow', 'lightpink']
    for i, z0 in enumerate(z_positions):
        color = colors[i % len(colors)]
        pearl_z = circle_y + z0
        ax_xz.plot(circle_x, pearl_z, 'b-', linewidth=2)
        ax_xz.fill(circle_x, pearl_z, color, alpha=0.5)

    # Draw connecting lines
    for i in range(num_pearls - 1):
        ax_xz.plot([0, 0], [z_positions[i] + radius, z_positions[i+1] - radius],
                   'k-', linewidth=1, alpha=0.5)

    total_length = z_positions[-1] - z_positions[0] + 2*radius if num_pearls > 1 else 2*radius
    ax_xz.set_xlim(-radius*2, radius*2)
    ax_xz.set_ylim(z_positions[0] - radius*1.5, z_positions[-1] + radius*1.5)
    ax_xz.set_xlabel('X (Å)')
    ax_xz.set_ylabel('Z (Å)')
    ax_xz.set_title(f'XZ Cross-section ({num_pearls} pearls)')
    ax_xz.grid(True, alpha=0.3)

    # YZ plane - same as XZ
    for i, z0 in enumerate(z_positions):
        color = colors[i % len(colors)]
        pearl_z = circle_y + z0
        ax_yz.plot(circle_x, pearl_z, 'g-', linewidth=2)
        ax_yz.fill(circle_x, pearl_z, color, alpha=0.5)

    for i in range(num_pearls - 1):
        ax_yz.plot([0, 0], [z_positions[i] + radius, z_positions[i+1] - radius],
                   'k-', linewidth=1, alpha=0.5)

    ax_yz.set_xlim(-radius*2, radius*2)
    ax_yz.set_ylim(z_positions[0] - radius*1.5, z_positions[-1] + radius*1.5)
    ax_yz.set_xlabel('Y (Å)')
    ax_yz.set_ylabel('Z (Å)')
    ax_yz.set_title(f'YZ Cross-section ({num_pearls} pearls)')
    ax_yz.grid(True, alpha=0.3)
single = False

source = ["lib/sas_3j1x_x.c", "linear_pearls.c"]

def random():
    """Return a random parameter set for the model."""
    radius = 10**np.random.uniform(1, 3) # 1 - 1000
    edge_sep = 10**np.random.uniform(0, 3)  # 1 - 1000
    num_pearls = np.round(10**np.random.uniform(0.3, 3)) # 2 - 1000
    pars = dict(
        radius=radius,
        edge_sep=edge_sep,
        num_pearls=num_pearls,
    )
    return pars

_ = """
Tests temporarily disabled, until single-double precision accuracy issue solved.

tests = [
    # Accuracy tests based on content in test/utest_model_pearlnecklace.py
    [{'radius':      20.0,
      'num_pearls':   2.0,
      'sld':    1.0,
      'sld_solvent':  6.3,
      'edge_sep':   400.0,
     }, 0.001, 185.135],

    # Additional tests with larger range of parameters
    [{'radius':     120.0,
      'num_pearls':   5.0,
      'sld':    2.0,
      'sld_solvent':  2.3,
      'edge_sep':   100.0,
     }, 0.01, 45.4984],

    [{'radius':       7.0,
      'num_pearls':   2.0,
      'sld':   10.0,
      'sld_solvent': 99.3,
      'edge_sep':    20.0,
     }, 1.0, 0.632811],
    ]
"""
