r"""
Definition
----------

The form factor for this bent disc is essentially that of a hyperbolic
paraboloid and calculated as

.. math::

    P(q) = (\Delta \rho )^2 V \int^{\pi/2}_0 d\psi \sin{\psi} sinc^2
    \left( \frac{qd\cos{\psi}}{2} \right)
    \left[ \left( S^2_0+C^2_0\right) + 2\sum_{n=1}^{\infty}
     \left( S^2_n+C^2_n\right) \right]

where

.. math::

    C_n = \frac{1}{r^2}\int^{R}_{0} r dr\cos(qr^2\alpha \cos{\psi})
    J_n\left( qr^2\beta \cos{\psi}\right)
    J_{2n}\left( qr \sin{\psi}\right)

.. math::

    S_n = \frac{1}{r^2}\int^{R}_{0} r dr\sin(qr^2\alpha \cos{\psi})
    J_n\left( qr^2\beta \cos{\psi}\right)
    J_{2n}\left( qr \sin{\psi}\right)

and $\Delta\rho\text{ is }\rho_{pringle}-\rho_{solvent}$, $V$ is the volume of
the disc, $\psi$ is the angle between the normal to the disc and the q vector,
$d$ and $R$ are the "pringle" thickness and radius respectively, $\alpha$ and
$\beta$ are the two curvature parameters, and $J_n$ is the n\ :sup:`th` order
Bessel function of the first kind.

.. figure:: img/pringles_fig1.png

    Schematic of model shape (Graphic from Matt Henderson, matt@matthen.com)

Reference
---------

#. Karen Edler, Universtiy of Bath, Private Communication. 2012.
   Derivation by Stefan Alexandru Rautu.
#. L. Onsager, *Ann. New York Acad. Sci.*, 51 (1949) 627-659

Authorship and Verification
----------------------------

* **Author:** Andrew Jackson **Date:** 2008
* **Last Modified by:** Wojciech Wpotrzebowski **Date:** March 20, 2016
* **Last Reviewed by:** Andrew Jackson **Date:** September 26, 2016
"""

import numpy as np
from numpy import inf

name = "pringle"
title = "The Pringle model provides the form factor, $P(q)$, for a 'pringle' \
or 'saddle-shaped' disc that is bent in two directions."
description = """\

"""
category = "shape:cylinder"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["radius",      "Ang",         60.0,   [0, inf],    "volume", "Pringle radius"],
    ["thickness",   "Ang",         10.0,   [0, inf],    "volume", "Thickness of pringle"],
    ["alpha",       "",            0.001,  [-inf, inf], "volume", "Curvature parameter alpha"],
    ["beta",        "",            0.02,   [-inf, inf], "volume", "Curvature paramter beta"],
    ["sld", "1e-6/Ang^2",  1.0,    [-inf, inf], "sld", "Pringle sld"],
    ["sld_solvent", "1e-6/Ang^2",  6.3,    [-inf, inf], "sld", "Solvent sld"]
    ]
# pylint: enable=bad-whitespace, line-too-long


source = ["lib/polevl.c", "lib/sas_J0.c", "lib/sas_J1.c",
          "lib/sas_JN.c", "lib/gauss76.c", "pringle.c"]
radius_effective_modes = [
    "equivalent cylinder excluded volume",
    "equivalent volume sphere",
    "radius"]
has_shape_visualization = True

def create_shape_mesh(params, resolution=50):
    import numpy as np
    radius = params.get('radius', 60)
    thickness = params.get('thickness', 10)
    alpha = params.get('alpha', 0.001)
    beta = params.get('beta', 0.02)

    # Create saddle surface (hyperbolic paraboloid approximation)
    r = np.linspace(0, radius, resolution//2)
    theta = np.linspace(0, 2*np.pi, resolution)
    r_mesh, theta_mesh = np.meshgrid(r, theta)

    x = r_mesh * np.cos(theta_mesh)
    y = r_mesh * np.sin(theta_mesh)

    # Saddle shape: z = alpha * x^2 - beta * y^2
    z = alpha * x**2 - beta * y**2

    # Top and bottom surfaces (offset by thickness)
    z_top = z + thickness/2
    z_bottom = z - thickness/2

    return {
        'top_surface': (x, y, z_top),
        'bottom_surface': (x, y, z_bottom)
    }

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    import numpy as np
    radius = params.get('radius', 60)
    thickness = params.get('thickness', 10)
    alpha = params.get('alpha', 0.001)
    beta = params.get('beta', 0.02)

    # XY plane (top view) - circular outline
    theta = np.linspace(0, 2*np.pi, 100)
    circle_x = radius * np.cos(theta)
    circle_y = radius * np.sin(theta)

    ax_xy.plot(circle_x, circle_y, 'b-', linewidth=2, label='Pringle edge')
    ax_xy.fill(circle_x, circle_y, 'lightblue', alpha=0.3)
    ax_xy.set_xlim(-radius*1.3, radius*1.3)
    ax_xy.set_ylim(-radius*1.3, radius*1.3)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title('XY Cross-section (Top View)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)

    # XZ plane (side view) - parabolic curve
    x_profile = np.linspace(-radius, radius, 100)
    z_profile = alpha * x_profile**2

    ax_xz.plot(x_profile, z_profile + thickness/2, 'r-', linewidth=2, label='Top surface')
    ax_xz.plot(x_profile, z_profile - thickness/2, 'r-', linewidth=2, label='Bottom surface')
    ax_xz.fill_between(x_profile, z_profile - thickness/2, z_profile + thickness/2,
                      alpha=0.3, color='lightcoral')

    max_z = alpha * radius**2 + thickness
    ax_xz.set_xlim(-radius*1.2, radius*1.2)
    ax_xz.set_ylim(-max_z*1.5, max_z*1.5)
    ax_xz.set_xlabel('X (Å)')
    ax_xz.set_ylabel('Z (Å)')
    ax_xz.set_title('XZ Cross-section (Parabolic)')
    ax_xz.grid(True, alpha=0.3)
    ax_xz.text(0, max_z*1.3, f'α = {alpha:.4f}', ha='center', fontsize=9)

    # YZ plane (side view) - parabolic curve (downward)
    y_profile = np.linspace(-radius, radius, 100)
    z_profile_y = -beta * y_profile**2

    ax_yz.plot(y_profile, z_profile_y + thickness/2, 'g-', linewidth=2, label='Top surface')
    ax_yz.plot(y_profile, z_profile_y - thickness/2, 'g-', linewidth=2, label='Bottom surface')
    ax_yz.fill_between(y_profile, z_profile_y - thickness/2, z_profile_y + thickness/2,
                      alpha=0.3, color='lightgreen')

    max_z_y = beta * radius**2 + thickness
    ax_yz.set_xlim(-radius*1.2, radius*1.2)
    ax_yz.set_ylim(-max_z_y*1.5, max_z_y*1.5)
    ax_yz.set_xlabel('Y (Å)')
    ax_yz.set_ylabel('Z (Å)')
    ax_yz.set_title('YZ Cross-section (Hyperbolic)')
    ax_yz.grid(True, alpha=0.3)
    ax_yz.text(0, -max_z_y*1.3, f'β = {beta:.4f}', ha='center', fontsize=9)

def random():
    """Return a random parameter set for the model."""
    alpha, beta = 10**np.random.uniform(-1, 1, size=2)
    radius = 10**np.random.uniform(1, 3)
    thickness = 10**np.random.uniform(0.7, 2)
    pars = dict(
        radius=radius,
        thickness=thickness,
        alpha=alpha,
        beta=beta,
    )
    return pars

tests = [
    [{'scale' : 1.0,
      'radius': 60.0,
      'thickness': 10.0,
      'alpha': 0.001,
      'beta': 0.02,
      'sld': 1.0,
      'sld_solvent': 6.3,
      'background': 0.001,
     }, 0.1, 9.87676],

    [{'scale' : 1.0,
      'radius': 60.0,
      'thickness': 10.0,
      'alpha': 0.001,
      'beta': 0.02,
      'sld': 1.0,
      'sld_solvent': 6.3,
      'background': 0.001,
     }, 0.01, 290.56723],

    [{'scale' : 1.0,
      'radius': 60.0,
      'thickness': 10.0,
      'alpha': 0.001,
      'beta': 0.02,
      'sld': 1.0,
      'sld_solvent': 6.3,
      'background': 0.001,
     }, 0.001, 317.40847],
]
