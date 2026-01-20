# Note: model title and parameter table are inserted automatically
r"""
This model provides the scattering intensity, $I(q)$, for a lyotropic lamellar
phase where a random distribution in solution are assumed. The SLD of the head
region is taken to be different from the SLD of the tail region.

Definition
----------

The scattering intensity $I(q)$ is

.. math::

   I(q) = 2\pi\frac{\text{scale}}{2(\delta_H + \delta_T)}  P(q) \frac{1}{q^2}

The form factor $P(q)$ is

.. math::

    P(q) = \frac{4}{q^2}
        \left\lbrace
            \Delta \rho_H
            \left[\sin[q(\delta_H + \delta_T)\ - \sin(q\delta_T)\right]
            + \Delta\rho_T\sin(q\delta_T)
        \right\rbrace^2

where $\delta_T$ is *length_tail*, $\delta_H$ is *length_head*,
$\Delta\rho_H$ is the head contrast (*sld_head* $-$ *sld_solvent*),
and $\Delta\rho_T$ is tail contrast (*sld* $-$ *sld_solvent*).

The total thickness of the lamellar sheet is
$a_H + \delta_T + \delta_T + \delta_H$. Note that in a non aqueous solvent
the chemical "head" group may be the "Tail region" and vice-versa.

The 2D scattering intensity is calculated in the same way as 1D, where
the $q$ vector is defined as

.. math:: q = \sqrt{q_x^2 + q_y^2}


References
----------

#. F Nallet, R Laversanne, and D Roux, *J. Phys. II France*, 3, (1993) 487-502
#. J Berghausen, J Zipfel, P Lindner, W Richtering,
   *J. Phys. Chem. B*, 105, (2001) 11081-11088

Authorship and Verification
----------------------------

* **Author:**
* **Last Modified by:**
* **Last Reviewed by:** S King and P Butler **Date** April 17, 2014
"""

import numpy as np
from numpy import inf

name = "lamellar_hg"
title = "Random lamellar phase with Head and Tail Groups"
description = """\
    [Random lamellar phase with Head and Tail Groups]
        I(q)= 2*pi*P(q)/(2(H+T)*q^(2)), where
        P(q)= see manual
        layer thickness =(H+T+T+H) = 2(Head+Tail)
        sld = Tail scattering length density
        sld_head = Head scattering length density
        sld_solvent = solvent scattering length density
        background = incoherent background
        scale = scale factor
"""
category = "shape:lamellae"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["length_tail", "Ang",       15,   [0, inf],  "volume",  "Tail thickness ( total = H+T+T+H)"],
              ["length_head", "Ang",       10,   [0, inf],  "volume",  "Head thickness"],
              ["sld",         "1e-6/Ang^2", 0.4, [-inf,inf], "sld",    "Tail scattering length density"],
              ["sld_head",    "1e-6/Ang^2", 3.0, [-inf,inf], "sld",    "Head scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 6,   [-inf,inf], "sld",    "Solvent scattering length density"]]
# pylint: enable=bad-whitespace, line-too-long

# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return 1.0;
    """

source = ["lamellar_hg.c"]

has_shape_visualization = True

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for lamellar with head groups visualization."""
    import numpy as np
    length_tail = params.get('length_tail', 15)
    length_head = params.get('length_head', 10)

    total_thickness = 2 * (length_head + length_tail)  # H+T+T+H
    sheet_size = total_thickness * 3

    # Create meshgrid for surfaces
    x_grid = np.linspace(-sheet_size/2, sheet_size/2, resolution)
    y_grid = np.linspace(-sheet_size/2, sheet_size/2, resolution)
    x_mesh, y_mesh = np.meshgrid(x_grid, y_grid)

    # Outer surfaces (head groups)
    z_top_outer = np.full_like(x_mesh, total_thickness/2)
    z_bottom_outer = np.full_like(x_mesh, -total_thickness/2)

    # Inner surfaces (between head and tail)
    z_top_inner = np.full_like(x_mesh, length_tail)  # Top head/tail interface
    z_bottom_inner = np.full_like(x_mesh, -length_tail)  # Bottom head/tail interface

    return {
        'top_head': (x_mesh, y_mesh, z_top_outer),
        'top_interface': (x_mesh, y_mesh, z_top_inner),
        'bottom_interface': (x_mesh, y_mesh, z_bottom_inner),
        'bottom_head': (x_mesh, y_mesh, z_bottom_outer)
    }

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of the lamellar with head groups."""
    length_tail = params.get('length_tail', 15)
    length_head = params.get('length_head', 10)

    total_thickness = 2 * (length_head + length_tail)
    sheet_size = total_thickness * 3

    # XY plane (top view)
    rect_x = [-sheet_size/2, sheet_size/2, sheet_size/2, -sheet_size/2, -sheet_size/2]
    rect_y = [-sheet_size/2, -sheet_size/2, sheet_size/2, sheet_size/2, -sheet_size/2]

    ax_xy.plot(rect_x, rect_y, 'b-', linewidth=2)
    ax_xy.fill(rect_x, rect_y, 'lightblue', alpha=0.3)
    ax_xy.set_xlim(-sheet_size*0.7, sheet_size*0.7)
    ax_xy.set_ylim(-sheet_size*0.7, sheet_size*0.7)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title('XY Cross-section (Top View)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)

    # XZ plane (side view) - shows head and tail regions
    # Top head group
    head_top = [[-sheet_size/2, sheet_size/2, sheet_size/2, -sheet_size/2, -sheet_size/2],
                [length_tail, length_tail, total_thickness/2, total_thickness/2, length_tail]]
    # Tail region
    tail = [[-sheet_size/2, sheet_size/2, sheet_size/2, -sheet_size/2, -sheet_size/2],
            [-length_tail, -length_tail, length_tail, length_tail, -length_tail]]
    # Bottom head group
    head_bottom = [[-sheet_size/2, sheet_size/2, sheet_size/2, -sheet_size/2, -sheet_size/2],
                   [-total_thickness/2, -total_thickness/2, -length_tail, -length_tail, -total_thickness/2]]

    ax_xz.fill(head_top[0], head_top[1], 'coral', alpha=0.5, label='Head')
    ax_xz.fill(tail[0], tail[1], 'lightblue', alpha=0.5, label='Tail')
    ax_xz.fill(head_bottom[0], head_bottom[1], 'coral', alpha=0.5)
    ax_xz.plot([-sheet_size/2, sheet_size/2], [total_thickness/2, total_thickness/2], 'r-', linewidth=2)
    ax_xz.plot([-sheet_size/2, sheet_size/2], [-total_thickness/2, -total_thickness/2], 'r-', linewidth=2)
    ax_xz.set_xlim(-sheet_size*0.7, sheet_size*0.7)
    ax_xz.set_ylim(-total_thickness, total_thickness)
    ax_xz.set_xlabel('X (Å)')
    ax_xz.set_ylabel('Z (Å)')
    ax_xz.set_title('XZ Cross-section (Side View)')
    ax_xz.legend(loc='upper right', fontsize=8)
    ax_xz.grid(True, alpha=0.3)

    # YZ plane (front view)
    ax_yz.fill(head_top[0], head_top[1], 'coral', alpha=0.5)
    ax_yz.fill(tail[0], tail[1], 'lightblue', alpha=0.5)
    ax_yz.fill(head_bottom[0], head_bottom[1], 'coral', alpha=0.5)
    ax_yz.plot([-sheet_size/2, sheet_size/2], [total_thickness/2, total_thickness/2], 'r-', linewidth=2)
    ax_yz.plot([-sheet_size/2, sheet_size/2], [-total_thickness/2, -total_thickness/2], 'r-', linewidth=2)
    ax_yz.set_xlim(-sheet_size*0.7, sheet_size*0.7)
    ax_yz.set_ylim(-total_thickness, total_thickness)
    ax_yz.set_xlabel('Y (Å)')
    ax_yz.set_ylabel('Z (Å)')
    ax_yz.set_title('YZ Cross-section (Front View)')
    ax_yz.grid(True, alpha=0.3)

def random():
    """Return a random parameter set for the model."""
    thickness = 10**np.random.uniform(1, 4)
    length_head = thickness * np.random.uniform(0, 1)
    length_tail = thickness - length_head
    pars = dict(
        length_head=length_head,
        length_tail=length_tail,
    )
    return pars

#
tests = [
    [{'scale': 1.0, 'background': 0.0, 'length_tail': 15.0, 'length_head': 10.0,
      'sld': 0.4, 'sld_head': 3.0, 'sld_solvent': 6.0},
     [0.001], [653143.9209]],
]
# ADDED by: RKH  ON: 18Mar2016  converted from sasview previously, now renaming everything & sorting the docs
