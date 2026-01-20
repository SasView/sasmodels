r"""
Polydispersity in the bilayer thickness can be applied from the GUI.

Definition
----------

The scattering intensity $I(q)$ for dilute, randomly oriented,
"infinitely large" sheets or lamellae is

.. math::

    I(q) = \text{scale}\frac{2\pi P(q)}{q^2\delta} + \text{background}


The form factor is

.. math::

   P(q) = \frac{2\Delta\rho^2}{q^2}(1-\cos(q\delta))
        = \frac{4\Delta\rho^2}{q^2}\sin^2\left(\frac{q\delta}{2}\right)

where $\delta$ is the total layer thickness and $\Delta\rho$ is the
scattering length density difference.

This is the limiting form for a spherical shell of infinitely large radius.
Note that the division by $\delta$ means that $scale$ in sasview is the
volume fraction of sheet, $\phi = S\delta$ where $S$ is the area of sheet
per unit volume. $S$ is half the Porod surface area per unit volume of a
thicker layer (as that would include both faces of the sheet).

The 2D scattering intensity is calculated in the same way as 1D, where
the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


References
----------

#. F Nallet, R Laversanne, and D Roux, *J. Phys. II France*, 3, (1993) 487-502
#. J Berghausen, J Zipfel, P Lindner, W Richtering,
   *J. Phys. Chem. B*, 105, (2001) 11081-11088

Authorship and Verification
----------------------------

* **Author:**
* **Last Modified by:**
* **Last Reviewed by:**
"""

import numpy as np
from numpy import inf

name = "lamellar"
title = "Lyotropic lamellar phase with uniform SLD and random distribution"
description = """\
    [Dilute Lamellar Form Factor](from a lyotropic lamellar phase)
        I(q)= 2*pi*P(q)/(delta *q^(2)), where
        P(q)=2*(contrast/q)^(2)*(1-cos(q*delta))^(2))
        thickness = layer thickness
        sld = layer scattering length density
        sld_solvent = solvent scattering length density
        background = incoherent background
        scale = scale factor
"""
category = "shape:lamellae"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["thickness",          "Ang", 50, [0, inf],    "volume", "total layer thickness" ],
    ["sld",         "1e-6/Ang^2",  1, [-inf, inf], "sld",    "Layer scattering length density" ],
    ["sld_solvent", "1e-6/Ang^2",  6, [-inf, inf], "sld",    "Solvent scattering length density" ],
    ]
# pylint: enable=bad-whitespace, line-too-long

# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized? - it is!
form_volume = """
    return 1.0;
    """

Iq = """
    const double sub = sld - sld_solvent;
    const double qsq = q*q;
    // Original expression
    //return 4.0e-4*M_PI*sub*sub/qsq * (1.0-cos(q*thickness)) / (thickness*qsq);
    // const double alpha = fmod(q*thickness+0.1, 2.0*M_PI)-0.1;
    // Use small angle fix 1-cos(theta) = 2 sin^2(theta/2)
    const double sinq2 = sin(0.5*q*thickness);
    return 4.0e-4*M_PI*sub*sub/qsq * 2.0*sinq2*sinq2 / (thickness*qsq);
    """

has_shape_visualization = True

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for lamellar (flat sheet) visualization."""
    import numpy as np
    thickness = params.get('thickness', 50)

    # Create a flat sheet (lamella) - represent as a thin rectangular slab
    sheet_size = thickness * 4  # Make sheet visually larger than thickness

    # Top surface
    x_top = np.array([[-sheet_size/2, sheet_size/2, sheet_size/2, -sheet_size/2, -sheet_size/2]])
    y_top = np.array([[-sheet_size/2, -sheet_size/2, sheet_size/2, sheet_size/2, -sheet_size/2]])
    z_top = np.full_like(x_top, thickness/2)

    # Bottom surface
    z_bottom = np.full_like(x_top, -thickness/2)

    # Create meshgrid for surfaces
    x_grid = np.linspace(-sheet_size/2, sheet_size/2, resolution)
    y_grid = np.linspace(-sheet_size/2, sheet_size/2, resolution)
    x_mesh, y_mesh = np.meshgrid(x_grid, y_grid)
    z_top_mesh = np.full_like(x_mesh, thickness/2)
    z_bottom_mesh = np.full_like(x_mesh, -thickness/2)

    return {
        'top_surface': (x_mesh, y_mesh, z_top_mesh),
        'bottom_surface': (x_mesh, y_mesh, z_bottom_mesh)
    }

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of the lamellar structure."""
    thickness = params.get('thickness', 50)
    sheet_size = thickness * 4

    # XY plane (top view) - shows the sheet as a square
    rect_x = [-sheet_size/2, sheet_size/2, sheet_size/2, -sheet_size/2, -sheet_size/2]
    rect_y = [-sheet_size/2, -sheet_size/2, sheet_size/2, sheet_size/2, -sheet_size/2]

    ax_xy.plot(rect_x, rect_y, 'b-', linewidth=2, label='Lamella')
    ax_xy.fill(rect_x, rect_y, 'lightblue', alpha=0.3)
    ax_xy.set_xlim(-sheet_size*0.7, sheet_size*0.7)
    ax_xy.set_ylim(-sheet_size*0.7, sheet_size*0.7)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title('XY Cross-section (Top View)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)

    # XZ plane (side view) - shows the thickness
    rect_x_side = [-sheet_size/2, sheet_size/2, sheet_size/2, -sheet_size/2, -sheet_size/2]
    rect_z_side = [-thickness/2, -thickness/2, thickness/2, thickness/2, -thickness/2]

    ax_xz.plot(rect_x_side, rect_z_side, 'r-', linewidth=2, label='Lamella')
    ax_xz.fill(rect_x_side, rect_z_side, 'lightcoral', alpha=0.3)
    ax_xz.set_xlim(-sheet_size*0.7, sheet_size*0.7)
    ax_xz.set_ylim(-thickness*2, thickness*2)
    ax_xz.set_xlabel('X (Å)')
    ax_xz.set_ylabel('Z (Å)')
    ax_xz.set_title(f'XZ Cross-section (thickness={thickness:.1f}Å)')
    ax_xz.grid(True, alpha=0.3)

    # YZ plane (front view) - same as XZ
    rect_y_front = [-sheet_size/2, sheet_size/2, sheet_size/2, -sheet_size/2, -sheet_size/2]

    ax_yz.plot(rect_y_front, rect_z_side, 'g-', linewidth=2, label='Lamella')
    ax_yz.fill(rect_y_front, rect_z_side, 'lightgreen', alpha=0.3)
    ax_yz.set_xlim(-sheet_size*0.7, sheet_size*0.7)
    ax_yz.set_ylim(-thickness*2, thickness*2)
    ax_yz.set_xlabel('Y (Å)')
    ax_yz.set_ylabel('Z (Å)')
    ax_yz.set_title('YZ Cross-section (Front View)')
    ax_yz.grid(True, alpha=0.3)

def random():
    """Return a random parameter set for the model."""
    thickness = 10**np.random.uniform(1, 4)
    pars = dict(
        thickness=thickness,
    )
    return pars

#  [(qx1, qy1), (qx2, qy2), ...], [I(qx1,qy1), I(qx2,qy2), ...]],
tests = [
    [{'scale': 1.0, 'background': 0.0, 'thickness': 50.0,
      'sld': 1.0, 'sld_solvent': 6.3, 'thickness_pd': 0.0},
     [0.001], [882289.54309]]
]
# ADDED by: converted by PAK? (or RKH?)
# ON: 16Mar2016 - RKH adding unit tests from sasview to early 2015 conversion
