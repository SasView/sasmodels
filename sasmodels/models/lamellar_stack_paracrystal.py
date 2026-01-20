# Note: model title and parameter table are inserted automatically
r"""
This model calculates the scattering from a stack of repeating lamellar
structures. The stacks of lamellae (infinite in lateral dimension) are
treated as a paracrystal to account for the repeating spacing. The repeat
distance is further characterized by a Gaussian polydispersity. **This model
can be used for large multilamellar vesicles.**

Definition
----------

In the equations below,

- *scale* is used instead of the mass per area of the bilayer $\Gamma_m$
  (this corresponds to the volume fraction of the material in the bilayer,
  *not* the total excluded volume of the paracrystal),

- *sld* $-$ *sld_solvent* is the contrast $\Delta \rho$,

- *thickness* is the layer thickness $t$,

- *Nlayers* is the number of layers $N$,

- *d_spacing* is the average distance between adjacent layers
  $\langle D \rangle$, and

- *sigma_d* is the relative standard deviation of the Gaussian
  layer distance distribution $\sigma_D / \langle D \rangle$.


The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = 2\pi\Delta\rho^2\Gamma_m\frac{P_\text{bil}(q)}{q^2} Z_N(q)

The form factor of the bilayer is approximated as the cross section of an
infinite, planar bilayer of thickness $t$ (compare the equations for the
lamellar model).

.. math::

    P_\text{bil}(q) = \left(\frac{\sin(qt/2)}{qt/2}\right)^2

$Z_N(q)$ describes the interference effects for aggregates
consisting of more than one bilayer. The equations used are (3-5)
from the Bergstrom reference:

.. math::


    Z_N(q) = \frac{1 - w^2}{1 + w^2 - 2w \cos(q \langle D \rangle)}
        + x_N S_N + (1 - x_N) S_{N+1}

where

.. math::

    S_N(q) = \frac{a_N}{N[1 + w^2 - 2 w \cos(q \langle D \rangle)]^2}

and

.. math::

    a_N &= 4w^2 - 2(w^3 + w) \cos(q \langle D \rangle) \\
        &\quad - 4w^{N+2}\cos(Nq \langle D \rangle)
        + 2 w^{N+3}\cos[(N-1)q \langle D \rangle]
        + 2w^{N+1}\cos[(N+1)q \langle D \rangle]

for the layer spacing distribution $w = \exp(-\sigma_D^2 q^2/2)$.

Non-integer numbers of stacks are calculated as a linear combination of
the lower and higher values

.. math::

    N_L = x_N N + (1 - x_N)(N+1)

The 2D scattering intensity is the same as 1D, regardless of the orientation
of the $q$ vector which is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


Reference
---------

#. M Bergstrom, J S Pedersen, P Schurtenberger, S U Egelhaaf,
   *J. Phys. Chem. B*, 103 (1999) 9888-9897

Authorship and Verification
----------------------------

* **Author:**
* **Last Modified by:**
* **Last Reviewed by:** Oliver Hammond, **Date:**  January 9, 2025
"""
import numpy as np
from numpy import inf

name = "lamellar_stack_paracrystal"
title = "Random lamellar sheet with paracrystal structure factor"
description = """\
    [Random lamellar phase with paracrystal structure factor]
        randomly oriented stacks of infinite sheets
        with paracrystal S(Q), having polydisperse spacing.
        sld = sheet scattering length density
        sld_solvent = solvent scattering length density
        background = incoherent background
        scale = scale factor
"""
category = "shape:lamellae"

has_shape_visualization = True
single = False

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for stacked lamellar paracrystal visualization."""
    import numpy as np
    thickness = params.get('thickness', 33)
    Nlayers = int(params.get('Nlayers', 20))
    d_spacing = params.get('d_spacing', 250)
    
    # Limit layers for visualization
    n_vis = min(Nlayers, 5)
    sheet_size = d_spacing * 2
    
    x_grid = np.linspace(-sheet_size/2, sheet_size/2, resolution)
    y_grid = np.linspace(-sheet_size/2, sheet_size/2, resolution)
    x_mesh, y_mesh = np.meshgrid(x_grid, y_grid)
    
    mesh_data = {}
    
    # Create stacked layers
    total_height = (n_vis - 1) * d_spacing
    start_z = -total_height / 2
    
    for i in range(n_vis):
        z_center = start_z + i * d_spacing
        z_top = np.full_like(x_mesh, z_center + thickness/2)
        z_bottom = np.full_like(x_mesh, z_center - thickness/2)
        mesh_data[f'layer_{i}_top'] = (x_mesh, y_mesh, z_top)
        mesh_data[f'layer_{i}_bottom'] = (x_mesh, y_mesh, z_bottom)
    
    return mesh_data

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of stacked lamellar structure."""
    import numpy as np
    thickness = params.get('thickness', 33)
    Nlayers = int(params.get('Nlayers', 20))
    d_spacing = params.get('d_spacing', 250)
    
    n_vis = min(Nlayers, 5)
    sheet_size = d_spacing * 2
    
    # XY plane (top view of one layer)
    rect_x = [-sheet_size/2, sheet_size/2, sheet_size/2, -sheet_size/2, -sheet_size/2]
    rect_y = [-sheet_size/2, -sheet_size/2, sheet_size/2, sheet_size/2, -sheet_size/2]
    
    ax_xy.plot(rect_x, rect_y, 'b-', linewidth=2)
    ax_xy.fill(rect_x, rect_y, 'lightblue', alpha=0.3)
    ax_xy.set_xlim(-sheet_size*0.7, sheet_size*0.7)
    ax_xy.set_ylim(-sheet_size*0.7, sheet_size*0.7)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title(f'XY Cross-section ({Nlayers} layers)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)
    
    # XZ/YZ planes - show stacked layers
    total_height = (n_vis - 1) * d_spacing
    start_z = -total_height / 2
    colors = ['lightblue', 'lightcoral', 'lightgreen', 'lightyellow', 'plum']
    
    for i in range(n_vis):
        z_center = start_z + i * d_spacing
        rect_z = [z_center - thickness/2, z_center - thickness/2, 
                  z_center + thickness/2, z_center + thickness/2, 
                  z_center - thickness/2]
        color = colors[i % len(colors)]
        
        ax_xz.plot(rect_x, rect_z, 'b-', linewidth=1)
        ax_xz.fill(rect_x, rect_z, color, alpha=0.5)
        
        ax_yz.plot(rect_y, rect_z, 'b-', linewidth=1)
        ax_yz.fill(rect_y, rect_z, color, alpha=0.5)
    
    max_z = total_height/2 + thickness * 2
    ax_xz.set_xlim(-sheet_size*0.7, sheet_size*0.7)
    ax_xz.set_ylim(-max_z, max_z)
    ax_xz.set_xlabel('X (Å)')
    ax_xz.set_ylabel('Z (Å)')
    ax_xz.set_title(f'XZ Cross-section (d={d_spacing:.0f}Å)')
    ax_xz.grid(True, alpha=0.3)
    
    ax_yz.set_xlim(-sheet_size*0.7, sheet_size*0.7)
    ax_yz.set_ylim(-max_z, max_z)
    ax_yz.set_xlabel('Y (Å)')
    ax_yz.set_ylabel('Z (Å)')
    ax_yz.set_title('YZ Cross-section')
    ax_yz.grid(True, alpha=0.3)

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["thickness", "Ang", 33.0, [0, inf], "volume",
               "sheet thickness"],
              ["Nlayers", "", 20, [1, inf], "",
               "Number of layers"],
              ["d_spacing", "Ang", 250., [0.0, inf], "",
               "lamellar spacing of paracrystal stack"],
              ["sigma_d", "Ang", 0.0, [0.0, inf], "",
               "Sigma (polydispersity) of the lamellar spacing"],
              ["sld", "1e-6/Ang^2", 1.0, [-inf, inf], "sld",
               "layer scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 6.34, [-inf, inf], "sld",
               "Solvent scattering length density"],
             ]


source = ["lamellar_stack_paracrystal.c"]

form_volume = """
    return 1.0;
"""

def random():
    """Return a random parameter set for the model."""
    total_thickness = 10**np.random.uniform(2, 4.7)
    Nlayers = np.random.randint(2, 200)
    d_spacing = total_thickness / Nlayers
    thickness = d_spacing * np.random.uniform(0, 1)
    # Let polydispersity peak around 15%; 95% < 0.4; max=100%
    sigma_d = np.random.beta(1.5, 7)
    pars = dict(
        thickness=thickness,
        Nlayers=Nlayers,
        d_spacing=d_spacing,
        sigma_d=sigma_d,
    )
    return pars

#
tests = [
    [{'scale': 1.0, 'background': 0.0, 'thickness': 33., 'Nlayers': 20.0,
      'd_spacing': 250., 'sigma_d': 0.2, 'sld': 1.0,
      'sld_solvent': 6.34, 'thickness_pd': 0.0, 'thickness_pd_n': 40},
     [0.001, 0.215268], [21829.3, 0.00487686]],
]
# ADDED by: RKH  ON: 18Mar2016  converted from sasview previously,
# now renaming everything & sorting the docs
