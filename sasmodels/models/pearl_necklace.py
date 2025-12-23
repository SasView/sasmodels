r"""
This model provides the form factor for a pearl necklace composed of two
elements: *N* pearls (homogeneous spheres of radius *R*) freely jointed by *M*
rods (like strings - with a total mass *Mw* = *M* \* *m*\ :sub:`r` + *N* \* *m*\
:sub:`s`, and the string segment length (or edge separation) *l*
(= *A* - 2\ *R*)). *A* is the center-to-center pearl separation distance.

.. figure:: img/pearl_necklace_geometry.jpg

    Pearl Necklace schematic

Definition
----------

The output of the scattering intensity function for the pearl_necklace is
given by (Schweins, 2004)

.. math::

    I(q)=\frac{ \text{scale} }{V} \cdot \frac{(S_{ss}(q)+S_{ff}(q)+S_{fs}(q))}
        {(M \cdot m_f + N \cdot m_s)^2} + \text{bkg}

where

.. math::

    S_{ss}(q) &= 2m_s^2\psi^2(q)\left[\frac{N}{1-sin(qA)/qA}-\frac{N}{2}-
        \frac{1-(sin(qA)/qA)^N}{(1-sin(qA)/qA)^2}\cdot\frac{sin(qA)}{qA}\right] \\
    S_{ff}(q) &= m_r^2\left[M\left\{2\Lambda(q)-\left(\frac{sin(ql/2)}{ql/2}\right)\right\}+
        \frac{2M\beta^2(q)}{1-sin(qA)/qA}-2\beta^2(q)\cdot
        \frac{1-(sin(qA)/qA)^M}{(1-sin(qA)/qA)^2}\right] \\
    S_{fs}(q) &= m_r \beta (q) \cdot m_s \psi (q) \cdot 4\left[
        \frac{N-1}{1-sin(qA)/qA}-\frac{1-(sin(qA)/qA)^{N-1}}{(1-sin(qA)/qA)^2}
        \cdot \frac{sin(qA)}{qA}\right] \\
    \psi(q) &= 3 \cdot \frac{sin(qR)-(qR)\cdot cos(qR)}{(qR)^3} \\
    \Lambda(q) &= \frac{\int_0^{ql}\frac{sin(t)}{t}dt}{ql} \\
    \beta(q) &= \frac{\int_{qR}^{q(A-R)}\frac{sin(t)}{t}dt}{ql}

where the mass *m*\ :sub:`i` is (SLD\ :sub:`i` - SLD\ :sub:`solvent`) \*
(volume of the *N* pearls/rods). *V* is the total volume of the necklace.

.. note::

   *num_pearls* must be an integer.

The 2D scattering intensity is the same as $P(q)$ above, regardless of the
orientation of the *q* vector.

References
----------

#. R Schweins and K Huber, *Particle Scattering Factor of Pearl Necklace Chains*,
   *Macromol. Symp.* 211 (2004) 25-42 2004

#. L. Onsager, *Ann. New York Acad. Sci.*, 51 (1949) 627-659

Authorship and Verification
----------------------------

* **Author:**
* **Last Modified by:** Andrew Jackson **Date:** March 28, 2019
* **Last Reviewed by:** Steve King **Date:** March 28, 2019
"""

import numpy as np
from numpy import inf

name = "pearl_necklace"
title = "Colloidal spheres chained together with no preferential orientation"
description = """
Calculate form factor for Pearl Necklace Model
[Macromol. Symp. 2004, 211, 25-42]
Parameters:
background:background
scale: scale factor
sld: the SLD of the pearl spheres
sld_string: the SLD of the strings
sld_solvent: the SLD of the solvent
num_pearls: number of the pearls
radius: the radius of a pearl
edge_sep: the length of string segment; surface to surface
thick_string: thickness (ie, diameter) of the string
"""
category = "shape:cylinder"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["radius", "Ang", 80.0, [0, inf], "volume",
               "Mean radius of the chained spheres"],
              ["edge_sep", "Ang", 350.0, [0, inf], "volume",
               "Mean separation of chained particles"],
              ["thick_string", "Ang", 2.5, [0, inf], "volume",
               "Thickness of the chain linkage"],
              ["num_pearls", "none", 3, [1, inf], "volume",
               "Number of pearls in the necklace (must be integer)"],
              ["sld", "1e-6/Ang^2", 1.0, [-inf, inf], "sld",
               "Scattering length density of the chained spheres"],
              ["sld_string", "1e-6/Ang^2", 1.0, [-inf, inf], "sld",
               "Scattering length density of the chain linkage"],
              ["sld_solvent", "1e-6/Ang^2", 6.3, [-inf, inf], "sld",
               "Scattering length density of the solvent"],
             ]

source = ["lib/sas_Si.c", "lib/sas_3j1x_x.c", "pearl_necklace.c"]
valid = "thick_string < radius && num_pearls > 0.0"
single = False  # use double precision unless told otherwise
radius_effective_modes = ["equivalent volume sphere"]
has_shape_visualization = True

def create_shape_mesh(params, resolution=40):
    import numpy as np
    radius = params.get('radius', 80.0)
    edge_sep = params.get('edge_sep', 350.0)  # surface-to-surface distance
    thick_string = params.get('thick_string', 2.5)
    num_pearls = int(round(params.get('num_pearls', 3)))
    num_pearls = max(num_pearls, 1)

    # Spacing between pearl centers
    center_step = 2 * radius + edge_sep
    z_positions = [
        (i - (num_pearls - 1) / 2.0) * center_step
        for i in range(num_pearls)
    ]

    # Sphere (pearl) mesh
    phi = np.linspace(0, np.pi, resolution // 2)
    theta = np.linspace(0, 2 * np.pi, resolution)
    phi_mesh, theta_mesh = np.meshgrid(phi, theta)

    pearls = {}
    for i, z0 in enumerate(z_positions):
        x = radius * np.sin(phi_mesh) * np.cos(theta_mesh)
        y = radius * np.sin(phi_mesh) * np.sin(theta_mesh)
        z = radius * np.cos(phi_mesh) + z0
        pearls[f'pearl_{i}'] = (x, y, z)

    # String segments as thin cylinders between neighboring pearls
    string_radius = thick_string / 2.0
    strings = {}
    if num_pearls > 1 and string_radius > 0:
        theta_c = np.linspace(0, 2 * np.pi, resolution)
        for i in range(num_pearls - 1):
            z_start = z_positions[i] + radius
            z_end = z_positions[i + 1] - radius
            z_seg = np.linspace(z_start, z_end, resolution // 2)
            theta_mesh_c, z_seg_mesh = np.meshgrid(theta_c, z_seg)
            x_c = string_radius * np.cos(theta_mesh_c)
            y_c = string_radius * np.sin(theta_mesh_c)
            strings[f'string_{i}'] = (x_c, y_c, z_seg_mesh)

    mesh = {}
    mesh.update(pearls)
    mesh.update(strings)
    return mesh

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    import numpy as np
    radius = params.get('radius', 80.0)
    edge_sep = params.get('edge_sep', 350.0)
    thick_string = params.get('thick_string', 2.5)
    num_pearls = int(round(params.get('num_pearls', 3)))
    num_pearls = max(num_pearls, 1)

    center_step = 2 * radius + edge_sep
    z_positions = np.array(
        [(i - (num_pearls - 1) / 2.0) * center_step for i in range(num_pearls)]
    )

    # XY: top view of a single pearl + string cross-section
    theta = np.linspace(0, 2 * np.pi, 200)
    pearl_x = radius * np.cos(theta)
    pearl_y = radius * np.sin(theta)
    ax_xy.plot(pearl_x, pearl_y, 'b-', linewidth=2, label='Pearl')
    ax_xy.fill(pearl_x, pearl_y, 'lightblue', alpha=0.4)

    if thick_string > 0:
        string_r = thick_string / 2.0
        sx = string_r * np.cos(theta)
        sy = string_r * np.sin(theta)
        ax_xy.plot(sx, sy, 'r--', linewidth=1, label='String')
        ax_xy.fill(sx, sy, 'lightcoral', alpha=0.3)

    ax_xy.set_xlim(-radius * 1.4, radius * 1.4)
    ax_xy.set_ylim(-radius * 1.4, radius * 1.4)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title('XY Cross-section (Single pearl + string)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)
    ax_xy.legend()

    # XZ: chain of circles along Z
    for z0 in z_positions:
        circle_z = radius * np.sin(theta)
        circle_x = radius * np.cos(theta)
        ax_xz.plot(z0 + circle_z, circle_x, 'b-', alpha=0.7)

    # Draw string as line along centers
    ax_xz.plot(z_positions, np.zeros_like(z_positions), 'r-', linewidth=2, label='String axis')
    ax_xz.set_xlabel('Z (Å)')
    ax_xz.set_ylabel('X (Å)')
    ax_xz.set_title('XZ Cross-section (Chain of pearls)')
    ax_xz.grid(True, alpha=0.3)
    ax_xz.legend()

    # YZ: same as XZ but in Y
    for z0 in z_positions:
        circle_z = radius * np.sin(theta)
        circle_y = radius * np.cos(theta)
        ax_yz.plot(z0 + circle_z, circle_y, 'g-', alpha=0.7)
    ax_yz.plot(z_positions, np.zeros_like(z_positions), 'r-', linewidth=2, label='String axis')
    ax_yz.set_xlabel('Z (Å)')
    ax_yz.set_ylabel('Y (Å)')
    ax_yz.set_title('YZ Cross-section (Chain of pearls)')
    ax_yz.grid(True, alpha=0.3)
    ax_yz.legend()

def random():
    """Return a random parameter set for the model."""
    radius = 10**np.random.uniform(1, 3) # 1 - 1000
    thick_string = 10**np.random.uniform(0, np.log10(radius)-1) # 1 - radius/10
    edge_sep = 10**np.random.uniform(0, 3)  # 1 - 1000
    num_pearls = np.round(10**np.random.uniform(0.3, 3)) # 2 - 1000
    pars = dict(
        radius=radius,
        edge_sep=edge_sep,
        thick_string=thick_string,
        num_pearls=num_pearls,
    )
    return pars

# ER function is not being used here, not that it is likely very sensible to
# include an S(Q) with this model, the default in sasview 5.0 would be to the
# "unconstrained" radius_effective.
#tests = [[{}, 0.001, 17380.245], [{}, 'ER', 115.39502]]
tests = [[{}, 0.001, 17380.245]]
