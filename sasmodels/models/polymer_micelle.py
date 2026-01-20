r"""

This model provides the form factor, $P(q)$, for a micelle with a spherical
core and Gaussian polymer chains attached to the surface, thus may be applied
to block copolymer micelles. Take special care before use and read both of the
provided references carefully; this model typically only works well when the
Gaussian chains are much smaller than the core, which is often not the case.

Definition
----------

The 1D scattering intensity for this model is calculated according to
the equations given by Pedersen (1996, 2000), summarised briefly here.

The micelle core is imagined as $N$ = *n_aggreg* polymer heads, each of volume
$V_\text{core}$, which then defines a micelle core of radius $r$ = *r_core*,
which is a separate parameter even though it could be directly determined.
The Gaussian random coil tails, of gyration radius $R_g$, are imagined
uniformly distributed around the spherical core, centred at a distance
$r + d \cdot R_g$ from the micelle centre, where $d$ = *d_penetration* is
of order unity. A volume $V_\text{corona}$ is defined for each coil. The
model in detail seems to separately parameterize the terms for the shape
of $I(Q)$ and the relative intensity of each term, so use with caution
and check parameters for consistency. The spherical core is monodisperse,
so it's intensity and the cross terms may have sharp oscillations (use $q$
resolution smearing if needs be to help remove them).

.. math::
    P(q) &= N^2\beta^2_s\Phi(qr)^2 + N\beta^2_cP_c(q)
            + 2N^2\beta_s\beta_cS_{sc}(q) + N(N-1)\beta_c^2S_{cc}(q) \\
    \beta_s &= V_\text{core}(\rho_\text{core} - \rho_\text{solvent}) \\
    \beta_c &= V_\text{corona}(\rho_\text{corona} - \rho_\text{solvent})

where $\rho_\text{core}$, $\rho_\text{corona}$ and $\rho_\text{solvent}$ are
the scattering length densities *sld_core*, *sld_corona* and *sld_solvent*.
For the spherical core of radius $r$

.. math::
   \Phi(qr)= \frac{\sin(qr) - qr\cos(qr)}{(qr)^3}

whilst for the Gaussian coils

.. math::

   P_c(q) &= 2 [\exp(-Z) + Z - 1] / Z^2 \\
   Z &= (q R_g)^2

The sphere to coil (core to corona) and coil to coil (corona to corona) cross
terms are approximated by:

.. math::

   S_{sc}(q) &= \Phi(qr)\psi(Z)
       \frac{\sin(q(r+d \cdot R_g))}{q(r+d \cdot R_g)} \\
   S_{cc}(q) &= \psi(Z)^2
       \left[\frac{\sin(q(r+d \cdot R_g))}{q(r+d \cdot R_g)} \right]^2 \\
   \psi(Z) &= \frac{[1-\exp^{-Z}]}{Z}

Validation
----------

$P(q)$ above is multiplied by *ndensity*, and a units conversion of $10^{-13}$,
so *scale* is likely 1.0 if the scattering data is in absolute units. This
model has not yet been independently validated.


References
----------

#.  J Pedersen and M Gerstenberg, *Macromolecules*, 29 (1996) 1363–1365
#.  J Pedersen, *J. Appl. Cryst.*, 33 (2000) 637-640


Authorship and Verification
----------------------------

* **Translated by   :** Richard Heenan **Date:** March 20, 2016
* **Last modified by:** Paul Kienzle **Date:** November 29, 2017
* **Last reviewed by:** Oliver Hammond **Date:** January 13, 2024
"""

import numpy as np
from numpy import inf, pi

name = "polymer_micelle"
title = "Polymer micelle model"
description = """
This model provides the form factor, $P(q)$, for a micelle with a spherical
core and Gaussian polymer chains attached to the surface, thus may be applied
to block copolymer micelles. Take special care before use and read both of the
provided references carefully; this model typically only works well when the
Gaussian chains are much smaller than the core, which is often not the case.
    """


category = "shape:sphere"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["ndensity",      "1e15/cm^3",  8.94, [0.0, inf], "", "Number density of micelles"],
    ["v_core",        "Ang^3",  62624.0,  [0.0, inf], "", "Core volume "],
    ["v_corona",      "Ang^3",  61940.0,  [0.0, inf], "", "Corona volume"],
    ["sld_solvent",   "1e-6/Ang^2", 6.4,  [0.0, inf], "sld", "Solvent scattering length density"],
    ["sld_core",      "1e-6/Ang^2", 0.34, [0.0, inf], "sld", "Core scattering length density"],
    ["sld_corona",    "1e-6/Ang^2", 0.8,  [0.0, inf], "sld", "Corona scattering length density"],
    ["radius_core",   "Ang",       45.0,  [0.0, inf], "", "Radius of core ( must be >> rg )"],
    ["rg",    "Ang",       20.0,  [0.0, inf], "", "Radius of gyration of chains in corona"],
    ["d_penetration", "",           1.0,  [-inf, inf], "", "Factor to mimic non-penetration of Gaussian chains"],
    ["n_aggreg",      "",           6.0,  [-inf, inf], "", "Aggregation number of the micelle"],
    ]
# pylint: enable=bad-whitespace, line-too-long

has_shape_visualization = True
single = False

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for polymer micelle visualization."""
    import numpy as np
    radius_core = params.get('radius_core', 45)
    rg = params.get('rg', 20)  # radius of gyration of corona chains
    d_penetration = params.get('d_penetration', 1.0)
    
    phi = np.linspace(0, np.pi, resolution//2)
    theta = np.linspace(0, 2*np.pi, resolution)
    phi_mesh, theta_mesh = np.meshgrid(phi, theta)
    
    mesh_data = {}
    
    # Core sphere
    x_core = radius_core * np.sin(phi_mesh) * np.cos(theta_mesh)
    y_core = radius_core * np.sin(phi_mesh) * np.sin(theta_mesh)
    z_core = radius_core * np.cos(phi_mesh)
    mesh_data['core'] = (x_core, y_core, z_core)
    
    # Corona (fuzzy outer region) - represented as a larger transparent sphere
    corona_radius = radius_core + d_penetration * rg * 2
    x_corona = corona_radius * np.sin(phi_mesh) * np.cos(theta_mesh)
    y_corona = corona_radius * np.sin(phi_mesh) * np.sin(theta_mesh)
    z_corona = corona_radius * np.cos(phi_mesh)
    mesh_data['corona'] = (x_corona, y_corona, z_corona)
    
    return mesh_data

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of polymer micelle."""
    import numpy as np
    radius_core = params.get('radius_core', 45)
    rg = params.get('rg', 20)
    d_penetration = params.get('d_penetration', 1.0)
    
    theta = np.linspace(0, 2*np.pi, 100)
    corona_radius = radius_core + d_penetration * rg * 2
    max_r = corona_radius * 1.3
    
    # Core circle
    core_x = radius_core * np.cos(theta)
    core_y = radius_core * np.sin(theta)
    
    # Corona circle
    corona_x = corona_radius * np.cos(theta)
    corona_y = corona_radius * np.sin(theta)
    
    # XY plane
    ax_xy.fill(corona_x, corona_y, 'lightyellow', alpha=0.3, label='Corona')
    ax_xy.plot(corona_x, corona_y, 'orange', linewidth=2, linestyle='--')
    ax_xy.fill(core_x, core_y, 'lightblue', alpha=0.5, label='Core')
    ax_xy.plot(core_x, core_y, 'b-', linewidth=2)
    ax_xy.set_xlim(-max_r, max_r)
    ax_xy.set_ylim(-max_r, max_r)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title('XY Cross-section')
    ax_xy.set_aspect('equal')
    ax_xy.legend(loc='upper right', fontsize=8)
    ax_xy.grid(True, alpha=0.3)
    
    # XZ plane
    ax_xz.fill(corona_x, corona_y, 'lightyellow', alpha=0.3)
    ax_xz.plot(corona_x, corona_y, 'orange', linewidth=2, linestyle='--')
    ax_xz.fill(core_x, core_y, 'lightblue', alpha=0.5)
    ax_xz.plot(core_x, core_y, 'b-', linewidth=2)
    ax_xz.set_xlim(-max_r, max_r)
    ax_xz.set_ylim(-max_r, max_r)
    ax_xz.set_xlabel('X (Å)')
    ax_xz.set_ylabel('Z (Å)')
    ax_xz.set_title('XZ Cross-section')
    ax_xz.set_aspect('equal')
    ax_xz.grid(True, alpha=0.3)
    
    # YZ plane
    ax_yz.fill(corona_x, corona_y, 'lightyellow', alpha=0.3)
    ax_yz.plot(corona_x, corona_y, 'orange', linewidth=2, linestyle='--')
    ax_yz.fill(core_x, core_y, 'lightblue', alpha=0.5)
    ax_yz.plot(core_x, core_y, 'b-', linewidth=2)
    ax_yz.set_xlim(-max_r, max_r)
    ax_yz.set_ylim(-max_r, max_r)
    ax_yz.set_xlabel('Y (Å)')
    ax_yz.set_ylabel('Z (Å)')
    ax_yz.set_title('YZ Cross-section')
    ax_yz.set_aspect('equal')
    ax_yz.grid(True, alpha=0.3)

source = ["lib/sas_3j1x_x.c", "polymer_micelle.c"]

def random():
    """Return a random parameter set for the model."""
    radius_core = 10**np.random.uniform(1, 3)
    rg = radius_core * 10**np.random.uniform(-2, -0.3)
    d_penetration = np.random.randn()*0.05 + 1
    n_aggreg = np.random.randint(3, 30)
    # volume of head groups is the core volume over the number of groups,
    # with a correction for packing fraction of the head groups.
    v_core = 4*pi/3*radius_core**3/n_aggreg * 0.68
    # Rg^2 for gaussian coil is a^2n/6 => a^2 = 6 Rg^2/n
    # a=2r => r = Rg sqrt(3/2n)
    # v = 4/3 pi r^3 n => v = 4/3 pi Rg^3 (3/2n)^(3/2) n = pi Rg^3 sqrt(6/n)
    tail_segments = np.random.randint(6, 30)
    v_corona = pi * rg**3 * np.sqrt(6/tail_segments)
    V = 4*pi/3*(radius_core + rg)**3
    pars = dict(
        background=0,
        scale=1e7/V,
        ndensity=8.94,
        v_core=v_core,
        v_corona=v_corona,
        radius_core=radius_core,
        rg=rg,
        d_penetration=d_penetration,
        n_aggreg=n_aggreg,
    )
    return pars

tests = [
    [{}, 0.01, 15.3532],
    ]
# RKH 20Mar2016 - need to check whether the core & corona volumes are per
#                 monomer ??? and how aggregation number works!
# renamed from micelle_spherical_core to polymer_micelle,
# moved from shape-independent to spheres section.
# Ought to be able to add polydisp to core? And add ability to x by S(Q) ?
