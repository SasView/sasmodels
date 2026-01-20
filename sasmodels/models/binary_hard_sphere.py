r"""
Definition
----------

The binary hard sphere model provides the scattering intensity, for binary
mixture of hard spheres including hard sphere interaction between those
particles, using rhw Percus-Yevick closure. The calculation is an exact
multi-component solution that properly accounts for the 3 partial structure
factors as follows:

.. math::

    I(q) = (1-x)f_1^2(q) S_{11}(q) + 2[x(1-x)]^{1/2} f_1(q)f_2(q)S_{12}(q) +
    x\,f_2^2(q)S_{22}(q)

where $S_{ij}$ are the partial structure factors and $f_i$ are the scattering
amplitudes of the particles. The subscript 1 is for the smaller particle and 2
is for the larger. The number fraction of the larger particle,
($x = n2/(n1+n2)$, where $n$ = the number density) is internally calculated
based on the diameter ratio and the volume fractions.

.. math::
    :nowrap:

    \begin{align*}
    x &= \frac{(\phi_2 / \phi)\alpha^3}{(1-(\phi_2/\phi) + (\phi_2/\phi)
    \alpha^3)} \\
    \phi &= \phi_1 + \phi_2 = \text{total volume fraction} \\
    \alpha &= R_1/R_2 = \text{size ratio}
    \end{align*}

The 2D scattering intensity is the same as 1D, regardless of the orientation of
the *q* vector which is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


**NOTE 1:** The volume fractions and the scattering contrasts are loosely
correlated, so holding as many parameters fixed to known values during fitting
will improve the robustness of the fit.

**NOTE 2:** Since the calculation uses the Percus-Yevick closure, all of the
limitations of that closure relation apply here. Specifically, one should be
wary of results for (total) volume fractions greater than approximately 40%.
Depending on the size ratios or number fractions, the limit on total volume
fraction may be lower.

**NOTE 3:** The heavy arithmetic operations also mean that at present the
function is poorly behaved at very low qr.  In some cases very large qr may
also be poorly behaved.  These should however be outside any useful region of
qr.

The code for this model is based originally on a c-library implementation by the
NIST Center for Neutron Research (Kline, 2006).

See the references for details.

References
----------

#. N W Ashcroft and D C Langreth, *Physical Review*, 156 (1967) 685-692
   [Errata found in *Phys. Rev.* 166 (1968) 934]

#. S R Kline, *J Appl. Cryst.*, 39 (2006) 895

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Paul Butler **Date:** March 20, 2016
* **Last Reviewed by:** Paul Butler **Date:** March 20, 2016
"""

import numpy as np
from numpy import inf

category = "shape:sphere"
has_shape_visualization = True
single = False  # double precision only!

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for binary hard sphere visualization."""
    import numpy as np
    radius_lg = params.get('radius_lg', 100)
    radius_sm = params.get('radius_sm', 25)

    # Create spheres
    phi = np.linspace(0, np.pi, resolution//2)
    theta = np.linspace(0, 2*np.pi, resolution)
    phi_mesh, theta_mesh = np.meshgrid(phi, theta)

    # Large sphere (centered at origin)
    x_lg = radius_lg * np.sin(phi_mesh) * np.cos(theta_mesh)
    y_lg = radius_lg * np.sin(phi_mesh) * np.sin(theta_mesh)
    z_lg = radius_lg * np.cos(phi_mesh)

    # Small sphere (offset to the side)
    offset = radius_lg * 1.5
    x_sm = radius_sm * np.sin(phi_mesh) * np.cos(theta_mesh) + offset
    y_sm = radius_sm * np.sin(phi_mesh) * np.sin(theta_mesh)
    z_sm = radius_sm * np.cos(phi_mesh)

    return {
        'large_sphere': (x_lg, y_lg, z_lg),
        'small_sphere': (x_sm, y_sm, z_sm)
    }

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of binary hard spheres."""
    import numpy as np
    radius_lg = params.get('radius_lg', 100)
    radius_sm = params.get('radius_sm', 25)

    theta = np.linspace(0, 2*np.pi, 100)
    offset = radius_lg * 1.5

    # Large sphere circle
    lg_x = radius_lg * np.cos(theta)
    lg_y = radius_lg * np.sin(theta)

    # Small sphere circle (offset)
    sm_x = radius_sm * np.cos(theta) + offset
    sm_y = radius_sm * np.sin(theta)

    max_extent = offset + radius_sm * 1.5

    # XY plane
    ax_xy.plot(lg_x, lg_y, 'b-', linewidth=2, label=f'Large (R={radius_lg:.0f}Å)')
    ax_xy.fill(lg_x, lg_y, 'lightblue', alpha=0.3)
    ax_xy.plot(sm_x, sm_y, 'r-', linewidth=2, label=f'Small (R={radius_sm:.0f}Å)')
    ax_xy.fill(sm_x, sm_y, 'lightcoral', alpha=0.3)
    ax_xy.set_xlim(-radius_lg*1.5, max_extent)
    ax_xy.set_ylim(-radius_lg*1.5, radius_lg*1.5)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title('XY Cross-section')
    ax_xy.set_aspect('equal')
    ax_xy.legend(loc='upper left', fontsize=8)
    ax_xy.grid(True, alpha=0.3)

    # XZ plane
    ax_xz.plot(lg_x, lg_y, 'b-', linewidth=2)
    ax_xz.fill(lg_x, lg_y, 'lightblue', alpha=0.3)
    ax_xz.plot(sm_x, sm_y, 'r-', linewidth=2)
    ax_xz.fill(sm_x, sm_y, 'lightcoral', alpha=0.3)
    ax_xz.set_xlim(-radius_lg*1.5, max_extent)
    ax_xz.set_ylim(-radius_lg*1.5, radius_lg*1.5)
    ax_xz.set_xlabel('X (Å)')
    ax_xz.set_ylabel('Z (Å)')
    ax_xz.set_title('XZ Cross-section')
    ax_xz.set_aspect('equal')
    ax_xz.grid(True, alpha=0.3)

    # YZ plane
    ax_yz.plot(lg_y, lg_y, 'b-', linewidth=2)
    ax_yz.fill(lg_y, lg_y, 'lightblue', alpha=0.3)
    # Small sphere at origin in YZ plane
    ax_yz.plot(sm_y - offset, sm_y, 'r-', linewidth=2)
    ax_yz.fill(sm_y - offset, sm_y, 'lightcoral', alpha=0.3)
    ax_yz.set_xlim(-radius_lg*1.5, radius_lg*1.5)
    ax_yz.set_ylim(-radius_lg*1.5, radius_lg*1.5)
    ax_yz.set_xlabel('Y (Å)')
    ax_yz.set_ylabel('Z (Å)')
    ax_yz.set_title('YZ Cross-section')
    ax_yz.set_aspect('equal')
    ax_yz.grid(True, alpha=0.3)

name = "binary_hard_sphere"
title = "binary mixture of hard spheres with hard sphere interactions."
description = """Describes the scattering from a mixture of two distinct
monodisperse, hard sphere particles.
        [Parameters];
        radius_lg: large radius of binary hard sphere,
        radius_sm: small radius of binary hard sphere,
        volfraction_lg: volume fraction of large spheres,
        volfraction_sm: volume fraction of small spheres,
        sld_lg: large sphere  scattering length density,
        sld_sm: small sphere scattering length density,
        sld_solvent: solvent scattering length density.
"""
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["radius_lg", "Ang", 100, [0, inf], "",
               "radius of large particle"],
              ["radius_sm", "Ang", 25, [0, inf], "",
               "radius of small particle"],
              ["volfraction_lg", "", 0.1, [0, 1], "",
               "volume fraction of large particle"],
              ["volfraction_sm", "", 0.2, [0, 1], "",
               "volume fraction of small particle"],
              ["sld_lg", "1e-6/Ang^2", 3.5, [-inf, inf], "sld",
               "scattering length density of large particle"],
              ["sld_sm", "1e-6/Ang^2", 0.5, [-inf, inf], "sld",
               "scattering length density of small particle"],
              ["sld_solvent", "1e-6/Ang^2", 6.36, [-inf, inf], "sld",
               "Solvent scattering length density"],
             ]

source = ["lib/sas_3j1x_x.c", "binary_hard_sphere.c"]

def random():
    """Return a random parameter set for the model."""
    # TODO: binary_hard_sphere fails at low qr
    radius_lg = 10**np.random.uniform(2, 4.7)
    radius_sm = 10**np.random.uniform(2, 4.7)
    volfraction_lg = 10**np.random.uniform(-3, -0.3)
    volfraction_sm = 10**np.random.uniform(-3, -0.3)
    # TODO: Get slightly different results if large and small are swapped
    # modify the model so it doesn't care which is which
    if radius_lg < radius_sm:
        radius_lg, radius_sm = radius_sm, radius_lg
        volfraction_lg, volfraction_sm = volfraction_sm, volfraction_lg
    pars = dict(
        radius_lg=radius_lg,
        radius_sm=radius_sm,
        volfraction_lg=volfraction_lg,
        volfraction_sm=volfraction_sm,
    )
    return pars

# NOTE: test results taken from values returned by SasView 3.1.2
tests = [[{}, 0.001, 25.8927262013]]
