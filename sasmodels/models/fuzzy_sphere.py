r"""
For information about polarised and magnetic scattering, see
the :ref:`magnetism` documentation.

Definition
----------

The scattering intensity $I(q)$ is calculated as:

.. math::

    I(q) = \frac{\text{scale}}{V}(\Delta \rho)^2 A^2(q) S(q)
           + \text{background}


where the amplitude $A(q)$ is given as the typical sphere scattering convoluted
with a Gaussian to get a gradual drop-off in the scattering length density:

.. math::

    A(q) = \frac{3\left[\sin(qR) - qR \cos(qR)\right]}{(qR)^3}
           \exp\left(\frac{-(\sigma_\text{fuzzy}q)^2}{2}\right)

Here $A(q)^2$ is the form factor, $P(q)$. The scale is equivalent to the
volume fraction of spheres, each of volume, $V$. Contrast $(\Delta \rho)$
is the difference of scattering length densities of the sphere and the
surrounding solvent.

Poly-dispersion in radius and in fuzziness is provided for, though the
fuzziness must be kept much smaller than the sphere radius for meaningful
results.

From the reference:

  The "fuzziness" of the interface is defined by the parameter
  $\sigma_\text{fuzzy}$. The particle radius $R$ represents the radius of the
  particle where the scattering length density profile decreased to 1/2 of the
  core density. $\sigma_\text{fuzzy}$ is the width of the smeared particle
  surface; i.e., the standard deviation from the average height of the fuzzy
  interface. The inner regions of the microgel that display a higher density
  are described by the radial box profile extending to a radius of
  approximately $R_\text{box} \sim R - 2 \sigma$. The profile approaches
  zero as $R_\text{sans} \sim R + 2\sigma$.

For 2D data: The 2D scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math:: q = \sqrt{{q_x}^2 + {q_y}^2}

References
----------

#. M Stieger, J. S Pedersen, P Lindner, W Richtering,
   *Langmuir*, 20 (2004) 7283-7292

Authorship and Verification
----------------------------

* **Author:**
* **Last Modified by:**
* **Last Reviewed by:**
"""

import numpy as np
from numpy import inf

name = "fuzzy_sphere"
title = "Scattering from spherical particles with a fuzzy surface."
description = """\
scale: scale factor times volume fraction,
or just volume fraction for absolute scale data
radius: radius of the solid sphere
fuzziness = the standard deviation of the fuzzy interfacial
thickness (ie., so-called interfacial roughness)
sld: the SLD of the sphere
solvend_sld: the SLD of the solvent
background: incoherent background
Note: By definition, this function works only when fuzziness << radius.
"""
category = "shape:sphere"

# pylint: disable=bad-whitespace,line-too-long
# ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld",         "1e-6/Ang^2",  1, [-inf, inf], "sld",    "Particle scattering length density"],
              ["sld_solvent", "1e-6/Ang^2",  3, [-inf, inf], "sld",    "Solvent scattering length density"],
              ["radius",      "Ang",        60, [0, inf],    "volume", "Sphere radius"],
              ["fuzziness",   "Ang",        10, [0, inf],    "volume",       "std deviation of Gaussian convolution for interface (must be << radius)"],
             ]
# pylint: enable=bad-whitespace,line-too-long

source = ["lib/sas_3j1x_x.c", "fuzzy_sphere.c"]
has_shape_visualization = True

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for fuzzy sphere visualization."""
    import numpy as np
    radius = params.get('radius', 60)
    fuzziness = params.get('fuzziness', 10)

    phi = np.linspace(0, np.pi, resolution//2)
    theta = np.linspace(0, 2*np.pi, resolution)
    phi_mesh, theta_mesh = np.meshgrid(phi, theta)

    # Core sphere (where density is ~full)
    core_radius = radius - 2 * fuzziness
    core_radius = max(core_radius, radius * 0.5)  # Ensure visible core
    x_core = core_radius * np.sin(phi_mesh) * np.cos(theta_mesh)
    y_core = core_radius * np.sin(phi_mesh) * np.sin(theta_mesh)
    z_core = core_radius * np.cos(phi_mesh)

    # Outer sphere (where density approaches zero)
    outer_radius = radius + 2 * fuzziness
    x_outer = outer_radius * np.sin(phi_mesh) * np.cos(theta_mesh)
    y_outer = outer_radius * np.sin(phi_mesh) * np.sin(theta_mesh)
    z_outer = outer_radius * np.cos(phi_mesh)

    return {
        'core': (x_core, y_core, z_core),
        'fuzzy_boundary': (x_outer, y_outer, z_outer)
    }

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of the fuzzy sphere."""
    import numpy as np
    radius = params.get('radius', 60)
    fuzziness = params.get('fuzziness', 10)

    theta = np.linspace(0, 2*np.pi, 100)

    # Core, half-density, and outer radii
    core_radius = max(radius - 2 * fuzziness, radius * 0.5)
    outer_radius = radius + 2 * fuzziness

    core_x = core_radius * np.cos(theta)
    core_y = core_radius * np.sin(theta)
    mid_x = radius * np.cos(theta)
    mid_y = radius * np.sin(theta)
    outer_x = outer_radius * np.cos(theta)
    outer_y = outer_radius * np.sin(theta)

    # XY plane
    ax_xy.fill(outer_x, outer_y, 'lightyellow', alpha=0.3, label='Fuzzy interface')
    ax_xy.fill(mid_x, mid_y, 'lightblue', alpha=0.4, label='Half-density (R)')
    ax_xy.fill(core_x, core_y, 'blue', alpha=0.3, label='Dense core')
    ax_xy.plot(outer_x, outer_y, 'y--', linewidth=1.5, alpha=0.7)
    ax_xy.plot(mid_x, mid_y, 'b-', linewidth=2)
    ax_xy.plot(core_x, core_y, 'b:', linewidth=1.5)
    ax_xy.set_xlim(-outer_radius*1.2, outer_radius*1.2)
    ax_xy.set_ylim(-outer_radius*1.2, outer_radius*1.2)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title('XY Cross-section')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)
    ax_xy.legend(fontsize=8)

    # XZ plane
    ax_xz.fill(outer_x, outer_y, 'lightyellow', alpha=0.3)
    ax_xz.fill(mid_x, mid_y, 'lightcoral', alpha=0.4)
    ax_xz.fill(core_x, core_y, 'red', alpha=0.3)
    ax_xz.plot(outer_x, outer_y, 'y--', linewidth=1.5, alpha=0.7)
    ax_xz.plot(mid_x, mid_y, 'r-', linewidth=2)
    ax_xz.plot(core_x, core_y, 'r:', linewidth=1.5)
    ax_xz.set_xlim(-outer_radius*1.2, outer_radius*1.2)
    ax_xz.set_ylim(-outer_radius*1.2, outer_radius*1.2)
    ax_xz.set_xlabel('X (Å)')
    ax_xz.set_ylabel('Z (Å)')
    ax_xz.set_title('XZ Cross-section')
    ax_xz.set_aspect('equal')
    ax_xz.grid(True, alpha=0.3)

    # Add annotations
    ax_xz.annotate('', xy=(radius, 0), xytext=(0, 0),
                  arrowprops=dict(arrowstyle='->', color='black'))
    ax_xz.text(radius/2, -outer_radius*0.15, f'R={radius:.0f}Å', ha='center', fontsize=9)
    ax_xz.text(outer_radius*0.7, outer_radius*0.7, f'σ={fuzziness:.0f}Å', fontsize=9)

    # YZ plane
    ax_yz.fill(outer_x, outer_y, 'lightyellow', alpha=0.3)
    ax_yz.fill(mid_x, mid_y, 'lightgreen', alpha=0.4)
    ax_yz.fill(core_x, core_y, 'green', alpha=0.3)
    ax_yz.plot(outer_x, outer_y, 'y--', linewidth=1.5, alpha=0.7)
    ax_yz.plot(mid_x, mid_y, 'g-', linewidth=2)
    ax_yz.plot(core_x, core_y, 'g:', linewidth=1.5)
    ax_yz.set_xlim(-outer_radius*1.2, outer_radius*1.2)
    ax_yz.set_ylim(-outer_radius*1.2, outer_radius*1.2)
    ax_yz.set_xlabel('Y (Å)')
    ax_yz.set_ylabel('Z (Å)')
    ax_yz.set_title('YZ Cross-section')
    ax_yz.set_aspect('equal')
    ax_yz.grid(True, alpha=0.3)
have_Fq = True
radius_effective_modes = ["radius", "radius + fuzziness"]

def random():
    """Return a random parameter set for the model."""
    radius = 10**np.random.uniform(1, 4.7)
    fuzziness = 10**np.random.uniform(-2, -0.5)*radius  # 1% to 31% fuzziness
    pars = dict(
        radius=radius,
        fuzziness=fuzziness,
    )
    return pars

tests = [
    # Accuracy tests based on content in test/utest_models_new1_3.py
    #[{'background': 0.001}, 1.0, 0.001],

    [{}, 0.00301005, 359.2315],

    ]
