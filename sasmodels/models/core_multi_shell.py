r"""
Definition
----------

This model is a trivial extension of the CoreShell function to a larger number
of shells. The scattering length density profile for the default sld values
(w/ 4 shells).

.. figure:: img/core_multi_shell_sld_default_profile.jpg

    SLD profile of the core_multi_shell object from the center of sphere out
    for the default SLDs.*

The 2D scattering intensity is the same as $P(q)$ above, regardless of the
orientation of the $\vec q$ vector which is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}

.. note:: **Be careful!** The SLDs and scale can be highly correlated. Hold as
         many of these parameters fixed as possible.

.. note:: The outer most radius (= *radius* + *thickness*) is used as the
          effective radius for $S(Q)$ when $P(Q)*S(Q)$ is applied.

For information about polarised and magnetic scattering, see
the :ref:`magnetism` documentation.

Our model uses the form factor calculations implemented in a C-library provided
by the NIST Center for Neutron Research [#Kline2006]_.

References
----------

Also see the :ref:`core-shell-sphere` model documentation and [#Feigin1987]_

.. [#Kline2006] S R Kline, *J Appl. Cryst.*, 39 (2006) 895

.. [#Feigin1987] L A Feigin and D I Svergun, *Structure Analysis by
   Small-Angle X-Ray and Neutron Scattering*, Plenum Press, New York, 1987.

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Paul Kienzle **Date:** September 12, 2016
* **Last Reviewed by:** Paul Kienzle **Date:** September 12, 2016
"""

import numpy as np
from numpy import inf

name = "core_multi_shell"
title = "This model provides the scattering from a spherical core with 1 to 10 \
 concentric shell structures. The SLDs of the core and each shell are \
 individually specified."

description = """\
Form factor for a core muti-shell sphere normalized by the volume.
Each of up to 10 shells can have a unique thickness and sld.

	background:background,
	rad_core0: radius of sphere(core)
	thick_shell#:the thickness of the shell#
	sld_core0: the SLD of the sphere
	sld_solv: the SLD of the solvent
	sld_shell: the SLD of the shell#
	A_shell#: the coefficient in the exponential function


    scale: 1.0 if data is on absolute scale
    volfraction: volume fraction of spheres
    radius: the radius of the core
    sld: the SLD of the core
    thick_shelli: the thickness of the i'th shell from the core
    sld_shelli: the SLD of the i'th shell from the core
    sld_solvent: the SLD of the solvent
    background: incoherent background

"""

category = "shape:sphere"


#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld_core", "1e-6/Ang^2", 1.0, [-inf, inf], "sld",
               "Core scattering length density"],
              ["radius", "Ang", 200., [0, inf], "volume",
               "Radius of the core"],
              ["sld_solvent", "1e-6/Ang^2", 6.4, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["n", "", 1, [0, 10], "volume",
               "number of shells"],
              ["sld[n]", "1e-6/Ang^2", 1.7, [-inf, inf], "sld",
               "scattering length density of shell k"],
              ["thickness[n]", "Ang", 40., [0, inf], "volume",
               "Thickness of shell k"],
             ]

source = ["lib/sas_3j1x_x.c", "core_multi_shell.c"]
has_shape_visualization = True
have_Fq = True

def create_shape_mesh(params, resolution=50):
    """Create 3D mesh for core multi-shell sphere visualization."""
    import numpy as np
    radius = params.get('radius', 200)
    n = int(params.get('n', 1))

    phi = np.linspace(0, np.pi, resolution//2)
    theta = np.linspace(0, 2*np.pi, resolution)
    phi_mesh, theta_mesh = np.meshgrid(phi, theta)

    mesh_data = {}

    # Core sphere
    x_core = radius * np.sin(phi_mesh) * np.cos(theta_mesh)
    y_core = radius * np.sin(phi_mesh) * np.sin(theta_mesh)
    z_core = radius * np.cos(phi_mesh)
    mesh_data['core'] = (x_core, y_core, z_core)

    # Add shells
    current_radius = radius
    for i in range(n):
        thickness = params.get(f'thickness{i+1}', params.get('thickness', 40))
        if isinstance(thickness, (list, np.ndarray)):
            thickness = thickness[i] if i < len(thickness) else thickness[-1]
        current_radius += thickness

        x_shell = current_radius * np.sin(phi_mesh) * np.cos(theta_mesh)
        y_shell = current_radius * np.sin(phi_mesh) * np.sin(theta_mesh)
        z_shell = current_radius * np.cos(phi_mesh)
        mesh_data[f'shell_{i+1}'] = (x_shell, y_shell, z_shell)

    return mesh_data

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of core multi-shell sphere."""
    import numpy as np
    radius = params.get('radius', 200)
    n = int(params.get('n', 1))

    theta = np.linspace(0, 2*np.pi, 100)
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']

    # Calculate all radii
    radii = [radius]
    current_radius = radius
    for i in range(n):
        thickness = params.get(f'thickness{i+1}', params.get('thickness', 40))
        if isinstance(thickness, (list, np.ndarray)):
            thickness = thickness[i] if i < len(thickness) else thickness[-1]
        current_radius += thickness
        radii.append(current_radius)

    max_r = radii[-1] * 1.2

    # XY plane
    for i, r in enumerate(radii):
        circle_x = r * np.cos(theta)
        circle_y = r * np.sin(theta)
        color = colors[i % len(colors)]
        label = 'Core' if i == 0 else f'Shell {i}'
        ax_xy.plot(circle_x, circle_y, '-', color=color, linewidth=2, label=label)
        if i == 0:
            ax_xy.fill(circle_x, circle_y, color=color, alpha=0.2)

    ax_xy.set_xlim(-max_r, max_r)
    ax_xy.set_ylim(-max_r, max_r)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title('XY Cross-section')
    ax_xy.set_aspect('equal')
    ax_xy.legend(loc='upper right', fontsize=7)
    ax_xy.grid(True, alpha=0.3)

    # XZ plane
    for i, r in enumerate(radii):
        circle_x = r * np.cos(theta)
        circle_z = r * np.sin(theta)
        color = colors[i % len(colors)]
        ax_xz.plot(circle_x, circle_z, '-', color=color, linewidth=2)
        if i == 0:
            ax_xz.fill(circle_x, circle_z, color=color, alpha=0.2)

    ax_xz.set_xlim(-max_r, max_r)
    ax_xz.set_ylim(-max_r, max_r)
    ax_xz.set_xlabel('X (Å)')
    ax_xz.set_ylabel('Z (Å)')
    ax_xz.set_title('XZ Cross-section')
    ax_xz.set_aspect('equal')
    ax_xz.grid(True, alpha=0.3)

    # YZ plane
    for i, r in enumerate(radii):
        circle_y = r * np.cos(theta)
        circle_z = r * np.sin(theta)
        color = colors[i % len(colors)]
        ax_yz.plot(circle_y, circle_z, '-', color=color, linewidth=2)
        if i == 0:
            ax_yz.fill(circle_y, circle_z, color=color, alpha=0.2)

    ax_yz.set_xlim(-max_r, max_r)
    ax_yz.set_ylim(-max_r, max_r)
    ax_yz.set_xlabel('Y (Å)')
    ax_yz.set_ylabel('Z (Å)')
    ax_yz.set_title('YZ Cross-section')
    ax_yz.set_aspect('equal')
    ax_yz.grid(True, alpha=0.3)
radius_effective_modes = ["outer radius", "core radius"]

def random():
    """Return a random parameter set for the model."""
    num_shells = np.minimum(np.random.poisson(3)+1, 10)
    total_radius = 10**np.random.uniform(1.7, 4)
    thickness = np.random.exponential(size=num_shells+1)
    thickness *= total_radius/np.sum(thickness)
    pars = dict(
        #background=0,
        n=num_shells,
        radius=thickness[0],
    )
    for k, v in enumerate(thickness[1:]):
        pars['thickness%d'%(k+1)] = v
    return pars

def profile(sld_core, radius, sld_solvent, n, sld, thickness):
    """
    Returns the SLD profile *r* (Ang), and *rho* (1e-6/Ang^2).
    """
    n = int(n+0.5)
    z = []
    rho = []

    # add in the core
    z.append(0)
    rho.append(sld_core)
    z.append(radius)
    rho.append(sld_core)

    # add in the shells
    for k in range(int(n)):
        # Left side of each shells
        z.append(z[-1])
        rho.append(sld[k])
        z.append(z[-1] + thickness[k])
        rho.append(sld[k])
    # add in the solvent
    z.append(z[-1])
    rho.append(sld_solvent)
    z.append(z[-1]*1.25)
    rho.append(sld_solvent)

    return np.asarray(z), np.asarray(rho)
