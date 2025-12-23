# superball model
# Note: model title and parameter table are inserted automatically
r"""
Definition
----------

.. figure:: img/superball_realSpace.png

   Superball visualisation for varied values of the parameter p.

This model calculates the scattering of a superball, which represents a cube
with rounded edges. It can be used to describe nanoparticles that deviate from
the perfect cube shape as it is often observed experimentally
[#WetterskogSuperball]_. The shape is described by

.. math::

    x^{2p} + y^{2p} + z^{2p} \leq \biggl( \frac{a}{2} \biggr)^{2p}

with $a$ the cube edge length of the superball and $p$ a parameter that
describes the roundness of the edges. In the limiting cases $p=1$ the superball
corresponds to a sphere with radius $R = a/2$ and for $p = \infty$ to a cube
with edge length $a$. The exponent $p$ is related to $a$ and the face diagonal
$d$ via

.. math::
    p = \frac{1}{1 + 2 \mathrm{log}_2 (a/d)}.

.. figure:: img/superball_geometry2d.png

    Cross-sectional view of a superball showing the principal axis length $a$,
    the face-diagonal $d$ and the superball radius $R$.

The oriented form factor is determined by solving

.. math::
    p_o(\vec{q}) =& \int_{V} \mathrm{d} \vec{r} e^{i \vec{q} \cdot \vec{r}}\\
        =& \frac{a^3}{8} \int_{-1}^{1} \mathrm{d} x \int_{-\gamma}^{\gamma}
            \mathrm{d} y \int_{-\zeta}^{\zeta} \mathrm{d} z
            e^{i a (q_x x + q_y y + q_z z) / 2}\\
        =& \frac{a^2}{2 q_z} \int_{-1}^{1} \mathrm{d} x \int_{-\gamma}^{\gamma}
            \mathrm{d} y  e^{i a(q_x x + q_y y)/2}
            \sin(q_z a \zeta / 2),

with

.. math::
    \gamma =& \sqrt[2p]{1-x^{2p}}, \\
    \zeta =& \sqrt[2p]{1-x^{2p} -y^{2p}}.

The integral can be transformed to

.. math::
    p_o(\vec{q}) = \frac{2 a^2}{q_z} \int_{0}^{1} \mathrm{d} x \, \cos
        \biggl(\frac{a q_x x}{2} \biggr) \int_{0}^{\gamma} \mathrm{d} y \,
        \cos \biggl( \frac{a q_y y}{2} \biggr) \sin
        \biggl( \frac{a q_z \zeta}{2} \biggr),

which can be solved numerically.

The orientational average is then obtained by calculating

.. math::
    P(q) = \int_0^{\tfrac{\pi}{2}} \mathrm{d} \varphi \int_0^{\tfrac{\pi}{2}}
        \mathrm{d} \theta \, \sin (\theta) | p_o(\vec{q}) |^2

with

.. math::
    \vec{q} &= q \begin{pmatrix} \cos (\varphi) \sin (\theta)\\
    \sin (\varphi) \sin(\theta)\\
    \cos (\theta)\end{pmatrix}

The implemented orientationally averaged superball model is then fully given by
[#DresenSuperball]_

.. math::
    I(q) = \mathrm{scale} (\Delta \rho)^2 P(q) + \mathrm{background}.


FITTING NOTES
~~~~~~~~~~~~~

Validation
----------

    The code is validated by reproducing the spherical form factor implemented
    in SasView for $p = 1$ and the parallelepiped form factor with $a = b = c$ for
    $p = 1000$. The form factors match in the first order oscillation with a
    precision in the order of $10^{-4}$. The agreement worsens for higher order
    oscillations and beyond the third order oscillations a higher order Gauss
    quadrature rule needs to be used to keep the agreement below $10^{-3}$.
    This is however avoided in this implementation to keep the computation time
    fast.

References
----------

.. [#WetterskogSuperball] E. Wetterskog, A. Klapper, S. Disch, E. Josten, R. P. Hermann, U. Rücker, T. Brückel, L. Bergström and G. Salazar-Alvarez, *Nanoscale*, 8 (2016) 15571

.. [#DresenSuperball] D. Dresen, A. Qdemat, S. Ulusoy, F. Mees, D. Zakutna, E. Wetterskog, E. Kentzinger, G. Salazar-Alvarez, S. Disch, *J. Phys. Chem. C* (2021), doi: 10.1021/acs.jpcc.1c06082

Source
------

`superball.py <https://github.com/SasView/sasmodels/blob/master/sasmodels/models/superball.py>`_

`superball.c <https://github.com/SasView/sasmodels/blob/master/sasmodels/models/superball.c>`_

Authorship and Verification
----------------------------

* **Author:** Dominique Dresen **Date:** March 27, 2019
* **Last Modified by:** Dominique Dresen **Date:** March 27, 2019
* **Last Reviewed by:** Dirk Honecker **Date:** November 05, 2021
* **Source added by :** Dominique Dresen **Date:** March 27, 2019"""

import numpy as np
from numpy import inf

# saved in utf-8 encoding for the German umlaut (üö)

name = "superball"
title = "Superball with uniform scattering length density."

description = """
    I(q)= scale*V*(sld - sld_solvent)^2*P(q)+background
        P(q) = (2/pi) * double integral from 0 to pi/2 of ...
           AP^2(q)*sin(theta)*dtheta*dphi
        AP = integral from -1 to 1 integral from -g to g of ...
        cos(R qx x) cos(R qy y]) sin(R qz Ze)
        g = (1 - x^(2p))^(1/(2p))
        Ze = (1 - x^(2p) - y^(2p))^(1/(2p))
"""
category = "shape:sphere"
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 4, [-inf, inf], "sld",
               "Superball scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["length_a", "Ang", 50, [0, inf], "volume",
               "Cube edge length of the superball"],
              ["exponent_p", "", 2.5, [0, inf], "volume",
               "Exponent describing the roundness of the superball"],
              ["theta", "degrees", 0, [-360, 360], "orientation",
               "c axis to beam angle"],
              ["phi", "degrees", 0, [-360, 360], "orientation",
               "rotation about beam"],
              ["psi", "degrees", 0, [-360, 360], "orientation",
               "rotation about c axis"],
              ]
# lib/gauss76.c
# lib/gauss20.c
source = ["lib/gauss20.c", "lib/sas_gamma.c", "superball.c"]
has_shape_visualization = True

def create_shape_mesh(params, resolution=60):
    """Create 3D mesh for superball visualization."""
    import numpy as np
    length_a = params.get('length_a', 50)
    p = params.get('exponent_p', 2.5)

    # Superball equation: |x|^(2p) + |y|^(2p) + |z|^(2p) <= (a/2)^(2p)
    # Parametric surface approximation
    R = length_a / 2

    # Use spherical-like coordinates but with superellipse shape
    u = np.linspace(0, np.pi, resolution)
    v = np.linspace(0, 2*np.pi, resolution)
    u_mesh, v_mesh = np.meshgrid(u, v)

    # Superellipsoid parametric equations
    def sign_pow(x, n):
        return np.sign(x) * np.abs(x)**n

    cos_u = np.cos(u_mesh)
    sin_u = np.sin(u_mesh)
    cos_v = np.cos(v_mesh)
    sin_v = np.sin(v_mesh)

    # Exponent for shape (1/p gives the superellipsoid exponent)
    e = 1.0 / p

    x = R * sign_pow(sin_u, e) * sign_pow(cos_v, e)
    y = R * sign_pow(sin_u, e) * sign_pow(sin_v, e)
    z = R * sign_pow(cos_u, e)

    return {'superball': (x, y, z)}

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    """Plot 2D cross-sections of the superball."""
    import numpy as np
    length_a = params.get('length_a', 50)
    p = params.get('exponent_p', 2.5)
    R = length_a / 2

    # Superellipse: |x/R|^(2p) + |y/R|^(2p) = 1
    # Parametric: x = R*sign(cos(t))*|cos(t)|^(1/p), y = R*sign(sin(t))*|sin(t)|^(1/p)
    t = np.linspace(0, 2*np.pi, 200)

    def sign_pow(x, n):
        return np.sign(x) * np.abs(x)**n

    e = 1.0 / p
    superellipse_x = R * sign_pow(np.cos(t), e)
    superellipse_y = R * sign_pow(np.sin(t), e)

    # Reference shapes
    circle_x = R * np.cos(t)
    circle_y = R * np.sin(t)

    # XY plane
    ax_xy.fill(superellipse_x, superellipse_y, 'lightblue', alpha=0.5)
    ax_xy.plot(superellipse_x, superellipse_y, 'b-', linewidth=2, label=f'p={p:.1f}')
    ax_xy.plot(circle_x, circle_y, 'g--', linewidth=1, alpha=0.5, label='sphere (p=1)')
    # Square reference
    sq = R * np.array([-1, -1, 1, 1, -1])
    ax_xy.plot(sq, np.array([-1, 1, 1, -1, -1])*R, 'r--', linewidth=1, alpha=0.5, label='cube (p=∞)')

    ax_xy.set_xlim(-R*1.4, R*1.4)
    ax_xy.set_ylim(-R*1.4, R*1.4)
    ax_xy.set_xlabel('X (Å)')
    ax_xy.set_ylabel('Y (Å)')
    ax_xy.set_title(f'XY Cross-section (a={length_a:.0f}Å)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)
    ax_xy.legend(fontsize=8)

    # XZ plane
    ax_xz.fill(superellipse_x, superellipse_y, 'lightcoral', alpha=0.5)
    ax_xz.plot(superellipse_x, superellipse_y, 'r-', linewidth=2)
    ax_xz.plot(circle_x, circle_y, 'g--', linewidth=1, alpha=0.5)
    ax_xz.plot(sq, np.array([-1, 1, 1, -1, -1])*R, 'b--', linewidth=1, alpha=0.5)

    ax_xz.set_xlim(-R*1.4, R*1.4)
    ax_xz.set_ylim(-R*1.4, R*1.4)
    ax_xz.set_xlabel('X (Å)')
    ax_xz.set_ylabel('Z (Å)')
    ax_xz.set_title(f'XZ Cross-section (p={p:.2f})')
    ax_xz.set_aspect('equal')
    ax_xz.grid(True, alpha=0.3)

    # YZ plane
    ax_yz.fill(superellipse_x, superellipse_y, 'lightgreen', alpha=0.5)
    ax_yz.plot(superellipse_x, superellipse_y, 'g-', linewidth=2)
    ax_yz.set_xlim(-R*1.4, R*1.4)
    ax_yz.set_ylim(-R*1.4, R*1.4)
    ax_yz.set_xlabel('Y (Å)')
    ax_yz.set_ylabel('Z (Å)')
    ax_yz.set_title('YZ Cross-section')
    ax_yz.set_aspect('equal')
    ax_yz.grid(True, alpha=0.3)
have_Fq = True
radius_effective_modes = [
    "radius of gyration",
    "equivalent volume sphere",
    "half length_a",
]

def random():
    """Return a random parameter set for the model."""
    length = np.random.uniform(10, 500)
    exponent = np.random.uniform(1.5, 5)
    pars = dict(
        length_a=length,
        exponent_p=exponent)
    return pars

# parameters for demo
demo = dict(scale=1, background=0,
            sld=6.3, sld_solvent=1.0,
            length_a=100, exponent_p=2.5,
            theta=45, phi=30, psi=15,
            length_a_pd=0.1, length_a_pd_n=10,
            theta_pd=10, theta_pd_n=1,
            phi_pd=10, phi_pd_n=1,
            psi_pd=10, psi_pd_n=1)

tests = [
    [{}, 0.2, 0.76833],
    [{"length_a": 100., "exponent_p": 1, "sld": 6., "sld_solvent": 1.},
     0.2, 0.7263],
    [{"length_a": 100., "exponent_p": 1000, "sld": 6., "sld_solvent": 1.},
     0.2, 0.2714],
    [{"length_a": 100., "exponent_p": 2.5, "sld": 6., "sld_solvent": 1.},
     0.2, 0.2810],
    [{"length_a": 100., "exponent_p": 2.5, "sld": 6., "sld_solvent": 1.,
      "length_a_pd": 0.1, "length_a_pd_n": 10},
     0.2, 0.49551865],
]
