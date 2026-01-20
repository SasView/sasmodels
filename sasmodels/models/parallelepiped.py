# parallelepiped model
# Note: model title and parameter table are inserted automatically
r"""
Definition
----------

This model calculates the scattering from a rectangular solid
(:numref:`parallelepiped-image`).
If you need to apply polydispersity, see also :ref:`rectangular-prism`. For
information about polarised and magnetic scattering, see
the :ref:`magnetism` documentation.

.. _parallelepiped-image:


.. figure:: img/parallelepiped_geometry.jpg

   Parallelepiped with the corresponding definition of sides.

The three dimensions of the parallelepiped (strictly here a cuboid) may be
given in *any* size order as long as the particles are randomly oriented (i.e.
take on all possible orientations see notes on 2D below). To avoid multiple fit
solutions, especially with Monte-Carlo fit methods, it may be advisable to
restrict their ranges. There may be a number of closely similar "best fits", so
some trial and error, or fixing of some dimensions at expected values, may
help.

The form factor is normalized by the particle volume and the 1D scattering
intensity $I(q)$ is then calculated as:

.. Comment by Miguel Gonzalez:
   I am modifying the original text because I find the notation a little bit
   confusing. I think that in most textbooks/papers, the notation P(Q) is
   used for the form factor (adim, P(Q=0)=1), although F(q) seems also to
   be used. But here (as for many other models), P(q) is used to represent
   the scattering intensity (in cm-1 normally). It would be good to agree on
   a common notation.

.. math::

    I(q) = \frac{\text{scale}}{V} (\Delta\rho \cdot V)^2
           \left< P(q, \alpha, \beta) \right> + \text{background}

where the volume $V = A B C$, the contrast is defined as
$\Delta\rho = \rho_\text{p} - \rho_\text{solvent}$, $P(q, \alpha, \beta)$
is the form factor corresponding to a parallelepiped oriented
at an angle $\alpha$ (angle between the long axis C and $\vec q$), and $\beta$
(the angle between the projection of the particle in the $xy$ detector plane
and the $y$ axis) and the averaging $\left<\ldots\right>$ is applied over all
orientations.

Assuming $a = A/B < 1$, $b = B /B = 1$, and $c = C/B > 1$, the
form factor is given by (Mittelbach and Porod, 1961 [#Mittelbach1961]_)

.. math::

    P(q, \alpha) = \int_0^1 \phi_Q\left(\mu \sqrt{1-\sigma^2},a\right)
        \left[S(\mu c \sigma/2)\right]^2 d\sigma

with

.. math::

    \phi_Q(\mu,a) &= \int_0^1
        \left\{S\left[\frac{\mu}{2}\cos\left(\frac{\pi}{2}u\right)\right]
               S\left[\frac{\mu a}{2}\sin\left(\frac{\pi}{2}u\right)\right]
               \right\}^2 du \\
    S(x) &= \frac{\sin x}{x} \\
    \mu &= qB

where substitution of $\sigma = cos\alpha$ and $\beta = \pi/2 \ u$ have been
applied.

For **oriented** particles, the 2D scattering intensity, $I(q_x, q_y)$, is
given as:

.. math::

    I(q_x, q_y) = \frac{\text{scale}}{V} (\Delta\rho \cdot V)^2 P(q_x, q_y)
            + \text{background}

.. Comment by Miguel Gonzalez:
   This reflects the logic of the code, as in parallelepiped.c the call
   to _pkernel returns $P(q_x, q_y)$ and then this is multiplied by
   $V^2 * (\Delta \rho)^2$. And finally outside parallelepiped.c it will be
   multiplied by scale, normalized by $V$ and the background added. But
   mathematically it makes more sense to write
   $I(q_x, q_y) = \text{scale} V \Delta\rho^2 P(q_x, q_y) + \text{background}$,
   with scale being the volume fraction.

Where $P(q_x, q_y)$ for a given orientation of the form factor is calculated as

.. math::

    P(q_x, q_y) = \left[\frac{\sin(\tfrac{1}{2}qA\cos\alpha)}{(\tfrac{1}
                   {2}qA\cos\alpha)}\right]^2
                  \left[\frac{\sin(\tfrac{1}{2}qB\cos\beta)}{(\tfrac{1}
                   {2}qB\cos\beta)}\right]^2
                  \left[\frac{\sin(\tfrac{1}{2}qC\cos\gamma)}{(\tfrac{1}
                   {2}qC\cos\gamma)}\right]^2

with

.. math::

    \cos\alpha &= \hat A \cdot \hat q, \\
    \cos\beta  &= \hat B \cdot \hat q, \\
    \cos\gamma &= \hat C \cdot \hat q


FITTING NOTES
~~~~~~~~~~~~~

#. The 2nd virial coefficient of the parallelepiped is calculated based on
   the averaged effective radius, after appropriately sorting the three
   dimensions, to give an oblate or prolate particle, $(=\sqrt{AB/\pi})$ and
   length $(= C)$ values, and used as the effective radius for
   $S(q)$ when $P(q) \cdot S(q)$ is applied.

#. For 2d data the orientation of the particle is required, described using
   angles $\theta$, $\phi$ and $\Psi$ as in the diagrams below, where $\theta$
   and $\phi$ define the orientation of the director in the laboratory reference
   frame of the beam direction ($z$) and detector plane ($x-y$ plane), while
   the angle $\Psi$ is effectively the rotational angle around the particle
   $C$ axis. For $\theta = 0$ and $\phi = 0$, $\Psi = 0$ corresponds to the
   $B$ axis oriented parallel to the y-axis of the detector with $A$ along
   the x-axis. For other $\theta$, $\phi$ values, the order of rotations
   matters. In particular, the parallelepiped must first be rotated $\theta$
   degrees in the $x-z$ plane before rotating $\phi$ degrees around the $z$
   axis (in the $x-y$ plane). Applying orientational distribution to the
   particle orientation (i.e  `jitter` to one or more of these angles) can get
   more confusing as `jitter` is defined **NOT** with respect to the laboratory
   frame but the particle reference frame. It is thus highly recommended to
   read :ref:`orientation` for further details of the calculation and angular
   dispersions.

.. note:: For 2d, constraints must be applied during fitting to ensure that the
   order of sides chosen is not altered, and hence that the correct definition
   of angles is preserved. For the default choice shown here, that means
   ensuring that the inequality $A < B < C$ is not violated,  The calculation
   will not report an error, but the results may be not correct.

.. _parallelepiped-orientation:

.. figure:: img/parallelepiped_angle_definition.png

    Definition of the angles for oriented parallelepiped, shown with $A<B<C$.

.. figure:: img/parallelepiped_angle_projection.png

    Examples of the angles for an oriented parallelepiped against the
    detector plane.

.. Comment by Paul Butler
   I am commenting this section out as we are trying to minimize the amount of
   oritentational detail here and encourage the user to go to the full
   orientation documentation so that changes can be made in just one place.
   below is the commented paragraph:
   On introducing "Orientational Distribution" in the angles, "distribution of
   theta" and "distribution of phi" parameters will appear. These are actually
   rotations about axes $\delta_1$ and $\delta_2$ of the parallelepiped,
   perpendicular to the $a$ x $c$ and $b$ x $c$ faces. (When $\theta = \phi = 0$
   these are parallel to the $Y$ and $X$ axes of the instrument.) The third
   orientation distribution, in $\psi$, is about the $c$ axis of the particle,
   perpendicular to the $a$ x $b$ face. Some experimentation may be required to
   understand the 2d patterns fully as discussed in :ref:`orientation` .


Validation
----------

Validation of the code was done by comparing the output of the 1D calculation
to the angular average of the output of a 2D calculation over all possible
angles.

References
----------

See also Nayuk [#Nayuk2012]_ and Onsager [#Onsager1949]_.

.. [#Mittelbach1961] P Mittelbach and G Porod, *Acta Physica Austriaca*,
   14 (1961) 185-211
.. [#Nayuk2012] R Nayuk and K Huber, *Z. Phys. Chem.*, 226 (2012) 837-854
.. [#Onsager1949] L. Onsager, *Ann. New York Acad. Sci.*, 51 (1949) 627-659

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:**  Paul Kienzle **Date:** April 05, 2017
* **Last Reviewed by:**  Miguel Gonzales and Paul Butler **Date:** May 24,
  2018 - documentation updated
"""

import numpy as np
from numpy import inf

name = "parallelepiped"
title = "Rectangular parallelepiped with uniform scattering length density."
description = """
    I(q)= scale*V*(sld - sld_solvent)^2*P(q,alpha)+background
        P(q,alpha) = integral from 0 to 1 of ...
           phi(mu*sqrt(1-sigma^2),a) * S(mu*c*sigma/2)^2 * dsigma
        with
            phi(mu,a) = integral from 0 to 1 of ..
            (S((mu/2)*cos(pi*u/2))*S((mu*a/2)*sin(pi*u/2)))^2 * du
            S(x) = sin(x)/x
            mu = q*B
        V: Volume of the rectangular parallelepiped
        alpha: angle between the long axis of the
            parallelepiped and the q-vector for 1D
"""
category = "shape:parallelepiped"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 4, [-inf, inf], "sld",
               "Parallelepiped scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["length_a", "Ang", 35, [0, inf], "volume",
               "Shorter side of the parallelepiped"],
              ["length_b", "Ang", 75, [0, inf], "volume",
               "Second side of the parallelepiped"],
              ["length_c", "Ang", 400, [0, inf], "volume",
               "Larger side of the parallelepiped"],
              ["theta", "degrees", 60, [-360, 360], "orientation",
               "c axis to beam angle"],
              ["phi", "degrees", 60, [-360, 360], "orientation",
               "rotation about beam"],
              ["psi", "degrees", 60, [-360, 360], "orientation",
               "rotation about c axis"],
             ]

source = ["lib/gauss76.c", "parallelepiped.c"]
have_Fq = True
radius_effective_modes = [
    "equivalent cylinder excluded volume", "equivalent volume sphere",
    "half length_a", "half length_b", "half length_c",
    "equivalent circular cross-section", "half ab diagonal", "half diagonal",
    ]
has_shape_visualization = True

def create_shape_mesh(params, resolution=50):
    import numpy as np
    length_a = params.get('length_a', 35)
    length_b = params.get('length_b', 75)
    length_c = params.get('length_c', 400)

    # Create faces of the parallelepiped
    faces = {}

    # Front and back faces (x = ±length_a/2)
    y = np.linspace(-length_b/2, length_b/2, resolution//4)
    z = np.linspace(-length_c/2, length_c/2, resolution//2)
    y_mesh, z_mesh = np.meshgrid(y, z)

    faces['front'] = (np.full_like(y_mesh, length_a/2), y_mesh, z_mesh)
    faces['back'] = (np.full_like(y_mesh, -length_a/2), y_mesh, z_mesh)

    # Left and right faces (y = ±length_b/2)
    x = np.linspace(-length_a/2, length_a/2, resolution//4)
    z = np.linspace(-length_c/2, length_c/2, resolution//2)
    x_mesh, z_mesh = np.meshgrid(x, z)

    faces['right'] = (x_mesh, np.full_like(x_mesh, length_b/2), z_mesh)
    faces['left'] = (x_mesh, np.full_like(x_mesh, -length_b/2), z_mesh)

    # Top and bottom faces (z = ±length_c/2)
    x = np.linspace(-length_a/2, length_a/2, resolution//4)
    y = np.linspace(-length_b/2, length_b/2, resolution//4)
    x_mesh, y_mesh = np.meshgrid(x, y)

    faces['top'] = (x_mesh, y_mesh, np.full_like(x_mesh, length_c/2))
    faces['bottom'] = (x_mesh, y_mesh, np.full_like(x_mesh, -length_c/2))

    return faces

def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):
    length_a = params.get('length_a', 35)
    length_b = params.get('length_b', 75)
    length_c = params.get('length_c', 400)

    # XY plane (top view) - rectangle A x B
    rect_xy_x = [-length_a/2, -length_a/2, length_a/2, length_a/2, -length_a/2]
    rect_xy_y = [-length_b/2, length_b/2, length_b/2, -length_b/2, -length_b/2]

    ax_xy.plot(rect_xy_x, rect_xy_y, 'b-', linewidth=2, label='A × B face')
    ax_xy.fill(rect_xy_x, rect_xy_y, 'lightblue', alpha=0.3)
    ax_xy.set_xlim(-length_a/2*1.3, length_a/2*1.3)
    ax_xy.set_ylim(-length_b/2*1.3, length_b/2*1.3)
    ax_xy.set_xlabel('X (Å) - Length A')
    ax_xy.set_ylabel('Y (Å) - Length B')
    ax_xy.set_title('XY Cross-section (Top View)')
    ax_xy.set_aspect('equal')
    ax_xy.grid(True, alpha=0.3)

    # XZ plane (side view) - rectangle A x C
    rect_xz_x = [-length_a/2, -length_a/2, length_a/2, length_a/2, -length_a/2]
    rect_xz_z = [-length_c/2, length_c/2, length_c/2, -length_c/2, -length_c/2]

    ax_xz.plot(rect_xz_x, rect_xz_z, 'r-', linewidth=2, label='A × C face')
    ax_xz.fill(rect_xz_x, rect_xz_z, 'lightcoral', alpha=0.3)
    ax_xz.set_xlim(-length_a/2*1.3, length_a/2*1.3)
    ax_xz.set_ylim(-length_c/2*1.3, length_c/2*1.3)
    ax_xz.set_xlabel('X (Å) - Length A')
    ax_xz.set_ylabel('Z (Å) - Length C')
    ax_xz.set_title('XZ Cross-section (Side View)')
    ax_xz.grid(True, alpha=0.3)

    # YZ plane (front view) - rectangle B x C
    rect_yz_y = [-length_b/2, -length_b/2, length_b/2, length_b/2, -length_b/2]
    rect_yz_z = [-length_c/2, length_c/2, length_c/2, -length_c/2, -length_c/2]

    ax_yz.plot(rect_yz_y, rect_yz_z, 'g-', linewidth=2, label='B × C face')
    ax_yz.fill(rect_yz_y, rect_yz_z, 'lightgreen', alpha=0.3)
    ax_yz.set_xlim(-length_b/2*1.3, length_b/2*1.3)
    ax_yz.set_ylim(-length_c/2*1.3, length_c/2*1.3)
    ax_yz.set_xlabel('Y (Å) - Length B')
    ax_yz.set_ylabel('Z (Å) - Length C')
    ax_yz.set_title('YZ Cross-section (Front View)')
    ax_yz.grid(True, alpha=0.3)

    # Add dimension annotations
    # XY plane
    ax_xy.annotate('', xy=(-length_a/2, -length_b/2*1.4), xytext=(length_a/2, -length_b/2*1.4),
                  arrowprops=dict(arrowstyle='<->', color='black'))
    ax_xy.text(0, -length_b/2*1.5, f'A = {length_a:.0f} Å', ha='center', fontsize=10)

    ax_xy.annotate('', xy=(-length_a/2*1.4, -length_b/2), xytext=(-length_a/2*1.4, length_b/2),
                  arrowprops=dict(arrowstyle='<->', color='black'))
    ax_xy.text(-length_a/2*1.5, 0, f'B = {length_b:.0f} Å', ha='center', fontsize=10, rotation=90)

    # XZ plane
    ax_xz.annotate('', xy=(-length_a/2, -length_c/2*1.2), xytext=(length_a/2, -length_c/2*1.2),
                  arrowprops=dict(arrowstyle='<->', color='black'))
    ax_xz.text(0, -length_c/2*1.25, f'A = {length_a:.0f} Å', ha='center', fontsize=10)

    ax_xz.annotate('', xy=(-length_a/2*1.2, -length_c/2), xytext=(-length_a/2*1.2, length_c/2),
                  arrowprops=dict(arrowstyle='<->', color='black'))
    ax_xz.text(-length_a/2*1.25, 0, f'C = {length_c:.0f} Å', ha='center', fontsize=10, rotation=90)

def random():
    """Return a random parameter set for the model."""
    length = 10**np.random.uniform(1, 4.7, size=3)
    pars = dict(
        length_a=length[0],
        length_b=length[1],
        length_c=length[2],
    )
    return pars

# rkh 7/4/17 add random unit test for 2d, note make all params different,
# 2d values not tested against other codes or models
qx, qy = 0.2 * np.cos(np.pi/6.), 0.2 * np.sin(np.pi/6.)
tests = [[{}, 0.2, 0.17758004974],
         [{}, [0.2], [0.17758004974]],
         [{'theta':10.0, 'phi':20.0}, (qx, qy), 0.0089517140475],
         [{'theta':10.0, 'phi':20.0}, [(qx, qy)], [0.0089517140475]],
        ]
del qx, qy  # not necessary to delete, but cleaner
