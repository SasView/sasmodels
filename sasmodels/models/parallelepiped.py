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

Let $q_a$, $q_b$, and $q_c$ be the components of $\vec q$ along the particle
axes. The normalized form factor for a single orientation is the product of
three squared sinc factors (see [#Mittelbach1961]_ and [#Mao2017]_).

.. math::

    P(q_a,q_b,q_c) =
    S^2(Aq_a/2)\,S^2(Bq_b/2)\,S^2(Cq_c/2),
    \qquad S(x) = \frac{\sin x}{x}, \quad S(0) = 1.

Here $A$, $B$, and $C$ are full edge lengths. The factors of one half come
from the Fourier transform of a box centered at the origin.

For randomly oriented particles, the 1D intensity uses a normalized average
over the direction of $\vec q$ in the particle frame. Let $\vartheta$ be the
angle from the $C$ axis and $\varphi$ the azimuth from the $A$ axis.

.. math::

    (q_a,q_b,q_c) =
    q(\sin\vartheta\cos\varphi,\sin\vartheta\sin\varphi,\cos\vartheta).

The squared sinc factors are even, so all eight octants contribute equally.

.. math::

    \langle P \rangle(q) &=
    \frac{2}{\pi} \int_0^{\pi/2}\int_0^{\pi/2}
    P(q_a,q_b,q_c)\sin\vartheta\,d\varphi\,d\vartheta, \\
    I(q) &= \frac{\text{scale}}{V}(\Delta\rho V)^2
             \langle P \rangle(q) + \text{background}.

The volume is $V=ABC$ and the contrast is
$\Delta\rho=\rho_\text{p}-\rho_\text{solvent}$. The average is one at $q=0$.
The implementation integrates over $\sigma=\cos\vartheta$ and a rescaled
azimuth in $[0,1]$. Its substitution $\mu=qB/2$ includes the factor of one
half in each sinc argument. For example, the $C$ term uses
$\mu(C/B)\sigma=qC\sigma/2$.

For **oriented** particles, the 2D scattering intensity, $I(q_x, q_y)$, is
given as:

.. math::

    I(q_x, q_y) = \frac{\text{scale}}{V} (\Delta\rho \cdot V)^2 P(q_a,q_b,q_c)
            + \text{background}

The detector-plane vector $(q_x,q_y,0)$ is rotated into the particle frame
using the orientation angles below to obtain $(q_a,q_b,q_c)$.


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
.. [#Mao2017] Y Mao et al., *J. Phys. Chem. B*, 121 (2017) 1340-1351,
   Eq. (6), https://doi.org/10.1021/acs.jpcb.6b11425
.. [#Nayuk2012] R Nayuk and K Huber, *Z. Phys. Chem.*, 226 (2012) 837-854
.. [#Onsager1949] L Onsager, *Ann. New York Acad. Sci.*, 51 (1949) 627-659

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
description = r"""
    I(q) = scale*V*(sld - sld_solvent)^2*<P>(q) + background
        P(qa,qb,qc) = [sinc(A*qa/2)*sinc(B*qb/2)*sinc(C*qc/2)]^2
        sinc(x) = sin(x)/x, with sinc(0) = 1
        qa, qb, qc are the components of q in the particle frame.
        V = A*B*C, where A, B, C are the full edge lengths.
        <P>(q) averages P over directions of q for the 1D model.
        The oriented 2D model evaluates P for a single direction.
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
              ["theta", "degrees", 60, [-180, 180], "orientation",
               "c axis to beam angle"],
              ["phi", "degrees", 60, [-360, 360], "orientation",
               "rotation about beam"],
              ["psi", "degrees", 60, [-360, 360], "orientation",
               "rotation about c axis"],
             ]

source = ["lib/adaptive.c", "parallelepiped.c"]
have_Fq = True
radius_effective_modes = [
    "equivalent cylinder excluded volume", "equivalent volume sphere",
    "half length_a", "half length_b", "half length_c",
    "equivalent circular cross-section", "half ab diagonal", "half diagonal",
    ]

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
