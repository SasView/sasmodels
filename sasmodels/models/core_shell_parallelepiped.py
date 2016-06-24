# core_shell_parallelepiped model
# Note: model title and parameter table are inserted automatically
r"""
Calculates the form factor for a rectangular solid with a core-shell structure.
**The thickness and the scattering length density of the shell or "rim"
can be different on all three (pairs) of faces.**

The form factor is normalized by the particle volume $V$ such that

.. math::

    I(q) = \text{scale}\frac{\langle f^2 \rangle}{V} + \text{background}

where $\langle \ldots \rangle$ is an average over all possible orientations
of the rectangular solid.

An instrument resolution smeared version of the model is also provided.


Definition
----------

The function calculated is the form factor of the rectangular solid below.
The core of the solid is defined by the dimensions $A$, $B$, $C$ such that
$A < B < C$.

.. image:: img/core_shell_parallelepiped_geometry.jpg

There are rectangular "slabs" of thickness $t_A$ that add to the $A$ dimension
(on the $BC$ faces). There are similar slabs on the $AC$ $(=t_B)$ and $AB$
$(=t_C)$ faces. The projection in the $AB$ plane is then

.. image:: img/core_shell_parallelepiped_projection.jpg

The volume of the solid is

.. math::

    V = ABC + 2t_ABC + 2t_BAC + 2t_CAB

**meaning that there are "gaps" at the corners of the solid.**

The intensity calculated follows the :ref:`parallelepiped` model, with the core-shell
intensity being calculated as the square of the sum of the amplitudes of the
core and shell, in the same manner as a core-shell model.

**For the calculation of the form factor to be valid, the sides of the solid
MUST be chosen such that** $A < B < C$.
**If this inequality is not satisfied, the model will not report an error,
and the calculation will not be correct.**

FITTING NOTES
If the scale is set equal to the particle volume fraction, |phi|, the returned
value is the scattered intensity per unit volume, $I(q) = \phi P(q)$.
However, **no interparticle interference effects are included in this calculation.**

There are many parameters in this model. Hold as many fixed as possible with
known values, or you will certainly end up at a solution that is unphysical.

Constraints must be applied during fitting to ensure that the inequality
$A < B < C$ is not violated. The calculation will not report an error,
but the results will not be correct.

The returned value is in units of |cm^-1|, on absolute scale.

NB: The 2nd virial coefficient of the core_shell_parallelepiped is calculated
based on the the averaged effective radius $(=\sqrt{(A+2t_A)(B+2t_B)/\pi})$
and length $(C+2t_C)$ values, and used as the effective radius
for $S(Q)$ when $P(Q) * S(Q)$ is applied.

.. Comment by Miguel Gonzalez:
   The later seems to contradict the previous statement that interparticle interference
   effects are not included.

To provide easy access to the orientation of the parallelepiped, we define the
axis of the cylinder using three angles $\theta$, $\phi$ and $\Psi$.
(see :ref:`cylinder orientation <cylinder-angle-definition>`).
The angle $\Psi$ is the rotational angle around the *long_c* axis against the
$q$ plane. For example, $\Psi = 0$ when the *short_b* axis is parallel to the
*x*-axis of the detector.

.. figure:: img/parallelepiped_angle_definition.jpg

    Definition of the angles for oriented core-shell parallelepipeds.

.. figure:: img/parallelepiped_angle_projection.jpg

    Examples of the angles for oriented core-shell parallelepipeds against the
    detector plane.

Validation
----------

The model uses the form factor calculations implemented in a c-library provided
by the NIST Center for Neutron Research (Kline, 2006).

References
----------

P Mittelbach and G Porod, *Acta Physica Austriaca*, 14 (1961) 185-211
Equations (1), (13-14). (in German)

"""

import numpy as np
from numpy import pi, inf, sqrt

name = "core_shell_parallelepiped"
title = "Rectangular solid with a core-shell structure."
description = """
     P(q)= 
"""
category = "shape:parallelepiped"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld_core", "1e-6/Ang^2", 1, [-inf, inf], "",
               "Parallelepiped core scattering length density"],
              ["sld_a", "1e-6/Ang^2", 2, [-inf, inf], "",
               "Parallelepiped A rim scattering length density"],
              ["sld_b", "1e-6/Ang^2", 4, [-inf, inf], "",
               "Parallelepiped B rim scattering length density"],
              ["sld_c", "1e-6/Ang^2", 2, [-inf, inf], "",
               "Parallelepiped C rim scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 6, [-inf, inf], "",
               "Solvent scattering length density"],
              ["a_side", "Ang", 35, [0, inf], "volume",
               "Shorter side of the parallelepiped"],
              ["b_side", "Ang", 75, [0, inf], "volume",
               "Second side of the parallelepiped"],
              ["c_side", "Ang", 400, [0, inf], "volume",
               "Larger side of the parallelepiped"],
              ["arim_thickness", "Ang", 10, [0, inf], "volume",
               "Thickness of A rim"],
              ["brim_thickness", "Ang", 10, [0, inf], "volume",
               "Thickness of B rim"],
              ["crim_thickness", "Ang", 10, [0, inf], "volume",
               "Thickness of C rim"],
              ["theta", "degrees", 0, [-inf, inf], "orientation",
               "In plane angle"],
              ["phi", "degrees", 0, [-inf, inf], "orientation",
               "Out of plane angle"],
              ["psi", "degrees", 0, [-inf, inf], "orientation",
               "Rotation angle around its own c axis against q plane"],
             ]

source = ["lib/gauss76.c", "core_shell_parallelepiped.c"]


def ER(a_side, b_side, c_side, arim_thickness, brim_thickness, crim_thickness):
    """
        Return equivalent radius (ER)
    """

    # surface average radius (rough approximation)
    surf_rad = sqrt((a_side + 2.0*arim_thickness) * (b_side + 2.0*brim_thickness) / pi)

    height = c_side + 2.0*crim_thickness

    ddd = 0.75 * surf_rad * (2 * surf_rad * height + (height + surf_rad) * (height + pi * surf_rad))
    return 0.5 * (ddd) ** (1. / 3.)

# VR defaults to 1.0

# parameters for demo
demo = dict(scale=1, background=0.0,
            sld_core=1e-6, sld_a=2e-6, sld_b=4e-6,
            sld_c=2e-6, sld_solvent=6e-6,
            a_side=35, b_side=75, c_side=400,
            arim_thickness=10, brim_thickness=10, crim_thickness=10,
            theta=0, phi=0, psi=0,
            a_side_pd=0.1, a_side_pd_n=1,
            b_side_pd=0.1, b_side_pd_n=1,
            c_side_pd=0.1, c_side_pd_n=1,
            arim_thickness_pd=0.1, arim_thickness_pd_n=1,
            brim_thickness_pd=0.1, brim_thickness_pd_n=1,
            crim_thickness_pd=0.1, crim_thickness_pd_n=1,
            theta_pd=10, theta_pd_n=1,
            phi_pd=10, phi_pd_n=1,
            psi_pd=10, psi_pd_n=10)

qx, qy = 0.2 * np.cos(2.5), 0.2 * np.sin(2.5)
tests = [[{}, 0.2, 0.533149288477],
         [{}, [0.2], [0.533149288477]],
         [{'theta':10.0, 'phi':10.0}, (qx, qy), 0.032102135569],
         [{'theta':10.0, 'phi':10.0}, [(qx, qy)], [0.032102135569]],
        ]
del qx, qy  # not necessary to delete, but cleaner
