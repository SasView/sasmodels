r"""
Definition
----------

Calculates the form factor for a rectangular solid with a core-shell structure.
The thickness and the scattering length density of the shell or 
"rim" can be different on each (pair) of faces. However at this time
the 1D calculation does **NOT** actually calculate a c face rim despite the presence of
the parameter. Some other aspects of the 1D calculation may be wrong.

.. note::
   This model was originally ported from NIST IGOR macros. However, it is not
   yet fully understood by the SasView developers and is currently under review.

The form factor is normalized by the particle volume $V$ such that

.. math::

    I(q) = \text{scale}\frac{\langle f^2 \rangle}{V} + \text{background}

where $\langle \ldots \rangle$ is an average over all possible orientations
of the rectangular solid.


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

**meaning that there are "gaps" at the corners of the solid.**  Again note that
$t_C = 0$ currently. 

The intensity calculated follows the :ref:`parallelepiped` model, with the
core-shell intensity being calculated as the square of the sum of the
amplitudes of the core and shell, in the same manner as a core-shell model.

.. math::

    F_{a}(Q,\alpha,\beta)=
    \left[\frac{\sin(\tfrac{1}{2}Q(L_A+2t_A)\sin\alpha \sin\beta)}{\tfrac{1}{2}Q(L_A+2t_A)\sin\alpha\sin\beta}
    - \frac{\sin(\tfrac{1}{2}QL_A\sin\alpha \sin\beta)}{\tfrac{1}{2}QL_A\sin\alpha \sin\beta} \right]
    \left[\frac{\sin(\tfrac{1}{2}QL_B\sin\alpha \sin\beta)}{\tfrac{1}{2}QL_B\sin\alpha \sin\beta} \right]
    \left[\frac{\sin(\tfrac{1}{2}QL_C\sin\alpha \sin\beta)}{\tfrac{1}{2}QL_C\sin\alpha \sin\beta} \right]

.. note::

    Why does t_B not appear in the above equation?
    For the calculation of the form factor to be valid, the sides of the solid
    MUST (perhaps not any more?) be chosen such that** $A < B < C$.
    If this inequality is not satisfied, the model will not report an error,
    but the calculation will not be correct and thus the result wrong.

FITTING NOTES
If the scale is set equal to the particle volume fraction, $\phi$, the returned
value is the scattered intensity per unit volume, $I(q) = \phi P(q)$.
However, **no interparticle interference effects are included in this
calculation.**

There are many parameters in this model. Hold as many fixed as possible with
known values, or you will certainly end up at a solution that is unphysical.

Constraints must be applied during fitting to ensure that the inequality
$A < B < C$ is not violated. The calculation will not report an error,
but the results will not be correct.

The returned value is in units of |cm^-1|, on absolute scale.

NB: The 2nd virial coefficient of the core_shell_parallelepiped is calculated
based on the the averaged effective radius $(=\sqrt{(A+2t_A)(B+2t_B)/\pi})$
and length $(C+2t_C)$ values, after appropriately
sorting the three dimensions to give an oblate or prolate particle, to give an 
effective radius, for $S(Q)$ when $P(Q) * S(Q)$ is applied.

To provide easy access to the orientation of the parallelepiped, we define the
axis of the cylinder using three angles $\theta$, $\phi$ and $\Psi$.
(see :ref:`cylinder orientation <cylinder-angle-definition>`).
The angle $\Psi$ is the rotational angle around the *long_c* axis against the
$q$ plane. For example, $\Psi = 0$ when the *short_b* axis is parallel to the
*x*-axis of the detector.

.. figure:: img/parallelepiped_angle_definition.png

    Definition of the angles for oriented core-shell parallelepipeds.

.. figure:: img/parallelepiped_angle_projection.png

    Examples of the angles for oriented core-shell parallelepipeds against the
    detector plane.

References
----------

.. [#] P Mittelbach and G Porod, *Acta Physica Austriaca*, 14 (1961) 185-211
    Equations (1), (13-14). (in German)
.. [#] D Singh (2009). *Small angle scattering studies of self assembly in
   lipid mixtures*, John's Hopkins University Thesis (2009) 223-225. `Available
   from Proquest <http://search.proquest.com/docview/304915826?accountid
   =26379>`_

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Converted to sasmodels by:** Miguel Gonzales **Date:** February 26, 2016
* **Last Modified by:** Wojciech Potrzebowski **Date:** January 11, 2017
* **Currently Under review by:** Paul Butler
"""

import numpy as np
from numpy import pi, inf, sqrt, cos, sin

name = "core_shell_parallelepiped"
title = "Rectangular solid with a core-shell structure."
description = """
     P(q)= 
"""
category = "shape:parallelepiped"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld_core", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Parallelepiped core scattering length density"],
              ["sld_a", "1e-6/Ang^2", 2, [-inf, inf], "sld",
               "Parallelepiped A rim scattering length density"],
              ["sld_b", "1e-6/Ang^2", 4, [-inf, inf], "sld",
               "Parallelepiped B rim scattering length density"],
              ["sld_c", "1e-6/Ang^2", 2, [-inf, inf], "sld",
               "Parallelepiped C rim scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 6, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["length_a", "Ang", 35, [0, inf], "volume",
               "Shorter side of the parallelepiped"],
              ["length_b", "Ang", 75, [0, inf], "volume",
               "Second side of the parallelepiped"],
              ["length_c", "Ang", 400, [0, inf], "volume",
               "Larger side of the parallelepiped"],
              ["thick_rim_a", "Ang", 10, [0, inf], "volume",
               "Thickness of A rim"],
              ["thick_rim_b", "Ang", 10, [0, inf], "volume",
               "Thickness of B rim"],
              ["thick_rim_c", "Ang", 10, [0, inf], "volume",
               "Thickness of C rim"],
              ["theta", "degrees", 0, [-360, 360], "orientation",
               "c axis to beam angle"],
              ["phi", "degrees", 0, [-360, 360], "orientation",
               "rotation about beam"],
              ["psi", "degrees", 0, [-360, 360], "orientation",
               "rotation about c axis"],
             ]

source = ["lib/gauss76.c", "core_shell_parallelepiped.c"]


def ER(length_a, length_b, length_c, thick_rim_a, thick_rim_b, thick_rim_c):
    """
        Return equivalent radius (ER)
    """

    # surface average radius (rough approximation)
    surf_rad = sqrt((length_a + 2.0*thick_rim_a) * (length_b + 2.0*thick_rim_b) / pi)

    height = length_c + 2.0*thick_rim_c

    ddd = 0.75 * surf_rad * (2 * surf_rad * height + (height + surf_rad) * (height + pi * surf_rad))
    return 0.5 * (ddd) ** (1. / 3.)

# VR defaults to 1.0

# parameters for demo
demo = dict(scale=1, background=0.0,
            sld_core=1, sld_a=2, sld_b=4, sld_c=2, sld_solvent=6,
            length_a=35, length_b=75, length_c=400,
            thick_rim_a=10, thick_rim_b=10, thick_rim_c=10,
            theta=0, phi=0, psi=0,
            length_a_pd=0.1, length_a_pd_n=1,
            length_b_pd=0.1, length_b_pd_n=1,
            length_c_pd=0.1, length_c_pd_n=1,
            thick_rim_a_pd=0.1, thick_rim_a_pd_n=1,
            thick_rim_b_pd=0.1, thick_rim_b_pd_n=1,
            thick_rim_c_pd=0.1, thick_rim_c_pd_n=1,
            theta_pd=10, theta_pd_n=1,
            phi_pd=10, phi_pd_n=1,
            psi_pd=10, psi_pd_n=1)

# rkh 7/4/17 add random unit test for 2d, note make all params different, 2d values not tested against other codes or models
qx, qy = 0.2 * cos(pi/6.), 0.2 * sin(pi/6.)
tests = [[{}, 0.2, 0.533149288477],
         [{}, [0.2], [0.533149288477]],
         [{'theta':10.0, 'phi':20.0}, (qx, qy), 0.0853299803222],
         [{'theta':10.0, 'phi':20.0}, [(qx, qy)], [0.0853299803222]],
        ]
del qx, qy  # not necessary to delete, but cleaner
