r"""
Definition
----------

Stabilized power-law structure factor.

WARNING: Should not be used in combination with very anisotropic particle shapes.

It uses the equation:

.. math::

    S(q) = 1.0 + \mathrm{power_law_scale}\,(0.01/q)^{\mathrm{power}}

where *power_law_scale* is the scale of the power law and *power* is the exponent.

The first two parameters follow sasmodels structure-factor conventions for
:math:`P@S` products: *radius_effective* (effective scatterer radius, unused
here) and *volfraction* are **not** used in the formula above.

See also
Larsen, A. H., Pedersen, J. S., & Arleth, L. (2020). Assessment of
structure factors for analysis of small-angle scattering data from
desired or undesired aggregates. Applied Crystallography, 53(4), 991-1005.

Validation
----------


References
----------

#. Jan Skov Pedersen

Authorship and Verification
----------------------------

* **Author:** Jan Skov Pedersen

* **Last Modified by:** Jan Skov Pedersen, April 12, 2026
* **Last Reviewed by:** Reviewer Name Here **Date:** Date Here
"""

from numpy import inf

name = "stabilized_power_law"
title = "Stabilized power-law structure factor"
description = """\
S(q) = 1 + power_law_scale *(0.01/q)^power
    power_law_scale: scale of power law
    power: exponent of power law
"""
category = "structure-factor"
structure_factor = True

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["radius_effective", "Ang", 50.0, [0.0, inf], "",
     "effective scatterer radius; unused in S(q), required for P@S products"],
    ["volfraction", "", 0.2, [0.0, 1.0], "",
     "unused in S(q); required for P@S products"],
    ["power_law_scale", "", 100, [0, inf], "", "scale of power law"],
    ["power", "", 2, [0, 6], "", "exponent of power law "],
]


def Iq(q, radius_effective, volfraction, power_law_scale, power):
    """Return S(q); *radius_effective* and *volfraction* are unused."""
    _ = (radius_effective, volfraction)
    return 1.0 + power_law_scale * (0.01 / q) ** power
