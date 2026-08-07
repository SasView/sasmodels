r"""
Definition
----------

Stabilized power-law structure factor.

WARNING: Should not be used in combination with very anisotropic particle shapes.

This structure factor is intended for describing large aggregates of glubular particles whose overall
dimensions exceed the length scales accessible in the experiment. In such
situations the aggregate size cannot be determined from the measured q-range,
and fitting a specific aggregate model is generally not warranted. However,
these large structures often give rise to a substantial increase in scattering
towards low q.

Over a limited q-range, the structure factor for such unresolved
aggregates is frequently well approximated by a power law,

.. math::

    S(q) = 1 + A q^{-p},

where :math:`p` is the power-law exponent and :math:`A` is a scale factor. The model therefore provides a
simple phenomenological description of excess low-q scattering.

The structure factor is implemented as

.. math::

    S(q)
    =
    1
    +
    A\left(\frac{q_0}{q}\right)^p,

where :math:`A` is *power_law_scale*, :math:`p` is *power*, and

.. math::

    q_0 = 0.01\;\mathrm{\AA}^{-1}.

The reference value :math:`q_0` has no physical significanc. It provides a more stable fit when it is chosen to
lie within the q-range commonly encountered in small-angle scattering
experiments so that the fitted amplitude has a direct interpretation. At the
reference point,

.. math::

    S(q_0)=1+A,

and therefore *power_law_scale* is simply the excess contribution from the
unresolved aggregates at :math:`q=0.01\;\mathrm{\AA}^{-1}`.

This parameterization is mathematically equivalent to

.. math::

    S(q)=1+Bq^{-p},

but, as mentioned, offers a practical advantage during fitting. In the conventional form,
the amplitude and exponent are often strongly correlated, since changing the
exponent changes the magnitude of the power-law contribution throughout the
entire fitted q-range. By introducing the reference value :math:`q_0`, the
amplitude is effectively defined at a q-value inside the experimental range.
The parameter *power_law_scale* is then determined primarily by the scattering
level near :math:`q_0`, whereas *power* is determined primarily by the slope,
which substantially reduces parameter covariance and improves fitting
robustness.

As :math:`q` increases, the unresolved aggregate contribution decreases and

.. math::

    S(q)\rightarrow 1,

so that the scattering becomes dominated by the underlying form factor.

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
