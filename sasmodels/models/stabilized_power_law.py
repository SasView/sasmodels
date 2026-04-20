r"""
Definition
----------

Stabilized power_law structure factor.

WARNING: Should not be used in combination with very anisotrpic particle shapes.

It uses rhe equation:

.. math::

    S(q) = 1.0 + amp (0.01/q)^pow

where *amp* is the scale of the power law and *pow* is the exponent.

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
I(q) = 1 + amp *(0.01/q)^pow
    amp: scale of power law
    pow: exponent of power law
"""
category = "structure-factor"
structure_factor = False

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["amp", "", 100, [0, inf], "", "scale of power law"],
              ["pow", "", 2, [0, 6], "", "exponent of power law "]
             ]

def Iq(q, amp=100, pow=2):
    """
    Parameters:
        q:      input scattering vectors, units 1/Ang
        amp:      amplitude of power law , default value=2
        pow:      exponent of power law , default value=2
    Returns:
        S(q):   1D scattering intensity at q, units none
    """
    iq = 1.0 + amp * (0.01/q) ** pow
    return iq

# include tests for your model
# tests = [
#     [{'m': 2.0, 'b' : 1.0}, [q1, q2], [expected_Iq1, expected_Iq2]],
#     [{'m': 3.0, 'b' : 5.0}, [q3, q4], [expected_Iq3, expected_Iq4]],
# ]
