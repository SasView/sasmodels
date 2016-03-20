r"""
This is a test function!

Definition
----------

Calculates a bessel function. Maybe...

References
----------

None
"""

from numpy import inf

name = "bessel"
title = "Bessel function testing"
description = """\
Leveraging current infrastructure to test Bessel function performance on
"""
category = "special_functions:bessel"

#             ["name", "units", default, [lower, upper], "type","description"],
#Bessel
parameters = [
    ["ignored", "", 0.0, [-inf, inf], "", "no parameterless functions"],
             ]

source = ["lib/sas_gamma.c"]

# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    """

Iq = """
    return sas_gamma(q);
    """

Iqxy = """
    // never called since no orientation or magnetic parameters.
    //return -1.0;
    """

# VR defaults to 1.0

demo = dict(scale=1, background=0,
            )
oldname = "Bessel"
oldpars = dict()
