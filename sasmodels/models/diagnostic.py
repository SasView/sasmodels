r"""
Gaussian function for diagnostic purposes
"""

import numpy as np
from numpy import inf, expm1, power

name = "diganostic"
title = "Diagnostic model"

description = """ """
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["mean", "1/Ang", 0, [-inf, inf], "volume", ""],
    ["std",  "1/Ang", 1, [0.0, inf], "", ""],
    ]
# pylint: enable=bad-whitespace, line-too-long

def Iq(q, mean, std):
    return np.exp(-0.5*((q - mean)/std)**2)

Iq.vectorized = True  # Iq accepts an array of q values


# these unit test values taken from SasView 3.1.2
tests = []
