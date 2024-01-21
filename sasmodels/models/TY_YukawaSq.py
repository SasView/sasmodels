r"""
Definition
----------

Calculates test2yv2.



References
----------

Authorship and Verification
---------------------------

* **Author:** --- **Date:** 2021YYY-05m-19d
* **Last Modified by:** --- **Date:** 2021YYY-05m-19d
* **Last Reviewed by:** --- **Date:** 2021YYY-05m-19d
"""

from sasmodels.special import *
from numpy import inf

name = "TY_YukawaSq"
title = "User model for two Yukawa structure factor (S(q))"
description = """"""

category = "structure-factor"
structure_factor = True

single = False

parameters = [ 
#   ["name", "units", default, [lower, upper], "type", "description"],
    ["radius_effective", "Ang", 50.0, [-inf, inf], '', ''],
    ['volfraction', '', 0.1, [-inf, inf], '', ''],
    ['k1', '', 6, [-inf, inf], '', ''],
    ['k2', '', -2.0, [-inf, inf], '', ''],
    ['z1', '', 10.0, [-inf, inf], '', ''],
    ['z2', '', 2.0, [-inf, inf], '', ''],
    ]
#def Iq(x, radius, volumefraction, k1, k2, z1, z2):
#    """Absolute scattering"""
#    q=x
#    
#    
#    return q


## uncomment the following if Iq works for vector x

source = ["TY_PairCorrelation.c", "TY_cpoly.c","TY_TwoYukawa.c", "TY_utility.c", "TY_YukawaSq.c"]

#source = ["test2yv2_inc.c"]

#haveFq = True   # for beta calculation
#single = False  # make sure that it has double digit precision
#opencl = False  # related with parallel computing

#Iq.vectorized = True

#def Iqxy(x, y, radius, volumefraction, k1, k2, z1, z2):
#    """Absolute scattering of oriented particles."""
#    ...
#    return oriented_form(x, y, args)
## uncomment the following if Iqxy works for vector x, y
#Iqxy.vectorized = True
