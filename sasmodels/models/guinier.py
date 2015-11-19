r"""
Guinier (Model)

Definition
----------

This model fits the Guinier function

.. math:: Q_1=\frac{1}{R_g}\sqrt{\frac{(m-s)(3-s)}{2}}

to the data directly without any need for linearisation (*cf*. Ln *I(q)* vs *q*\ :sup:`2`).

For 2D data: The 2D scattering intensity is calculated in the same way as 1D, where the *q* vector is defined as

.. math:: q=\sqrt{q_x^2 + q_y^2}

REFERENCE

A Guinier and G Fournet, *Small-Angle Scattering of X-Rays*, John Wiley & Sons, New York (1955)
"""

from numpy import inf

name = "guinier"
title = ""
description = """
 I(q) = scale exp ( - rg^2 q^2 / 3.0 )
 
    List of default parameters:
    scale = scale
    rg = Radius of gyration
"""
category = "shape-independent"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [
              ["rg", "Ang", 60.0, [0, inf], "", "Radius of Gyration"],
              ]

Iq = """
    double exponent = rg*rg*q*q/3.0;
    double value = exp(-exponent);
    return value;
"""

Iqxy = """
    return Iq(sqrt(qx*qx + qy*qy), rg);
    """

# parameters for demo
demo = dict(scale=1.0,rg=60.0)

# For testing against the old sasview models, include the converted parameter
# names and the target sasview model name.
oldname = 'GuinierModel'
oldpars = dict(rg='rg')

# parameters for unit tests
tests = [
         [{'rg' : 31.5}, 0.005, 0.991756]
         ]