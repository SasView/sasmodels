#power_law model
#conversion of PowerLawAbsModel.py

r"""
This model calculates a simple power law with a flat background.

Definition
----------

.. math::

    I(q) = \text{scale} \cdot q^{-m} + \text{background}

Note the minus sign in front of the exponent. The exponent *m* 
should therefore be entered as a **positive** number.

Also note that unlike many other models, *scale* in this model 
is NOT explicitly related to a volume fraction. Be careful if 
combining this model with other models.

.. image:: img/power_law_1d.jpg

*Figure. 1D plot using the default values (w/200 data point).*

REFERENCE

None.
"""

from numpy import inf, sqrt

name =  "power_law"
title = "Simple power law with a flat background"

description = """\
        Evaluates the function
        I(q) = scale * q^(-m) + bkgd
        NB: enter m as a positive number!
        """
category = "shape-independent"

#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["scale", "", 1.0, [-inf, inf], "", "Power law scale factor"],
              ["m", "", 4.0, [-inf, inf], "", "Power law exponent"],
              ["bkgd", "", 0.0, [-inf, inf], "", "Background level"]]

def Iq(scale,m,bkgd):
    inten = (scale * q ** m + bkgd)
    return inten
Iq.vectorized = True  # Iq accepts an array of Q values

def Iqxy(qx, qy, *args):
    return Iq(sqrt(qx ** 2 + qy ** 2), *args)
Iqxy.vectorized = True # Iqxy accepts an array of Qx, Qy values

demo = dict(scale=1.0,
            m=4.0,
            bkgd=0.0)

oldname = "PowerLawAbsModel"
oldpars = dict(scale='scale',
               m='m',
               bkgd='background')

tests = [
        [ {'scale': 1.0, 'm': 4.0, 'bkgd' : 0.0}, [0.0106939, 0.469418], [7.64644e+07, 20.5949]]
        ]
