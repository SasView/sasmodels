#power_law model
#conversion of PowerLawAbsModel.py
#converted by Steve King, Dec 2015

r"""
This model calculates a simple power law with a flat background.

Definition
----------

.. math::

    I(q) = \text{scale} \cdot q^{-\text{power}} + \text{background}

Note the minus sign in front of the exponent. The exponent *power*
should therefore be entered as a **positive** number for fitting.

Also note that unlike many other models, *scale* in this model
is NOT explicitly related to a volume fraction. Be careful if
combining this model with other models.


References
----------

None.
"""

from numpy import inf, sqrt

name = "power_law"
title = "Simple power law with a flat background"

description = """
    Evaluates the function
    I(q) = scale * q^(-power) + background
    NB: enter power as a positive number!
    """
category = "shape-independent"

#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["power", "", 4.0, [-inf, inf], "", "Power law exponent"]]

# NB: Scale and Background are implicit parameters on every model
def Iq(q, power):
    # pylint: disable=missing-docstring
    inten = (q**-power)
    return inten
Iq.vectorized = True  # Iq accepts an array of q values

def Iqxy(qx, qy, *args):
    # pylint: disable=missing-docstring
    return Iq(sqrt(qx ** 2 + qy ** 2), *args)
Iqxy.vectorized = True # Iqxy accepts an array of qx, qy values

demo = dict(scale=1.0,
            power=4.0,
            background=0.0)

oldname = "PowerLawAbsModel"
oldpars = dict(scale='scale',
               power='m',
               background='background')

tests = [
    [{'scale': 1.0, 'power': 4.0, 'background' : 0.0},
     [0.0106939, 0.469418], [7.64644e+07, 20.5949]],
    ]
