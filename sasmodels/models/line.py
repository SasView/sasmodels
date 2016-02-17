r"""
This model calculates intensity using simple linear function
Definition
----------

The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = A + B \cdot q

.. note::
    For 2D plots intensity has different definition than other shape independent models
.. math::
    I(q) = I(qx) \cdot I(qy)

.. figure:: None

References
----------

None.

"""
from numpy import inf
from numpy import cos
from numpy import sin

name = "line"
title = "Line model"
description = """\
      I(q) = A + B*q

      List of default parameters:
      A = intercept
      B = slope
      """
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["intercept",     "1/cm",        1.0, [-inf, inf], "", "intercept in linear model"],
              ["slope",     "Ang/cm",    1.0, [-inf, inf], "", "slope in linear model"],
              ]
# pylint: enable=bad-whitespace, line-too-long

def Iq(q, intercept, slope):
    """
    :param q:           Input q-value
    :param intercept:   Intrecept in linear model
    :param slope:       Slope in linear model
    :return:            Calculated Intensity
    """
    inten = intercept + slope*q
    return inten

Iq.vectorized = True # Iq accepts an array of q values

def Iqxy(qx, qy, *args):
    """
    :param qx:   Input q_x-value
    :param qy:   Input q_y-value
    :param args: Remaining arguments
    :return:     2D-Intensity
    """
    #TODO: Instrcution tels 2D has different deffinition than oher models
    #return  Iq(qy,*args)*Iq(qy,*args)
    return  Iq(qx, *args)*Iq(qy, *args)

Iqxy.vectorized = True # Iqxy accepts an array of qx, qy values

demo = dict(scale=1.0, background=0, intercept=1.0, slope=1.0)

oldname = "LineModel"
oldpars = dict(intercept='A', slope='B', background=None, scale=None)

tests = [

    [{'intercept':   1.0,
      'slope': 1.0,
     }, 0.4, 1.4],

    [{'intercept':   1.0,
      'slope': 1.0,
     }, 1.3, 2.3],

    [{'intercept':   1.0,
      'slope': 1.0,
     }, 0.5, 1.5],

    [{'intercept':   1.0,
      'slope': 1.0,
     }, [0.4, 0.5], [1.4, 1.5]],

    [{'intercept':   1.0,
      'slope': 1.0,
     }, (1.3, 1.57), 2.30238060425],
]
