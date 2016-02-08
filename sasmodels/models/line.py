r"""
This model calculates intensity using simple linear function
Definition
----------

The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = A + B*q

.. figure:: None

References
----------

None.

"""
from numpy import inf

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
parameters = [["intercept",     "cm^-1",        1.0, [-inf, inf], "", "intercept in linear model"],
              ["slope",     "Ang*cm^-1",    1.0, [-inf, inf], "", "slope in linear model"],
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
    #return Iq(sqrt(qx ** 2 + qy ** 2), *args)
    return  Iq(qx, *args)*Iq(qy, *args)

Iqxy.vectorized = True # Iqxy accepts an array of qx, qy values

demo = dict(intercept=1.0, slope=1.0)

oldname = "LineModel"
oldpars = dict(intercept='A', slope='B', scale=None, background=None)

tests = [
    # Accuracy tests based on content in test/utest_other_models.py
    [{'intercept':   1.0,
      'slope': 1.0,
     }, 0.4, 1.4],

    [{'intercept':   1.0,
      'slope': 1.0,
     }, 1.3, 2.3],

    [{'intercept':   1.0,
      'slope': 1.0,
     },0.5, 1.5],

    [{'intercept':   1.0,
      'slope': 1.0,
     }, [0.4,0.5], [1.4,1.5]],

    [{'intercept':   1.0,
      'slope': 1.0,
     }, [1.3,1.57], [2.3,2.57]],
]