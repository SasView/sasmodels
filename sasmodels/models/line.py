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
parameters = [["intercept", "1/cm",   1.0, [-inf, inf], "", "intercept in linear model"],
              ["slope",     "Ang/cm", 1.0, [-inf, inf], "", "slope in linear model"],
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
    # TODO: In SasView code additional formula for list has been specified.
    # if inten(x) = intercept + slope*x:
    # then if q is a list, Iq=inten(x[0]*math.cos(x[1]))*inten(x[0]*math.sin(x[1]))
    return inten

Iq.vectorized = True # Iq accepts an array of q values

def Iqxy(qx, qy, *args):
    """
    :param qx:   Input q_x-value
    :param qy:   Input q_y-value
    :param args: Remaining arguments
    :return:     2D-Intensity
    """
    # TODO: SasView documention lists 2D intensity as Iq(qx)*Iq(qy) but code says:
    # return self._line(x[1])
    return Iq(qx, *args)*Iq(qy, *args)

Iqxy.vectorized = True  # Iqxy accepts an array of qx, qy values

demo = dict(scale=1.0, background=0, intercept=1.0, slope=1.0)

tests = [

    [{'intercept':   1.0,
      'slope': 1.0,
     }, 1.0, 2.001],

    [{'intercept':   1.0,
      'slope': 1.0,
     }, 0.0, 1.001],

    [{'intercept':   1.0,
      'slope': 1.0,
     }, 0.4, 1.401],

    [{'intercept':   1.0,
      'slope': 1.0,
     }, 1.3, 2.301],

    [{'intercept':   1.0,
      'slope': 1.0,
     }, 0.5, 1.501],

    [{'intercept':   1.0,
      'slope': 1.0,
     }, [0.4, 0.5], [1.401, 1.501]],

    [{'intercept':   1.0,
      'slope': 1.0,
      'background': 0.0,
     }, (1.3, 1.57), 5.911],
]
