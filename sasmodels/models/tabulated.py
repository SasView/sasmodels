r"""
Definition
----------

References
----------

None.

Authorship and Verification
---------------------------
* **Author:** Jose Borreguero **Date:** August 26 2017
* **Last Modified by:** Jose Borreguero **Date:** July 27 2016
* **Last Reviewed by:** _REVIEWER_NAME_HERE **Date:** _DATE_HERE_
"""

import inspect
from scipy.interpolate import interp1d, interp2d
import numpy as np

name = "tabulated"
title = "An interpolator of a numeric structure factor"
description = """\
    I(q) = scale * interpolator(q|q', I')
    List of default parameters:
    scale = scaling factor"""
category = "shape-independent"


# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["scale",  "",      1.0,       [0, inf],     "",  "Scaling factor"],]
# pylint: enable=bad-whitespace, line-too-long


# Attributes of the model
_attr = {'itp': None,  # interpolator
         'q':   None,  # table of momentum transfer values
         'qx':  None,  # table of momentum transfer values along X-axis
         'qy':  None,  # table of momentum transfer values along Y-axis
         'I':   None}  # table of intensities


def update_interpolator(itp=None):
    """
    Update the interpolator with custom. 'None' for linear interpolation.
    :param itp: interpolator generator. See interp1D and interp2D of
    scipy.interpolate for examples.
    """
    if None not in (_attr['q'], _attr['I']):
        interpolator = itp if itp else interp1d
        _attr['itp'] = interpolator(_attr['q'], _attr['I'])
    elif None not in (_attr['qx'], _attr['qy'], _attr['I']):
        interpolator = itp if itp else interp2d
        _attr['itp'] = interpolator(_attr['qx'], _attr['qx'], _attr['I'])


def initialize(q=None, qx=None, qy=None, I=None):
    """
    Initialize any of the function attributes
    :param q:     sequence of momentum transfer values
    :param qx:    sequence of momentum transfer values along X-axis
    :param qy:    sequence of momentum transfer values along X-axis
    :param I:     sequence of intensities
    """
    frame = inspect.currentframe()
    attrs, _, _, values = inspect.getargvalues(frame)
    for attr in attrs:
        value = values[attr]
        if value is not None:
            _attr[attr] = np.array(value, dtype=np.float)
    update_interpolator()


def numpyze(calc_I):
    """
    Transform positional arguments into numpy arrays
    :param calc_I: Intensities calculator
    :return: function acting on the transformed positional arguments
    """
    def wrapper(*args, **kwargs):
        numpy_args = [np.array(arg, dtype=np.float) for arg in args]
        return calc_I(*numpy_args, **kwargs)


@numpyze
def Iq(q, scale=1.0):
    """
    :param q:     sequence of momentum transfer values
    :param scale: scaling factor
    :return:      calculated intensities
    """
    return scale * _attr['itp'](q)

Iq.vectorized = True  # Iq accepts an array of q values


@numpyze
def Iqxy(qx, qy, scale=1.0):
    """
    :param qx:    sequence of momentum transfer values along X-axis
    :param qy:    sequence of momentum transfer values along X-axis
    :param scale: scaling factor
    :return:      calculated intensities
    """
    return scale * _attr['itp'](qx, qy)
    pass

Iqxy.vectorized = True  # Iqxy accepts arrays of qx, qy values


