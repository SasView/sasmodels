#!/usr/bin/env python
"""
Small application to change theta, phi and psi from SasView 3.x models to the
new angle definition in SasView 4.x and above.

Usage: python explore/transform_angles.py theta phi psi
"""
from __future__ import print_function, division

import sys

import numpy as np
from numpy import pi, cos, sin, sqrt, exp, degrees, radians
from scipy.optimize import fmin

# Definition of rotation matrices comes from wikipedia:
#    https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations
def Rx(angle):
    """Construct a matrix to rotate points about *x* by *angle* degrees."""
    a = radians(angle)
    R = [[1, 0, 0],
         [0, +cos(a), -sin(a)],
         [0, +sin(a), +cos(a)]]
    return np.array(R)

def Ry(angle):
    """Construct a matrix to rotate points about *y* by *angle* degrees."""
    a = radians(angle)
    R = [[+cos(a), 0, +sin(a)],
         [0, 1, 0],
         [-sin(a), 0, +cos(a)]]
    return np.array(R)

def Rz(angle):
    """Construct a matrix to rotate points about *z* by *angle* degrees."""
    a = radians(angle)
    R = [[+cos(a), -sin(a), 0],
         [+sin(a), +cos(a), 0],
         [0, 0, 1]]
    return np.array(R)


def transform_angles(theta, phi, psi, qx=0.1, qy=0.1):
    Rold = Rz(-psi)@Rx(theta)@Ry(-(90 - phi))
    cost = lambda p: np.linalg.norm(Rz(-p[2])@Ry(-p[0])@Rz(-p[1]) - Rold)
    result = fmin(cost, (theta, phi, psi))
    theta_p, phi_p, psi_p = result
    Rnew = Rz(-psi_p)@Ry(-theta_p)@Rz(-phi_p)

    print("old: theta, phi, psi =", ", ".join(str(v) for v in (theta, phi, psi)))
    print("new: theta, phi, psi =", ", ".join(str(v) for v in result))
    try:
        point = np.array([qx, qy, [0]*len(qx)])
    except TypeError:
        point = np.array([[qx],[qy],[0]])
    for p in point.T:
        print("q abc old for", p, (Rold@p.T).T)
        print("q abc new for", p, (Rnew@p.T).T)

if __name__ == "__main__":
    theta, phi, psi = (float(v) for v in sys.argv[1:])
    #transform_angles(theta, phi, psi)
    transform_angles(theta, phi, psi, qx=-0.017, qy=0.035)
