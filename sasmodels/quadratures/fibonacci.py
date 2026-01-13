# Fibonacci quadrature
# number of points on the unit sphere is adjustable
r"""
This module provides a function to generate quasi-uniformly distributed points on the unit sphere
using the Fibonacci quadrature. These points can be used for numerical integration over the sphere,
which is particularly useful for orientational averaging in scattering calculations.
It also includes a plotting function to visualize the distribution of points on the sphere.
The only parameter is the number of points to generate. Around 500 points provide a good balance
between accuracy and computational efficiency for most applications.
"""


import matplotlib.pyplot as plt
import numpy as np


def fibonacci_sphere(npoints_fibonacci: int):
    """
    Generates npoints quasi-uniformly distributed on the unit sphere
    in Cartesian coordinates (x,y,z) and their associated weights.
    Parameters
    ----------
    npoints : int
        Number of points to generate.
    Returns
    -------
    points : ndarray, shape (npoints, 3)
        Cartesian coordinates of the points on the unit sphere.
    weights : ndarray, shape (npoints,)
        Weights associated with each point for integration on the sphere.
    """

    indices = np.arange(0, npoints_fibonacci, dtype=float) + 0.5
    phi = np.arccos(1.0 - 2.0 * indices / npoints_fibonacci)
    theta = 2.0 * np.pi * indices / ((1.0 + 5.0 ** 0.5) / 2.0)

    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)

    weights = np.full(len(x), 1.0 / npoints_fibonacci)

    return np.column_stack((x, y, z)), weights

def plot_fibonacci_sphere(npoints_fibonacci=500, figsize=(7, 7)):
    """
    3D representation of Fibonacci points on the unit sphere.
    Parameters
    ----------
    npoints : int
        Number of points to generate and display.
    figsize : tuple
        Size of the figure.

    """
    pts, w = fibonacci_sphere(npoints_fibonacci)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=10, alpha=0.6)

    # Axes
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(f"Fibonacci points on the unit sphere, ({npoints_fibonacci} total points)")
    ax.set_box_aspect([1, 1, 1])
    plt.show()
