"""Regenerate the oriented prism and tetrahedron reference intensities.

Run this script with NumPy and SciPy. It integrates exp(i q.r) over each
particle in real space, independently of the sasmodels kernels. Lengths are
in angstroms, q values are in inverse angstroms, and angles are in degrees.
"""

import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.spatial.transform import Rotation

Q_XY = np.array([[0.013, 0.021], [-0.027, 0.014], [0.035, -0.023], [-0.040, 0.005]])
VIEWS = [(0.0, 0.0, 0.0), (37.0, -22.0, 19.0)]


def unit_cube(order):
    """Return tensor-product Gauss-Legendre nodes and weights on [0, 1]^3."""
    nodes, weights = leggauss(order)
    nodes = (nodes + 1) / 2
    weights = weights / 2
    u, v, w = np.meshgrid(nodes, nodes, nodes, indexing="ij")
    wu, wv, ww = np.meshgrid(weights, weights, weights, indexing="ij")
    return u.ravel(), v.ravel(), w.ravel(), (wu * wv * ww).ravel()


def prism(order, sides=5, radius=24.0, length=110.0):
    """Map the cubature onto triangular wedges of an equal-area prism."""
    u, v, w, weights = unit_cube(order)
    angle = 2 * math.pi / sides
    circumradius = math.sqrt(2 * math.pi * radius**2 / (sides * math.sin(angle)))
    angles = np.arange(sides) * angle
    polygon = circumradius * np.column_stack((np.cos(angles), np.sin(angles)))
    points, volumes = [], []
    for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
        xy = u[:, None] * start + ((1 - u) * v)[:, None] * end
        points.append(np.column_stack((xy, length * (w - 0.5))))
        volumes.append(weights * (1 - u) * abs(np.linalg.det([start, end])) * length)
    return np.concatenate(points), np.concatenate(volumes)


def tetrahedron(order, radius=52.0):
    """Map the cubature onto a regular tetrahedron with the given circumradius."""
    u, v, w, weights = unit_cube(order)
    edges = (2 * radius / math.sqrt(3)) * np.array(
        [[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]]
    )
    barycentric = np.column_stack((u, (1 - u) * v, (1 - u) * (1 - v) * w))
    jacobian = (1 - u)**2 * (1 - v) * abs(np.linalg.det(edges))
    return barycentric @ edges, weights * jacobian


def intensity(points, weights, view):
    theta, phi, psi = view
    # Rotate particle coordinates into the laboratory frame, where q_z = 0.
    detector_xy = Rotation.from_euler("ZYZ", [phi, theta, psi], degrees=True).apply(points)[:, :2]
    amplitude = np.exp(1j * (Q_XY @ detector_xy.T)) @ weights
    volume = weights.sum()
    # SLDs are in 10^-6 inverse square angstroms. Intensity is in 1/cm.
    sld, solvent, scale, background = 2.8, 1.1, 0.7, 0.05
    return scale * 1e-4 * (sld - solvent)**2 * np.abs(amplitude)**2 / volume + background


def main():
    for shape in (prism, tetrahedron):
        for view in VIEWS:
            reference = intensity(*shape(28), view)
            check = intensity(*shape(16), view)
            np.testing.assert_allclose(check, reference, rtol=2e-12, atol=0.0)
            values = ", ".join(f"{value:.8f}" for value in reference)
            print(f"{shape.__name__} (theta, phi, psi) = {view}: [{values}]")


if __name__ == "__main__":
    main()
