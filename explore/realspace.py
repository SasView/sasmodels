from __future__ import division, print_function

import time
from copy import copy
import os
import argparse
from collections import OrderedDict

import numpy as np
from numpy import pi, radians, sin, cos, sqrt
from numpy.random import poisson, uniform, randn, rand, randint
from numpy.polynomial.legendre import leggauss
from scipy.integrate import simps
from scipy.special import j1 as J1

try:
    import numba
    USE_NUMBA = True
except ImportError:
    USE_NUMBA = False

# Definition of rotation matrices comes from wikipedia:
#    https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations
def Rx(angle):
    """Construct a matrix to rotate points about *x* by *angle* degrees."""
    a = radians(angle)
    R = [[1, 0, 0],
         [0, +cos(a), -sin(a)],
         [0, +sin(a), +cos(a)]]
    return np.matrix(R)

def Ry(angle):
    """Construct a matrix to rotate points about *y* by *angle* degrees."""
    a = radians(angle)
    R = [[+cos(a), 0, +sin(a)],
         [0, 1, 0],
         [-sin(a), 0, +cos(a)]]
    return np.matrix(R)

def Rz(angle):
    """Construct a matrix to rotate points about *z* by *angle* degrees."""
    a = radians(angle)
    R = [[+cos(a), -sin(a), 0],
         [+sin(a), +cos(a), 0],
         [0, 0, 1]]
    return np.matrix(R)

def rotation(theta, phi, psi):
    """
    Apply the jitter transform to a set of points.

    Points are stored in a 3 x n numpy matrix, not a numpy array or tuple.
    """
    return Rx(phi)*Ry(theta)*Rz(psi)

def apply_view(points, view):
    """
    Apply the view transform (theta, phi, psi) to a set of points.

    Points are stored in a 3 x n numpy array.

    View angles are in degrees.
    """
    theta, phi, psi = view
    return np.asarray((Rz(phi)*Ry(theta)*Rz(psi))*np.matrix(points.T)).T


def invert_view(qx, qy, view):
    """
    Return (qa, qb, qc) for the (theta, phi, psi) view angle at detector
    pixel (qx, qy).

    View angles are in degrees.
    """
    theta, phi, psi = view
    q = np.vstack((qx.flatten(), qy.flatten(), 0*qx.flatten()))
    return np.asarray((Rz(-psi)*Ry(-theta)*Rz(-phi))*np.matrix(q))


I3 = np.matrix([[1., 0, 0], [0, 1, 0], [0, 0, 1]])

class Shape:
    rotation = I3
    center = np.array([0., 0., 0.])[:, None]
    r_max = None
    lattice_size = np.array((1, 1, 1))
    lattice_spacing = np.array((1., 1., 1.))
    lattice_distortion = 0.0
    lattice_rotation = 0.0
    lattice_type = ""

    def volume(self):
        # type: () -> float
        raise NotImplementedError()

    def sample(self, density):
        # type: (float) -> np.ndarray[N], np.ndarray[N, 3]
        raise NotImplementedError()

    def dims(self):
        # type: () -> float, float, float
        raise NotImplementedError()

    def rotate(self, theta, phi, psi):
        if theta != 0. or phi != 0. or psi != 0.:
            self.rotation = rotation(theta, phi, psi) * self.rotation
        return self

    def shift(self, x, y, z):
        self.center = self.center + np.array([x, y, z])[:, None]
        return self

    def lattice(self, size=(1, 1, 1), spacing=(2, 2, 2), type="sc",
                distortion=0.0, rotation=0.0):
        self.lattice_size = np.asarray(size, 'i')
        self.lattice_spacing = np.asarray(spacing, 'd')
        self.lattice_type = type
        self.lattice_distortion = distortion
        self.lattice_rotation = rotation

    def _adjust(self, points):
        if self.rotation is I3:
            points = points.T + self.center
        else:
            points = np.asarray(self.rotation * np.matrix(points.T)) + self.center
        if self.lattice_type:
            points = self._apply_lattice(points)
        return points.T

    def r_bins(self, q, over_sampling=10, r_step=0.):
        if self.lattice_type:
            r_max = np.sqrt(np.sum(self.lattice_size*self.lattice_spacing*self.dims)**2)/2
        else:
            r_max = self.r_max
        #r_max = min(2 * pi / q[0], r_max)
        if r_step == 0.:
            r_step = 2 * pi / q[-1] / over_sampling
        #r_step = 0.01
        return np.arange(r_step, r_max, r_step)

    def _apply_lattice(self, points):
        """Spread points to different lattice positions"""
        size = self.lattice_size
        spacing = self.lattice_spacing
        shuffle = self.lattice_distortion
        rotate = self.lattice_rotation
        lattice = self.lattice_type

        if rotate != 0:
            # To vectorize the rotations we will need to unwrap the matrix multiply
            raise NotImplementedError("don't handle rotations yet")

        # Determine the number of lattice points in the lattice
        shapes_per_cell = 2 if lattice == "bcc" else 4 if lattice == "fcc" else 1
        number_of_lattice_points = np.prod(size) * shapes_per_cell

        # For each point in the original shape, figure out which lattice point
        # to translate it to.  This is both cell index (i*ny*nz + j*nz  + k) as
        # well as the point in the cell (corner, body center or face center).
        nsamples = points.shape[1]
        lattice_point = randint(number_of_lattice_points, size=nsamples)

        # Translate the cell index into the i,j,k coordinates of the senter
        cell_index = lattice_point // shapes_per_cell
        center = np.vstack((cell_index//(size[1]*size[2]),
                            (cell_index%(size[1]*size[2]))//size[2],
                            cell_index%size[2]))
        center = np.asarray(center, dtype='d')
        if lattice == "bcc":
            center[:, lattice_point % shapes_per_cell == 1] += [[0.5], [0.5], [0.5]]
        elif lattice == "fcc":
            center[:, lattice_point % shapes_per_cell == 1] += [[0.0], [0.5], [0.5]]
            center[:, lattice_point % shapes_per_cell == 2] += [[0.5], [0.0], [0.5]]
            center[:, lattice_point % shapes_per_cell == 3] += [[0.5], [0.5], [0.0]]

        # Each lattice point has its own displacement from the ideal position.
        # Not checking that shapes do not overlap if displacement is too large.
        offset = shuffle*(randn(3, number_of_lattice_points) if shuffle < 0.3
                          else rand(3, number_of_lattice_points))
        center += offset[:, cell_index]

        # Each lattice point has its own rotation.  Rotate the point prior to
        # applying any displacement.
        # rotation = rotate*(randn(size=(shapes, 3)) if shuffle < 30 else rand(size=(nsamples, 3)))
        # for k in shapes: points[k] = rotation[k]*points[k]
        points += center*(np.array([spacing])*np.array(self.dims)).T
        return points

class Composite(Shape):
    def __init__(self, shapes, center=(0, 0, 0), orientation=(0, 0, 0)):
        self.shapes = shapes
        self.rotate(*orientation)
        self.shift(*center)

        # Find the worst case distance between any two points amongst a set
        # of shapes independent of orientation.  This could easily be a
        # factor of two worse than necessary, e.g., a pair of thin rods
        # end-to-end vs the same pair side-by-side.
        distances = [((s1.r_max + s2.r_max)/2
                      + sqrt(np.sum((s1.center - s2.center)**2)))
                     for s1 in shapes
                     for s2 in shapes]
        self.r_max = max(distances + [s.r_max for s in shapes])
        self.volume = sum(shape.volume for shape in self.shapes)

    def sample(self, density):
        values, points = zip(*(shape.sample(density) for shape in self.shapes))
        return np.hstack(values), self._adjust(np.vstack(points))

class Box(Shape):
    def __init__(self, a, b, c,
                 value, center=(0, 0, 0), orientation=(0, 0, 0)):
        self.value = np.asarray(value)
        self.rotate(*orientation)
        self.shift(*center)
        self.a, self.b, self.c = a, b, c
        self._scale = np.array([a/2, b/2, c/2])[None, :]
        self.r_max = sqrt(a**2 + b**2 + c**2)
        self.dims = a, b, c
        self.volume = a*b*c

    def sample(self, density):
        num_points = poisson(density*self.volume)
        points = self._scale*uniform(-1, 1, size=(num_points, 3))
        values = self.value.repeat(points.shape[0])
        return values, self._adjust(points)

class EllipticalCylinder(Shape):
    def __init__(self, ra, rb, length,
                 value, center=(0, 0, 0), orientation=(0, 0, 0)):
        self.value = np.asarray(value)
        self.rotate(*orientation)
        self.shift(*center)
        self.ra, self.rb, self.length = ra, rb, length
        self._scale = np.array([ra, rb, length/2])[None, :]
        self.r_max = sqrt(4*max(ra, rb)**2 + length**2)
        self.dims = 2*ra, 2*rb, length
        self.volume = pi*ra*rb*length

    def sample(self, density):
        # randomly sample from a box of side length 2*r, excluding anything
        # not in the cylinder
        num_points = poisson(density*4*self.ra*self.rb*self.length)
        points = uniform(-1, 1, size=(num_points, 3))
        radius = points[:, 0]**2 + points[:, 1]**2
        points = points[radius <= 1]
        values = self.value.repeat(points.shape[0])
        return values, self._adjust(self._scale*points)

class EllipticalBicelle(Shape):
    def __init__(self, ra, rb, length,
                 thick_rim, thick_face,
                 value_core, value_rim, value_face,
                 center=(0, 0, 0), orientation=(0, 0, 0)):
        self.rotate(*orientation)
        self.shift(*center)
        self.value = value_core
        self.ra, self.rb, self.length = ra, rb, length
        self.thick_rim, self.thick_face = thick_rim, thick_face
        self.value_rim, self.value_face = value_rim, value_face

        # reset cylinder to outer dimensions for calculating scale, etc.
        ra = self.ra + self.thick_rim
        rb = self.rb + self.thick_rim
        length = self.length + 2*self.thick_face
        self._scale = np.array([ra, rb, length/2])[None, :]
        self.r_max = sqrt(4*max(ra, rb)**2 + length**2)
        self.dims = 2*ra, 2*rb, length
        self.volume = pi*ra*rb*length

    def sample(self, density):
        # randomly sample from a box of side length 2*r, excluding anything
        # not in the cylinder
        ra = self.ra + self.thick_rim
        rb = self.rb + self.thick_rim
        length = self.length + 2*self.thick_face
        num_points = poisson(density*4*ra*rb*length)
        points = uniform(-1, 1, size=(num_points, 3))
        radius = points[:, 0]**2 + points[:, 1]**2
        points = points[radius <= 1]
        # set all to core value first
        values = np.ones_like(points[:, 0])*self.value
        # then set value to face value if |z| > face/(length/2))
        values[abs(points[:, 2]) > self.length/(self.length + 2*self.thick_face)] = self.value_face
        # finally set value to rim value if outside the core ellipse
        radius = (points[:, 0]**2*(1 + self.thick_rim/self.ra)**2
                  + points[:, 1]**2*(1 + self.thick_rim/self.rb)**2)
        values[radius>1] = self.value_rim
        return values, self._adjust(self._scale*points)

class TriaxialEllipsoid(Shape):
    def __init__(self, ra, rb, rc,
                 value, center=(0, 0, 0), orientation=(0, 0, 0)):
        self.value = np.asarray(value)
        self.rotate(*orientation)
        self.shift(*center)
        self.ra, self.rb, self.rc = ra, rb, rc
        self._scale = np.array([ra, rb, rc])[None, :]
        self.r_max = 2*max(ra, rb, rc)
        self.dims = 2*ra, 2*rb, 2*rc
        self.volume = 4*pi/3 * ra * rb * rc

    def sample(self, density):
        # randomly sample from a box of side length 2*r, excluding anything
        # not in the ellipsoid
        num_points = poisson(density*8*self.ra*self.rb*self.rc)
        points = uniform(-1, 1, size=(num_points, 3))
        radius = np.sum(points**2, axis=1)
        points = self._scale*points[radius <= 1]
        values = self.value.repeat(points.shape[0])
        return values, self._adjust(points)

class Helix(Shape):
    def __init__(self, helix_radius, helix_pitch, tube_radius, tube_length,
                 value, center=(0, 0, 0), orientation=(0, 0, 0)):
        self.value = np.asarray(value)
        self.rotate(*orientation)
        self.shift(*center)
        helix_length = helix_pitch * tube_length/sqrt(helix_radius**2 + helix_pitch**2)
        total_radius = self.helix_radius + self.tube_radius
        self.helix_radius, self.helix_pitch = helix_radius, helix_pitch
        self.tube_radius, self.tube_length = tube_radius, tube_length
        self.r_max = sqrt(4*total_radius + (helix_length + 2*tube_radius)**2)
        self.dims = 2*total_radius, 2*total_radius, helix_length
        # small tube radius approximation; for larger tubes need to account
        # for the fact that the inner length is much shorter than the outer
        # length
        self.volume = pi*self.tube_radius**2*self.tube_length

    def points(self, density):
        num_points = poisson(density*4*self.tube_radius**2*self.tube_length)
        points = uniform(-1, 1, size=(num_points, 3))
        radius = points[:, 0]**2 + points[:, 1]**2
        points = points[radius <= 1]

        # Based on math stackexchange answer by Jyrki Lahtonen
        #     https://math.stackexchange.com/a/461637
        # with helix along z rather than x [so tuples in answer are (z, x, y)]
        # and with random points in the cross section (p1, p2) rather than
        # uniform points on the surface (cos u, sin u).
        a, R = self.tube_radius, self.helix_radius
        h = self.helix_pitch
        scale = 1/sqrt(R**2 + h**2)
        t = points[:, 3] * (self.tube_length * scale/2)
        cos_t, sin_t = cos(t), sin(t)

        # rx = R*cos_t
        # ry = R*sin_t
        # rz = h*t
        # nx = -a * cos_t * points[:, 1]
        # ny = -a * sin_t * points[:, 1]
        # nz = 0
        # bx = (a * h/scale) * sin_t * points[:, 2]
        # by = (-a * h/scale) * cos_t * points[:, 2]
        # bz = a*R/scale
        # x = rx + nx + bx
        # y = ry + ny + by
        # z = rz + nz + bz
        u, v = (R - a*points[:, 1]), (a * h/scale)*points[:, 2]
        x = u * cos_t + v * sin_t
        y = u * sin_t - v * cos_t
        z = a*R/scale + h * t

        points = np.hstack((x, y, z))
        values = self.value.repeat(points.shape[0])
        return values, self._adjust(points)

def csbox(a=10, b=20, c=30, da=1, db=2, dc=3, slda=1, sldb=2, sldc=3, sld_core=4):
    core = Box(a, b, c, sld_core)
    side_a = Box(da, b, c, slda, center=((a+da)/2, 0, 0))
    side_b = Box(a, db, c, sldb, center=(0, (b+db)/2, 0))
    side_c = Box(a, b, dc, sldc, center=(0, 0, (c+dc)/2))
    side_a2 = copy(side_a).shift(-a-da, 0, 0)
    side_b2 = copy(side_b).shift(0, -b-db, 0)
    side_c2 = copy(side_c).shift(0, 0, -c-dc)
    shape = Composite((core, side_a, side_b, side_c, side_a2, side_b2, side_c2))
    shape.dims = 2*da+a, 2*db+b, 2*dc+c
    return shape

def _Iqxy(values, x, y, z, qa, qb, qc):
    """I(q) = |sum V(r) rho(r) e^(1j q.r)|^2 / sum V(r)"""
    Iq = [abs(np.sum(values*np.exp(1j*(qa_k*x + qb_k*y + qc_k*z))))**2
            for qa_k, qb_k, qc_k in zip(qa.flat, qb.flat, qc.flat)]
    return Iq

if USE_NUMBA:
    # Override simple numpy solution with numba if available
    from numba import njit
    @njit("f8[:](f8[:],f8[:],f8[:],f8[:],f8[:],f8[:],f8[:])")
    def _Iqxy(values, x, y, z, qa, qb, qc):
        Iq = np.zeros_like(qa)
        for j in range(len(Iq)):
            total = 0. + 0j
            for k in range(len(values)):
                total += values[k]*np.exp(1j*(qa[j]*x[k] + qb[j]*y[k] + qc[j]*z[k]))
            Iq[j] = abs(total)**2
        return Iq

def calc_Iqxy(qx, qy, rho, points, volume=1.0, view=(0, 0, 0)):
    qx, qy = np.broadcast_arrays(qx, qy)
    qa, qb, qc = invert_view(qx, qy, view)
    rho, volume = np.broadcast_arrays(rho, volume)
    values = rho*volume
    x, y, z = points.T
    values, x, y, z, qa, qb, qc = [np.asarray(v, 'd')
                                   for v in (values, x, y, z, qa, qb, qc)]

    # I(q) = |sum V(r) rho(r) e^(1j q.r)|^2 / sum V(r)
    Iq = _Iqxy(values, x, y, z, qa.flatten(), qb.flatten(), qc.flatten())
    return np.asarray(Iq).reshape(qx.shape) / np.sum(volume)

def _calc_Pr_nonuniform(r, rho, points):
    # Make Pr a little be bigger than necessary so that only distances
    # min < d < max end up in Pr
    n_max = len(r)+1
    extended_Pr = np.zeros(n_max+1, 'd')
    # r refers to bin centers; find corresponding bin edges
    bins = bin_edges(r)
    t_next = time.clock() + 3
    for k, rho_k in enumerate(rho[:-1]):
        distance = sqrt(np.sum((points[k] - points[k+1:])**2, axis=1))
        weights = rho_k * rho[k+1:]
        index = np.searchsorted(bins, distance)
        # Note: indices may be duplicated, so "Pr[index] += w" will not work!!
        extended_Pr += np.bincount(index, weights, n_max+1)
        t = time.clock()
        if t > t_next:
            t_next = t + 3
            print("processing %d of %d"%(k, len(rho)-1))
    Pr = extended_Pr[1:-1]
    return Pr

def _calc_Pr_uniform(r, rho, points):
    # Make Pr a little be bigger than necessary so that only distances
    # min < d < max end up in Pr
    dr, n_max = r[0], len(r)
    extended_Pr = np.zeros(n_max+1, 'd')
    t0 = time.clock()
    t_next = t0 + 3
    for k, rho_k in enumerate(rho[:-1]):
        distances = sqrt(np.sum((points[k] - points[k+1:])**2, axis=1))
        weights = rho_k * rho[k+1:]
        index = np.minimum(np.asarray(distances/dr, 'i'), n_max)
        # Note: indices may be duplicated, so "Pr[index] += w" will not work!!
        extended_Pr += np.bincount(index, weights, n_max+1)
        t = time.clock()
        if t > t_next:
            t_next = t + 3
            print("processing %d of %d"%(k, len(rho)-1))
    #print("time py:", time.clock() - t0)
    Pr = extended_Pr[:-1]
    # Make Pr independent of sampling density.  The factor of 2 comes because
    # we are only accumulating the upper triangular distances.
    #Pr = Pr * 2 / len(rho)**2
    return Pr

    # Can get an additional 2x by going to C.  Cuda/OpenCL will allow even
    # more speedup, though still bounded by the n^2 cose.
    """
void pdfcalc(int n, const double *pts, const double *rho,
	     int nPr, double *Pr, double rstep)
{
  int i,j;

  for (i=0; i<n-2; i++) {
    for (j=i+1; j<=n-1; j++) {
      const double dxx=pts[3*i]-pts[3*j];
      const double dyy=pts[3*i+1]-pts[3*j+1];
      const double dzz=pts[3*i+2]-pts[3*j+2];
      const double d=sqrt(dxx*dxx+dyy*dyy+dzz*dzz);
      const int k=rint(d/rstep);
      if (k < nPr) Pr[k]+=rho[i]*rho[j];
    }
  }
}
"""

if USE_NUMBA:
    # Override simple numpy solution with numba if available
    @njit("f8[:](f8[:], f8[:], f8[:,:])")
    def _calc_Pr_uniform(r, rho, points):
        dr = r[0]
        n_max = len(r)
        Pr = np.zeros_like(r)
        for j in range(len(rho) - 1):
            x, y, z = points[j, 0], points[j, 1], points[j, 2]
            for k in range(j+1, len(rho)):
                distance = sqrt((x - points[k, 0])**2
                                + (y - points[k, 1])**2
                                + (z - points[k, 2])**2)
                index = int(distance/dr)
                if index < n_max:
                    Pr[index] += rho[j] * rho[k]
        # Make Pr independent of sampling density.  The factor of 2 comes because
        # we are only accumulating the upper triangular distances.
        #Pr = Pr * 2 / len(rho)**2
        return Pr


def calc_Pr(r, rho, points):
    # P(r) with uniform steps in r is 3x faster; check if we are uniform
    # before continuing
    r, rho, points = [np.asarray(v, 'd') for v in (r, rho, points)]
    if np.max(np.abs(np.diff(r) - r[0])) > r[0]*0.01:
        Pr = _calc_Pr_nonuniform(r, rho, points)
    else:
        Pr = _calc_Pr_uniform(r, rho, points)
    return Pr / Pr.max()


def j0(x):
    return np.sinc(x/np.pi)

def calc_Iq(q, r, Pr):
    Iq = np.array([simps(Pr * j0(qk*r), r) for qk in q])
    #Iq = np.array([np.trapz(Pr * j0(qk*r), r) for qk in q])
    Iq /= Iq[0]
    return Iq

# NOTE: copied from sasmodels/resolution.py
def bin_edges(x):
    """
    Determine bin edges from bin centers, assuming that edges are centered
    between the bins.

    Note: this uses the arithmetic mean, which may not be appropriate for
    log-scaled data.
    """
    if len(x) < 2 or (np.diff(x) < 0).any():
        raise ValueError("Expected bins to be an increasing set")
    edges = np.hstack([
        x[0]  - 0.5*(x[1]  - x[0]),  # first point minus half first interval
        0.5*(x[1:] + x[:-1]),        # mid points of all central intervals
        x[-1] + 0.5*(x[-1] - x[-2]), # last point plus half last interval
        ])
    return edges

# -------------- plotters ----------------
def plot_calc(r, Pr, q, Iq, theory=None, title=None):
    import matplotlib.pyplot as plt
    plt.subplot(211)
    plt.plot(r, Pr, '-', label="Pr")
    plt.xlabel('r (A)')
    plt.ylabel('Pr (1/A^2)')
    if title is not None:
        plt.title(title)
    plt.subplot(212)
    plt.loglog(q, Iq, '-', label='from Pr')
    plt.xlabel('q (1/A')
    plt.ylabel('Iq')
    if theory is not None:
        plt.loglog(theory[0], theory[1]/theory[1][0], '-', label='analytic')
        plt.legend()

def plot_calc_2d(qx, qy, Iqxy, theory=None, title=None):
    import matplotlib.pyplot as plt
    qx, qy = bin_edges(qx), bin_edges(qy)
    #qx, qy = np.meshgrid(qx, qy)
    if theory is not None:
        plt.subplot(121)
    #plt.pcolor(qx, qy, np.log10(Iqxy))
    extent = [qx[0], qx[-1], qy[0], qy[-1]]
    plt.imshow(np.log10(Iqxy), extent=extent, interpolation="nearest",
               origin='lower')
    plt.xlabel('qx (1/A)')
    plt.ylabel('qy (1/A)')
    plt.axis('equal')
    plt.axis(extent)
    #plt.grid(True)
    if title is not None:
        plt.title(title)
    if theory is not None:
        plt.subplot(122)
        plt.imshow(np.log10(theory), extent=extent, interpolation="nearest",
                   origin='lower')
        plt.axis('equal')
        plt.axis(extent)
        plt.xlabel('qx (1/A)')

def plot_points(rho, points):
    import mpl_toolkits.mplot3d
    import matplotlib.pyplot as plt

    ax = plt.axes(projection='3d')
    try:
        ax.axis('square')
    except Exception:
        pass
    n = len(points)
    #print("len points", n)
    index = np.random.choice(n, size=500) if n > 500 else slice(None, None)
    ax.scatter(points[index, 0], points[index, 1], points[index, 2], c=rho[index])
    # make square axes
    minmax = np.array([points.min(), points.max()])
    ax.scatter(minmax, minmax, minmax, c='w')
    #low, high = points.min(axis=0), points.max(axis=0)
    #ax.axis([low[0], high[0], low[1], high[1], low[2], high[2]])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.autoscale(True)

# ----------- Analytic models --------------
def sas_sinx_x(x):
    with np.errstate(all='ignore'):
        retvalue = sin(x)/x
    retvalue[x == 0.] = 1.
    return retvalue

def sas_2J1x_x(x):
    with np.errstate(all='ignore'):
        retvalue = 2*J1(x)/x
    retvalue[x == 0] = 1.
    return retvalue

def sas_3j1x_x(x):
    """return 3*j1(x)/x"""
    with np.errstate(all='ignore'):
        retvalue = 3*(sin(x) - x*cos(x))/x**3
    retvalue[x == 0.] = 1.
    return retvalue

def cylinder_Iq(q, radius, length):
    z, w = leggauss(76)
    cos_alpha = (z+1)/2
    sin_alpha = sqrt(1.0 - cos_alpha**2)
    Iq = np.empty_like(q)
    for k, qk in enumerate(q):
        qab, qc = qk*sin_alpha, qk*cos_alpha
        Fq = sas_2J1x_x(qab*radius) * sas_sinx_x(qc*length/2)
        Iq[k] = np.sum(w*Fq**2)
    Iq = Iq
    return Iq

def cylinder_Iqxy(qx, qy, radius, length, view=(0, 0, 0)):
    qa, qb, qc = invert_view(qx, qy, view)
    qab = sqrt(qa**2 + qb**2)
    Fq = sas_2J1x_x(qab*radius) * sas_sinx_x(qc*length/2)
    Iq = Fq**2
    return Iq.reshape(qx.shape)

def sphere_Iq(q, radius):
    Iq = sas_3j1x_x(q*radius)**2
    return Iq

def box_Iq(q, a, b, c):
    z, w = leggauss(76)
    outer_sum = np.zeros_like(q)
    for cos_alpha, outer_w in zip((z+1)/2, w):
        sin_alpha = sqrt(1.0-cos_alpha*cos_alpha)
        qc = q*cos_alpha
        siC = c*sas_sinx_x(c*qc/2)
        inner_sum = np.zeros_like(q)
        for beta, inner_w in zip((z + 1)*pi/4, w):
            qa, qb = q*sin_alpha*sin(beta), q*sin_alpha*cos(beta)
            siA = a*sas_sinx_x(a*qa/2)
            siB = b*sas_sinx_x(b*qb/2)
            Fq = siA*siB*siC
            inner_sum += inner_w * Fq**2
        outer_sum += outer_w * inner_sum
    Iq = outer_sum / 4  # = outer*um*zm*8.0/(4.0*M_PI)
    return Iq

def box_Iqxy(qx, qy, a, b, c, view=(0, 0, 0)):
    qa, qb, qc = invert_view(qx, qy, view)
    sia = sas_sinx_x(qa*a/2)
    sib = sas_sinx_x(qb*b/2)
    sic = sas_sinx_x(qc*c/2)
    Fq = sia*sib*sic
    Iq = Fq**2
    return Iq.reshape(qx.shape)

def csbox_Iq(q, a, b, c, da, db, dc, slda, sldb, sldc, sld_core):
    z, w = leggauss(76)

    sld_solvent = 0
    overlapping = False
    dr0 = sld_core - sld_solvent
    drA, drB, drC = slda-sld_solvent, sldb-sld_solvent, sldc-sld_solvent
    tA, tB, tC = a + 2*da, b + 2*db, c + 2*dc

    outer_sum = np.zeros_like(q)
    for cos_alpha, outer_w in zip((z+1)/2, w):
        sin_alpha = sqrt(1.0-cos_alpha*cos_alpha)
        qc = q*cos_alpha
        siC = c*sas_sinx_x(c*qc/2)
        siCt = tC*sas_sinx_x(tC*qc/2)
        inner_sum = np.zeros_like(q)
        for beta, inner_w in zip((z + 1)*pi/4, w):
            qa, qb = q*sin_alpha*sin(beta), q*sin_alpha*cos(beta)
            siA = a*sas_sinx_x(a*qa/2)
            siB = b*sas_sinx_x(b*qb/2)
            siAt = tA*sas_sinx_x(tA*qa/2)
            siBt = tB*sas_sinx_x(tB*qb/2)
            if overlapping:
                Fq = (dr0*siA*siB*siC
                      + drA*(siAt-siA)*siB*siC
                      + drB*siAt*(siBt-siB)*siC
                      + drC*siAt*siBt*(siCt-siC))
            else:
                Fq = (dr0*siA*siB*siC
                      + drA*(siAt-siA)*siB*siC
                      + drB*siA*(siBt-siB)*siC
                      + drC*siA*siB*(siCt-siC))
            inner_sum += inner_w * Fq**2
        outer_sum += outer_w * inner_sum
    Iq = outer_sum / 4  # = outer*um*zm*8.0/(4.0*M_PI)
    return Iq/Iq[0]

def csbox_Iqxy(qx, qy, a, b, c, da, db, dc, slda, sldb, sldc, sld_core, view=(0,0,0)):
    qa, qb, qc = invert_view(qx, qy, view)

    sld_solvent = 0
    overlapping = False
    dr0 = sld_core - sld_solvent
    drA, drB, drC = slda-sld_solvent, sldb-sld_solvent, sldc-sld_solvent
    tA, tB, tC = a + 2*da, b + 2*db, c + 2*dc
    siA = a*sas_sinx_x(a*qa/2)
    siB = b*sas_sinx_x(b*qb/2)
    siC = c*sas_sinx_x(c*qc/2)
    siAt = tA*sas_sinx_x(tA*qa/2)
    siBt = tB*sas_sinx_x(tB*qb/2)
    siCt = tC*sas_sinx_x(tC*qc/2)
    Fq = (dr0*siA*siB*siC
          + drA*(siAt-siA)*siB*siC
          + drB*siA*(siBt-siB)*siC
          + drC*siA*siB*(siCt-siC))
    Iq = Fq**2
    return Iq.reshape(qx.shape)

def sasmodels_Iq(kernel, q, pars):
    from sasmodels.data import empty_data1D
    from sasmodels.direct_model import DirectModel
    data = empty_data1D(q)
    calculator = DirectModel(data, kernel)
    Iq = calculator(**pars)
    return Iq

def sasmodels_Iqxy(kernel, qx, qy, pars, view):
    from sasmodels.data import Data2D
    from sasmodels.direct_model import DirectModel
    Iq = 100 * np.ones_like(qx)
    data = Data2D(x=qx, y=qy, z=Iq, dx=None, dy=None, dz=np.sqrt(Iq))
    data.x_bins = qx[0, :]
    data.y_bins = qy[:, 0]
    data.filename = "fake data"

    calculator = DirectModel(data, kernel)
    pars_plus_view = pars.copy()
    pars_plus_view.update(theta=view[0], phi=view[1], psi=view[2])
    Iqxy = calculator(**pars_plus_view)
    return Iqxy.reshape(qx.shape)

def wrap_sasmodel(name, **pars):
    from sasmodels.core import load_model
    kernel = load_model(name)
    fn = lambda q: sasmodels_Iq(kernel, q, pars)
    fn_xy = lambda qx, qy, view: sasmodels_Iqxy(kernel, qx, qy, pars, view)
    return fn, fn_xy


# --------- Test cases -----------

def build_cylinder(radius=25, length=125, rho=2.):
    shape = EllipticalCylinder(radius, radius, length, rho)
    fn = lambda q: cylinder_Iq(q, radius, length)*rho**2
    fn_xy = lambda qx, qy, view: cylinder_Iqxy(qx, qy, radius, length, view=view)*rho**2
    return shape, fn, fn_xy

DEFAULT_SPHERE_RADIUS = 125
DEFAULT_SPHERE_CONTRAST = 2
def build_sphere(radius=DEFAULT_SPHERE_RADIUS, rho=DEFAULT_SPHERE_CONTRAST):
    shape = TriaxialEllipsoid(radius, radius, radius, rho)
    fn = lambda q: sphere_Iq(q, radius)*rho**2
    fn_xy = lambda qx, qy, view: sphere_Iq(np.sqrt(qx**2+qy**2), radius)*rho**2
    return shape, fn, fn_xy

def build_box(a=10, b=20, c=30, rho=2.):
    shape = Box(a, b, c, rho)
    fn = lambda q: box_Iq(q, a, b, c)*rho**2
    fn_xy = lambda qx, qy, view: box_Iqxy(qx, qy, a, b, c, view=view)*rho**2
    return shape, fn, fn_xy

def build_csbox(a=10, b=20, c=30, da=1, db=2, dc=3, slda=1, sldb=2, sldc=3, sld_core=4):
    shape = csbox(a, b, c, da, db, dc, slda, sldb, sldc, sld_core)
    fn = lambda q: csbox_Iq(q, a, b, c, da, db, dc, slda, sldb, sldc, sld_core)
    fn_xy = lambda qx, qy, view: csbox_Iqxy(qx, qy, a, b, c, da, db, dc,
                                            slda, sldb, sldc, sld_core, view=view)
    return shape, fn, fn_xy

def build_ellcyl(ra=25, rb=50, length=125, rho=2.):
    shape = EllipticalCylinder(ra, rb, length, rho)
    fn, fn_xy = wrap_sasmodel(
        'elliptical_cylinder',
        scale=1,
        background=0,
        radius_minor=ra,
        axis_ratio=rb/ra,
        length=length,
        sld=rho,
        sld_solvent=0,
    )
    return shape, fn, fn_xy

def build_cscyl(ra=30, rb=90, length=30, thick_rim=8, thick_face=14,
                sld_core=4, sld_rim=1, sld_face=7):
    shape = EllipticalBicelle(
        ra=ra, rb=rb, length=length,
        thick_rim=thick_rim, thick_face=thick_face,
        value_core=sld_core, value_rim=sld_rim, value_face=sld_face,
        )
    fn, fn_xy = wrap_sasmodel(
        'core_shell_bicelle_elliptical',
        scale=1,
        background=0,
        radius=ra,
        x_core=rb/ra,
        length=length,
        thick_rim=thick_rim,
        thick_face=thick_face,
        sld_core=sld_core,
        sld_face=sld_face,
        sld_rim=sld_rim,
        sld_solvent=0,
    )
    return shape, fn, fn_xy

def build_sc_lattice(shape, nx=1, ny=1, nz=1, dx=2, dy=2, dz=2,
                        shuffle=0, rotate=0):
    a, b, c = shape.dims
    corners= [copy(shape)
              .shift((ix+(randn() if shuffle < 0.3 else rand())*shuffle)*dx*a,
                     (iy+(randn() if shuffle < 0.3 else rand())*shuffle)*dy*b,
                     (iz+(randn() if shuffle < 0.3 else rand())*shuffle)*dz*c)
              .rotate(*((randn(3) if rotate < 30 else rand(3))*rotate))
              for ix in range(nx)
              for iy in range(ny)
              for iz in range(nz)]
    lattice = Composite(corners)
    return lattice

def build_bcc_lattice(shape, nx=1, ny=1, nz=1, dx=2, dy=2, dz=2,
                      shuffle=0, rotate=0):
    a, b, c = shape.dims
    corners = [copy(shape)
               .shift((ix+(randn() if shuffle < 0.3 else rand())*shuffle)*dx*a,
                      (iy+(randn() if shuffle < 0.3 else rand())*shuffle)*dy*b,
                      (iz+(randn() if shuffle < 0.3 else rand())*shuffle)*dz*c)
               .rotate(*((randn(3) if rotate < 30 else rand(3))*rotate))
               for ix in range(nx)
               for iy in range(ny)
               for iz in range(nz)]
    centers = [copy(shape)
               .shift((ix+0.5+(randn() if shuffle < 0.3 else rand())*shuffle)*dx*a,
                      (iy+0.5+(randn() if shuffle < 0.3 else rand())*shuffle)*dy*b,
                      (iz+0.5+(randn() if shuffle < 0.3 else rand())*shuffle)*dz*c)
               .rotate(*((randn(3) if rotate < 30 else rand(3))*rotate))
               for ix in range(nx)
               for iy in range(ny)
               for iz in range(nz)]
    lattice = Composite(corners + centers)
    return lattice

def build_fcc_lattice(shape, nx=1, ny=1, nz=1, dx=2, dy=2, dz=2,
                      shuffle=0, rotate=0):
    a, b, c = shape.dims
    corners = [copy(shape)
               .shift((ix+(randn() if shuffle < 0.3 else rand())*shuffle)*dx*a,
                      (iy+(randn() if shuffle < 0.3 else rand())*shuffle)*dy*b,
                      (iz+(randn() if shuffle < 0.3 else rand())*shuffle)*dz*c)
               .rotate(*((randn(3) if rotate < 30 else rand(3))*rotate))
               for ix in range(nx)
               for iy in range(ny)
               for iz in range(nz)]
    faces_a = [copy(shape)
               .shift((ix+0.0+(randn() if shuffle < 0.3 else rand())*shuffle)*dx*a,
                      (iy+0.5+(randn() if shuffle < 0.3 else rand())*shuffle)*dy*b,
                      (iz+0.5+(randn() if shuffle < 0.3 else rand())*shuffle)*dz*c)
               .rotate(*((randn(3) if rotate < 30 else rand(3))*rotate))
               for ix in range(nx)
               for iy in range(ny)
               for iz in range(nz)]
    faces_b = [copy(shape)
               .shift((ix+0.5+(randn() if shuffle < 0.3 else rand())*shuffle)*dx*a,
                      (iy+0.0+(randn() if shuffle < 0.3 else rand())*shuffle)*dy*b,
                      (iz+0.5+(randn() if shuffle < 0.3 else rand())*shuffle)*dz*c)
               .rotate(*((randn(3) if rotate < 30 else rand(3))*rotate))
               for ix in range(nx)
               for iy in range(ny)
               for iz in range(nz)]
    faces_c = [copy(shape)
               .shift((ix+0.5+(randn() if shuffle < 0.3 else rand())*shuffle)*dx*a,
                      (iy+0.5+(randn() if shuffle < 0.3 else rand())*shuffle)*dy*b,
                      (iz+0.0+(randn() if shuffle < 0.3 else rand())*shuffle)*dz*c)
               .rotate(*((randn(3) if rotate < 30 else rand(3))*rotate))
               for ix in range(nx)
               for iy in range(ny)
               for iz in range(nz)]
    lattice = Composite(corners + faces_a + faces_b + faces_c)
    return lattice

SHAPE_FUNCTIONS = OrderedDict([
    ("cyl", build_cylinder),
    ("ellcyl", build_ellcyl),
    ("sphere", build_sphere),
    ("box", build_box),
    ("csbox", build_csbox),
    ("cscyl", build_cscyl),
])
SHAPES = list(SHAPE_FUNCTIONS.keys())
LATTICE_FUNCTIONS = OrderedDict([
    ("sc", build_sc_lattice),
    ("bcc", build_bcc_lattice),
    ("fcc", build_fcc_lattice),
])
LATTICE_TYPES = list(LATTICE_FUNCTIONS.keys())

def check_shape(title, shape, fn=None, show_points=False,
                mesh=100, qmax=1.0, r_step=0.01, samples=5000):
    rho_solvent = 0
    qmin = qmax/100.
    q = np.logspace(np.log10(qmin), np.log10(qmax), mesh)
    r = shape.r_bins(q, r_step=r_step)
    sampling_density = samples / shape.volume
    print("sampling points")
    rho, points = shape.sample(sampling_density)
    print("calculating Pr")
    t0 = time.time()
    Pr = calc_Pr(r, rho-rho_solvent, points)
    print("calc Pr time", time.time() - t0)
    Iq = calc_Iq(q, r, Pr)
    theory = (q, fn(q)) if fn is not None else None

    import pylab
    if show_points:
        plot_points(rho, points); pylab.figure()
    plot_calc(r, Pr, q, Iq, theory=theory, title=title)
    pylab.gcf().canvas.set_window_title(title)
    pylab.show()

def check_shape_2d(title, shape, fn=None, view=(0, 0, 0), show_points=False,
                   mesh=100, qmax=1.0, samples=5000):
    rho_solvent = 0
    #qx = np.linspace(0.0, qmax, mesh)
    #qy = np.linspace(0.0, qmax, mesh)
    qx = np.linspace(-qmax, qmax, mesh)
    qy = np.linspace(-qmax, qmax, mesh)
    Qx, Qy = np.meshgrid(qx, qy)
    sampling_density = samples / shape.volume
    print("sampling points")
    t0 = time.time()
    rho, points = shape.sample(sampling_density)
    print("point generation time", time.time() - t0)
    t0 = time.time()
    Iqxy = calc_Iqxy(Qx, Qy, rho, points, view=view)
    print("calc Iqxy time", time.time() - t0)
    t0 = time.time()
    theory = fn(Qx, Qy, view) if fn is not None else None
    print("calc theory time", time.time() - t0)
    Iqxy += 0.001 * Iqxy.max()
    if theory is not None:
        theory += 0.001 * theory.max()

    import pylab
    if show_points:
        plot_points(rho, points); pylab.figure()
    plot_calc_2d(qx, qy, Iqxy, theory=theory, title=title)
    pylab.gcf().canvas.set_window_title(title)
    pylab.show()

def main():
    parser = argparse.ArgumentParser(
        description="Compute scattering from realspace sampling",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
    parser.add_argument('-d', '--dim', type=int, default=1,
                        help='dimension 1 or 2')
    parser.add_argument('-m', '--mesh', type=int, default=100,
                        help='number of mesh points')
    parser.add_argument('-s', '--samples', type=int, default=5000,
                        help="number of sample points")
    parser.add_argument('-q', '--qmax', type=float, default=0.5,
                        help='max q')
    parser.add_argument('-v', '--view', type=str, default='0,0,0',
                        help='theta,phi,psi angles')
    parser.add_argument('-n', '--lattice', type=str, default='1,1,1',
                        help='lattice size')
    parser.add_argument('-z', '--spacing', type=str, default='2,2,2',
                        help='lattice spacing (relative to shape)')
    parser.add_argument('-t', '--type', choices=LATTICE_TYPES,
                        default=LATTICE_TYPES[0],
                        help='lattice type')
    parser.add_argument('-r', '--rotate', type=float, default=0.,
                        help="rotation relative to lattice, gaussian < 30 degrees, uniform otherwise")
    parser.add_argument('-w', '--shuffle', type=float, default=0.,
                        help="position relative to lattice, gaussian < 0.3, uniform otherwise")
    parser.add_argument('-p', '--plot', action='store_true',
                        help='plot points')
    parser.add_argument('shape', choices=SHAPES, nargs='?', default=SHAPES[0],
                        help='oriented shape')
    parser.add_argument('pars', type=str, nargs='*', help='shape parameters')
    opts = parser.parse_args()
    pars = {key: float(value) for p in opts.pars for key, value in [p.split('=')]}
    nx, ny, nz = [int(v) for v in opts.lattice.split(',')]
    dx, dy, dz = [float(v) for v in opts.spacing.split(',')]
    distortion, rotation = opts.shuffle, opts.rotate
    shape, fn, fn_xy = SHAPE_FUNCTIONS[opts.shape](**pars)
    view = tuple(float(v) for v in opts.view.split(','))
    # If comparing a sphere in a cubic lattice, compare against the
    # corresponding paracrystalline model.
    if opts.shape == "sphere" and dx == dy == dz and nx*ny*nz > 1:
        radius = pars.get('radius', DEFAULT_SPHERE_RADIUS)
        model_name = opts.type + "_paracrystal"
        model_pars = {
            "scale": 1.,
            "background": 0.,
            "lattice_spacing": 2*radius*dx,
            "lattice_distortion": distortion,
            "radius": radius,
            "sld": pars.get('rho', DEFAULT_SPHERE_CONTRAST),
            "sld_solvent": 0.,
            "theta": view[0],
            "phi": view[1],
            "psi": view[2],
        }
        fn, fn_xy = wrap_sasmodel(model_name, **model_pars)
    if nx*ny*nz > 1:
        if rotation != 0:
            print("building %s lattice"%opts.type)
            build_lattice = LATTICE_FUNCTIONS[opts.type]
            shape = build_lattice(shape, nx, ny, nz, dx, dy, dz,
                                  distortion, rotation)
        else:
            shape.lattice(size=(nx, ny, nz), spacing=(dx, dy, dz),
                          type=opts.type,
                          rotation=rotation, distortion=distortion)

    title = "%s(%s)" % (opts.shape, " ".join(opts.pars))
    if opts.dim == 1:
        check_shape(title, shape, fn, show_points=opts.plot,
                    mesh=opts.mesh, qmax=opts.qmax, samples=opts.samples)
    else:
        check_shape_2d(title, shape, fn_xy, view=view, show_points=opts.plot,
                       mesh=opts.mesh, qmax=opts.qmax, samples=opts.samples)


if __name__ == "__main__":
    main()
