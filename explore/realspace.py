from __future__ import division, print_function

import time
from copy import copy

import numpy as np
from numpy import pi, radians, sin, cos, sqrt
from numpy.random import poisson, uniform
from numpy.polynomial.legendre import leggauss
from scipy.integrate import simps
from scipy.special import j1 as J1

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


class Shape:
    rotation = np.matrix([[1., 0, 0], [0, 1, 0], [0, 0, 1]])
    center = np.array([0., 0., 0.])[:, None]
    r_max = None

    def volume(self):
        # type: () -> float
        raise NotImplementedError()

    def sample(self, density):
        # type: (float) -> np.ndarray[N], np.ndarray[N, 3]
        raise NotImplementedError()

    def rotate(self, theta, phi, psi):
        self.rotation = rotation(theta, phi, psi) * self.rotation
        return self

    def shift(self, x, y, z):
        self.center = self.center + np.array([x, y, z])[:, None]
        return self

    def _adjust(self, points):
        points = np.asarray(self.rotation * np.matrix(points.T)) + self.center
        return points.T

    def r_bins(self, q, over_sampling=1, r_step=0.):
        r_max = min(2 * pi / q[0], self.r_max)
        if r_step == 0.:
            r_step = 2 * pi / q[-1] / over_sampling
        #r_step = 0.01
        return np.arange(r_step, r_max, r_step)

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

    def volume(self):
        return sum(shape.volume() for shape in self.shapes)

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

    def volume(self):
        return self.a*self.b*self.c

    def sample(self, density):
        num_points = poisson(density*self.a*self.b*self.c)
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

    def volume(self):
        return pi*self.ra*self.rb*self.length

    def sample(self, density):
        # density of the bounding box
        num_points = poisson(density*4*self.ra*self.rb*self.length)
        points = uniform(-1, 1, size=(num_points, 3))
        radius = points[:, 0]**2 + points[:, 1]**2
        points = self._scale*points[radius <= 1]
        values = self.value.repeat(points.shape[0])
        return values, self._adjust(points)

class TriaxialEllipsoid(Shape):
    def __init__(self, ra, rb, rc,
                 value, center=(0, 0, 0), orientation=(0, 0, 0)):
        self.value = np.asarray(value)
        self.rotate(*orientation)
        self.shift(*center)
        self.ra, self.rb, self.rc = ra, rb, rc
        self._scale = np.array([ra, rb, rc])[None, :]
        self.r_max = 2*max(ra, rb, rc)

    def volume(self):
        return 4*pi/3 * self.ra * self.rb * self.rc

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
        self.helix_radius, self.helix_pitch = helix_radius, helix_pitch
        self.tube_radius, self.tube_length = tube_radius, tube_length
        helix_length = helix_pitch * tube_length/sqrt(helix_radius**2 + helix_pitch**2)
        self.r_max = sqrt((2*helix_radius + 2*tube_radius)*2
                          + (helix_length + 2*tube_radius)**2)

    def volume(self):
        # small tube radius approximation; for larger tubes need to account
        # for the fact that the inner length is much shorter than the outer
        # length
        return pi*self.tube_radius**2*self.tube_length

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

NUMBA = False
if NUMBA:
    from numba import njit
    @njit("f8[:](f8[:],f8[:],f8[:],f8[:],f8[:],f8[:],f8[:])")
    def _Iqxy(values, x, y, z, qa, qb, qc):
        Iq = np.zeros_like(qa)
        for j in range(len(Iq)):
            total = 0. + 0j
            for k in range(len(Iq)):
                total += values[k]*np.exp(1j*(qa[j]*x[k] + qb[j]*y[k] + qc[j]*z[k]))
            Iq[j] = abs(total)**2
        return Iq

def calc_Iqxy(qx, qy, rho, points, volume=1.0, view=(0, 0, 0)):
    qx, qy = np.broadcast_arrays(qx, qy)
    qa, qb, qc = invert_view(qx, qy, view)
    rho, volume = np.broadcast_arrays(rho, volume)
    values = rho*volume
    x, y, z = points.T

    # I(q) = |sum V(r) rho(r) e^(1j q.r)|^2 / sum V(r)
    if NUMBA:
        Iq = _Iqxy(values, x, y, z, qa.flatten(), qb.flatten(), qc.flatten())
    else:
        Iq = [abs(np.sum(values*np.exp(1j*(qa_k*x + qb_k*y + qc_k*z))))**2
              for qa_k, qb_k, qc_k in zip(qa.flat, qb.flat, qc.flat)]
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

if NUMBA:
    @njit("f8[:](f8[:], f8[:], f8[:,:])")
    def _calc_Pr_uniform_jit(r, rho, points):
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
    if np.max(np.abs(np.diff(r) - r[0])) > r[0]*0.01:
        Pr = _calc_Pr_nonuniform(r, rho, points)
    else:
        if NUMBA:
            Pr = _calc_Pr_uniform_jit(r, rho, points)
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

def plot_calc(r, Pr, q, Iq, theory=None):
    import matplotlib.pyplot as plt
    plt.subplot(211)
    plt.plot(r, Pr, '-', label="Pr")
    plt.xlabel('r (A)')
    plt.ylabel('Pr (1/A^2)')
    plt.subplot(212)
    plt.loglog(q, Iq, '-', label='from Pr')
    plt.xlabel('q (1/A')
    plt.ylabel('Iq')
    if theory is not None:
        plt.loglog(theory[0], theory[1], '-', label='analytic')
        plt.legend()

def plot_calc_2d(qx, qy, Iqxy, theory=None):
    import matplotlib.pyplot as plt
    qx, qy = bin_edges(qx), bin_edges(qy)
    #qx, qy = np.meshgrid(qx, qy)
    if theory is not None:
        plt.subplot(121)
    plt.pcolormesh(qx, qy, np.log10(Iqxy))
    plt.xlabel('qx (1/A)')
    plt.ylabel('qy (1/A)')
    if theory is not None:
        plt.subplot(122)
        plt.pcolormesh(qx, qy, np.log10(theory))
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
    #low, high = points.min(axis=0), points.max(axis=0)
    #ax.axis([low[0], high[0], low[1], high[1], low[2], high[2]])
    ax.autoscale(True)

def check_shape(shape, fn=None):
    rho_solvent = 0
    q = np.logspace(-3, 0, 200)
    r = shape.r_bins(q, r_step=0.01)
    sampling_density = 6*5000 / shape.volume()
    rho, points = shape.sample(sampling_density)
    t0 = time.time()
    Pr = calc_Pr(r, rho-rho_solvent, points)
    print("calc Pr time", time.time() - t0)
    Iq = calc_Iq(q, r, Pr)
    theory = (q, fn(q)) if fn is not None else None

    import pylab
    #plot_points(rho, points); pylab.figure()
    plot_calc(r, Pr, q, Iq, theory=theory)
    pylab.show()

def check_shape_2d(shape, fn=None, view=(0, 0, 0)):
    rho_solvent = 0
    nq, qmax = 100, 1.0
    qx = np.linspace(0.0, qmax, nq)
    qy = np.linspace(0.0, qmax, nq)
    Qx, Qy = np.meshgrid(qx, qy)
    sampling_density = 50000 / shape.volume()
    #t0 = time.time()
    rho, points = shape.sample(sampling_density)
    #print("sample time", time.time() - t0)
    t0 = time.time()
    Iqxy = calc_Iqxy(Qx, Qy, rho, points, view=view)
    print("calc time", time.time() - t0)
    theory = fn(Qx, Qy) if fn is not None else None
    Iqxy += 0.001 * Iqxy.max()
    if theory is not None:
        theory += 0.001 * theory.max()

    import pylab
    #plot_points(rho, points); pylab.figure()
    plot_calc_2d(qx, qy, Iqxy, theory=theory)
    pylab.show()

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
        Fq = sas_2J1x_x(qab*radius) * j0(qc*length/2)
        Iq[k] = np.sum(w*Fq**2)
    Iq = Iq/Iq[0]
    return Iq

def cylinder_Iqxy(qx, qy, radius, length, view=(0, 0, 0)):
    qa, qb, qc = invert_view(qx, qy, view)
    qab = np.sqrt(qa**2 + qb**2)
    Fq = sas_2J1x_x(qab*radius) * j0(qc*length/2)
    Iq = Fq**2
    return Iq.reshape(qx.shape)

def sphere_Iq(q, radius):
    Iq = sas_3j1x_x(q*radius)**2
    return Iq/Iq[0]

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
        siC = c*j0(c*qc/2)
        siCt = tC*j0(tC*qc/2)
        inner_sum = np.zeros_like(q)
        for beta, inner_w in zip((z + 1)*pi/4, w):
            qa, qb = q*sin_alpha*sin(beta), q*sin_alpha*cos(beta)
            siA = a*j0(a*qa/2)
            siB = b*j0(b*qb/2)
            siAt = tA*j0(tA*qa/2)
            siBt = tB*j0(tB*qb/2)
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

def check_cylinder(radius=25, length=125, rho=2.):
    shape = EllipticalCylinder(radius, radius, length, rho)
    fn = lambda q: cylinder_Iq(q, radius, length)
    check_shape(shape, fn)

def check_cylinder_2d(radius=25, length=125, rho=2., view=(0, 0, 0)):
    shape = EllipticalCylinder(radius, radius, length, rho)
    fn = lambda qx, qy, view=view: cylinder_Iqxy(qx, qy, radius, length, view=view)
    check_shape_2d(shape, fn, view=view)

def check_cylinder_2d_lattice(radius=25, length=125, rho=2.,
                              view=(0, 0, 0)):
    nx, dx = 1, 2*radius
    ny, dy = 30, 2*radius
    nz, dz = 30, length
    dx, dy, dz = 2*dx, 2*dy, 2*dz
    def center(*args):
        sigma = 0.333
        space = 2
        return [(space*n+np.random.randn()*sigma)*x for n, x in args]
    shapes = [EllipticalCylinder(radius, radius, length, rho,
                                 #center=(ix*dx, iy*dy, iz*dz)
                                 orientation=np.random.randn(3)*0,
                                 center=center((ix, dx), (iy, dy), (iz, dz))
                                )
              for ix in range(nx)
              for iy in range(ny)
              for iz in range(nz)]
    shape = Composite(shapes)
    fn = lambda qx, qy, view=view: cylinder_Iqxy(qx, qy, radius, length, view=view)
    check_shape_2d(shape, fn, view=view)

def check_sphere(radius=125, rho=2):
    shape = TriaxialEllipsoid(radius, radius, radius, rho)
    fn = lambda q: sphere_Iq(q, radius)
    check_shape(shape, fn)

def check_csbox(a=10, b=20, c=30, da=1, db=2, dc=3, slda=1, sldb=2, sldc=3, sld_core=4):
    core = Box(a, b, c, sld_core)
    side_a = Box(da, b, c, slda, center=((a+da)/2, 0, 0))
    side_b = Box(a, db, c, sldb, center=(0, (b+db)/2, 0))
    side_c = Box(a, b, dc, sldc, center=(0, 0, (c+dc)/2))
    side_a2 = copy(side_a).shift(-a-da, 0, 0)
    side_b2 = copy(side_b).shift(0, -b-db, 0)
    side_c2 = copy(side_c).shift(0, 0, -c-dc)
    shape = Composite((core, side_a, side_b, side_c, side_a2, side_b2, side_c2))
    def fn(q):
        return csbox_Iq(q, a, b, c, da, db, dc, slda, sldb, sldc, sld_core)
    #check_shape(shape, fn)

    view = (20, 30, 40)
    def fn_xy(qx, qy):
        return csbox_Iqxy(qx, qy, a, b, c, da, db, dc,
                          slda, sldb, sldc, sld_core, view=view)
    check_shape_2d(shape, fn_xy, view=view)

if __name__ == "__main__":
    check_cylinder(radius=10, length=40)
    #check_cylinder_2d(radius=10, length=40, view=(90,30,0))
    #check_cylinder_2d_lattice(radius=10, length=50, view=(90,30,0))
    #check_sphere()
    #check_csbox()
    #check_csbox(da=100, db=200, dc=300)
