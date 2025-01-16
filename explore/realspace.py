from __future__ import division, print_function

import cmath
import time
from copy import copy
import os
import argparse
from collections import OrderedDict
from timeit import default_timer as timer
from typing import Tuple
from inspect import getfullargspec

import numpy as np
from numpy import pi, radians, sin, cos, sqrt, clip
from numpy.random import poisson, uniform, randn, rand
from numpy.polynomial.legendre import leggauss
from scipy.integrate import simps
from scipy.special import j1 as J1
from scipy.special import gamma

try:
    from numba import njit, prange
    # SAS_NUMBA: 0=None, 1=CPU, 2=GPU
    SAS_NUMBA = int(os.environ.get("SAS_NUMBA", "1"))
    USE_NUMBA = SAS_NUMBA > 0
    USE_CUDA = SAS_NUMBA > 1
except ImportError:
    # Identity decorator @njit or @njit(...)
    njit = lambda f, *args, **kw: f if callable(f) else (lambda k: k)
    USE_NUMBA = USE_CUDA = False

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

def pol2rec(r, theta, phi):
    """
    Convert from 3D polar coordinates to rectangular coordinates.
    """
    theta, phi = radians(theta), radians(phi)
    x = +r * sin(theta) * cos(phi)
    y = +r * sin(theta) * sin(phi)
    z = +r * cos(theta)
    return x, y, z

def jitter(theta, phi, psi):
    r"""
    Return the jitter transform to rotate a set of points.
    View is in degrees using nautical angles with roll $\psi$ around $c$,
    pitch $\theta$ around $b$ and yaw $\phi$ around $a$.

    **unused**
    """
    return Rx(phi) @ Ry(theta) @ Rz(psi)

def rotation(theta, phi, psi):
    r"""
    Return a rotation matrix to apply to a set of points.
    View is in degrees using a $z$-$y$-$z$ rotation sequence of Euler angles
    $\phi$-$\theta$-$\psi$.  The $c$-axis of the shape starts along $z$ and
    the $b$-axis starts along $y$.
    """
    return Rz(phi) @ Ry(theta) @ Rz(psi)

def invert_view(qx, qy, view):
    r"""
    Return $(q_a, q_b, q_c)$ for the $(\theta, \phi, \psi)$ view angle at
    detector pixel corresponding to $(q_x, q_y)$.  View is in degrees using
    a $z$-$y$-$z$ sequence of Euler angles $\phi$-$\theta$-$\psi$.
    """
    theta, phi, psi = view
    Rinv = Rz(-psi) @ Ry(-theta) @ Rz(-phi)
    q = np.vstack((qx.flatten(), qy.flatten(), 0*qx.flatten()))
    return Rinv @ q

def apply_view(points, view):
    r"""
    Return $(p_x, p_y, p_x)$ rotated by the $(\theta, \phi, \psi)$ view angle.
    View is in degrees using a $z$-$y$-$z$ sequence of Euler angles
    $\phi$-$\theta$-$\psi$.
    """
    R = rotation(*view)
    return points @ R.T

class Shape:
    rotation = np.eye(3)
    center = np.array([0., 0., 0.])[:, None]
    r_max = None
    is_magnetic = False

    def volume(self):
        # type: () -> float
        raise NotImplementedError()

    def sample(self, density):
        # type: (float) -> Tuple[np.ndarray, np.ndarray]
        """
        Returns arrays (rho[N], points[N, 3]).
        """
        raise NotImplementedError()

    def dims(self):
        # type: () -> Tuple[float, float, float]
        raise NotImplementedError()

    def rotate(self, theta, phi, psi):
        """See :func:`rotation` for details on the rotation matrix."""
        self.rotation = rotation(theta, phi, psi) @ self.rotation
        return self

    def shift(self, x, y, z):
        self.center = self.center + np.array([x, y, z])[:, None]
        return self

    def _adjust(self, points):
        points = self.rotation @ points.T + self.center
        return points.T

    def r_bins(self, q, over_sampling=1, r_step=None):
        return r_bins(q, r_max=self.r_max, r_step=r_step,
                      over_sampling=over_sampling)

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

class Superball(Shape):
    def __init__(self, a, p,
                 value, center=(0, 0, 0), orientation=(0, 0, 0)):
        self.value = np.asarray(value)
        self.rotate(*orientation)
        self.shift(*center)
        self.a, self.p = a, p
        self._scale = a/2
        # Solve for rounded corner radius x = y = z:
        #    x^2p + y^2p + z^2p = 3 x^2p = (a/2)^2p
        #    => x = a / 2 root[2p](3)
        #    => d = 2r = 2x root(3) = root(3)/root[2p](3) a = 3^(p-1)/2p a
        #self.r_max = 3**((p-1)/(2*p)) * a  # Too short---don't know why.
        self.r_max = sqrt(3)*a
        self.dims = a, a, a
        g1 = gamma(1.0 / (2.0 * p))
        g3 = gamma(3.0 / (2.0 * p))
        self.volume = a**3 / 12.0 / p**2 * g1**3 / g3

    def sample(self, density):
        # Sample from cube[-a/2, a/2]
        num_points = poisson(density*self.a**3)
        points = uniform(-1, 1, size=(num_points, 3))
        # Trim points outside maximum "squared radius", x^2p + y^2p + z^2p < 1 
        radius_sq = np.sum((points**2)**self.p, axis=1)
        points = points[radius_sq <= 1]
        values = self.value.repeat(points.shape[0])
        return values, self._adjust(self._scale*points)

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
        radius_sq = points[:, 0]**2 + points[:, 1]**2
        points = points[radius_sq <= 1]
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
        values = np.full_like(points[:, 0], self.value)
        # then set value to face value if |z| > face/(length/2))
        values[abs(points[:, 2]) > self.length/(self.length + 2*self.thick_face)] = self.value_face
        # finally set value to rim value if outside the core ellipse
        radius = (points[:, 0]**2*(1 + self.thick_rim/self.ra)**2
                  + points[:, 1]**2*(1 + self.thick_rim/self.rb)**2)
        values[radius>1] = self.value_rim
        return values, self._adjust(self._scale*points)

class TruncatedSphere(Shape):
    """
    Sphere of radius r, with points z < -h truncated.
    """
    def __init__(self, r, h, value, center=(0, 0, 0), orientation=(0, 0, 0)):
        self.value = np.asarray(value)
        self.rotate(*orientation)
        self.shift(*center)
        self.r, self.h = r, h
        # Max distance between points in the shape is the maximum diameter
        self.r_max = 2*r if h >= 0 else 2*sqrt(r**2 - h**2)
        self.dims = self.r_max, self.r_max, r+h
        self.volume = pi*(2*r**3/3 + r**2*h - h**3/3)
        #Vp = pi*(2*r**3/3 + r**2*h - h**3/3)
        #Vm = pi*(2*r**3/3 - r**2*h + h**3/3)
        #Vd = Vp + Vm - 4*pi*r**3/3

    def sample(self, density):
        num_points = poisson(density*np.prod(self.dims))
        points = uniform(-1, 1, size=(num_points, 3))
        # Translate U ~ [-1, 1] in x,y to [-r_trunc/r, r_trunc/r] when
        # truncation starts above the equator, otherwise leave it at [-1, 1].
        # This makes for more efficient sampling since we don't have to
        # consider the maximum diameter.  We already calculated r_max as
        # 2*r_trunc in this case, so just use the ratio of r_max to 2*r.
        points[:, 0:2] *= 0.5*self.r_max/self.r
        # Translate U ~ [-1, 1] in z to [-h/r, 1], with h representing
        # distance below equator.  So:
        #    (U + 1)/2 => [0, 1]
        #    [0, 1] * (1+h/r) => [0, 1+h/r]
        #    [0, 1+h/r] - h/r => [-h/r, 1]
        # Combining:
        #    (U + 1)/2 * (1+h/r) - h/r
        #       = U*(1+h/r)/2 + (1+h/r)/2 - h/r
        #       = U*(1/2 + h/2r) + 1/2 + h/2r - 2h/2r
        ratio = 0.5*self.h/self.r
        points[:, 2] *= (0.5 + ratio)
        points[:, 2] += (0.5 - ratio)
        radius = np.sum(points**2, axis=1)
        points = self.r*points[radius<=1]
        values = self.value.repeat(points.shape[0])
        return values, self._adjust(points)

class TriaxialEllipsoid(Shape):
    def __init__(self, ra, rb, rc,
                 value, center=(0, 0, 0), orientation=(0, 0, 0),
                 magnetism=None):
        self.is_magnetic = (magnetism is not None)
        self.value = np.asarray(value)
        self.magnetism = magnetism if self.is_magnetic else (0., 0., 0.)
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
        radius_sq = np.sum(points**2, axis=1)
        points = self._scale*points[radius_sq <= 1]
        values = self.value.repeat(points.shape[0])
        return values, self._adjust(points)

    def sample_magnetic(self, density):
        values, points = self.sample(density)
        magnetism = np.tile(self.magnetism, (points.shape[0], 1)).T
        return values, magnetism, points

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
        radius_sq = points[:, 0]**2 + points[:, 1]**2
        points = points[radius_sq <= 1]

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

def barbell(r=20, rbell=50, length=20, rho=2):
    h = sqrt(rbell**2 - r**2)
    top = TruncatedSphere(rbell, h, value=rho, center=(0, 0, length/2+h))
    rod = EllipticalCylinder(r, r, length, value=rho)
    bottom = TruncatedSphere(rbell, h, value=rho, center=(0, 0, -length/2-h),
                             orientation=(180, 0, 0))
    shape = Composite((top, rod, bottom))
    shape.dims = 2*rbell, 2*rbell, length+2*(rbell+h)
    # r_max should be total length?
    shape.r_max = (length + 2*(rbell + h))
    return shape

def capped_cylinder(r=20, rcap=50, length=20, rho=2):
    h = -sqrt(rcap**2 - r**2)
    top = TruncatedSphere(rcap, h, value=rho, center=(0, 0, length/2+h))
    rod = EllipticalCylinder(r, r, length, value=rho)
    bottom = TruncatedSphere(rcap, h, value=rho, center=(0, 0, -length/2-h),
                             orientation=(180, 0, 0))
    shape = Composite((top, rod, bottom))
    shape.dims = 2*r, 2*r, length+2*(rcap+h)
    # r_max is the length of the diagonal + height of the cap for safety.
    # This is a bit larger than necessary, but that's better than truncation.
    shape.r_max = sqrt(length**2 + 4*r**2) + (rcap + h)
    return shape

def _Iqabc(weight, x, y, z, qa, qb, qc):
    """I(q) = |sum V(r) rho(r) e^(1j q.r)|^2 / sum V(r)"""
    #print("calling python")
    Iq = [abs(np.sum(weight*np.exp(1j*(qa_k*x + qb_k*y + qc_k*z))))**2
          for qa_k, qb_k, qc_k in zip(qa.flat, qb.flat, qc.flat)]
    return np.asarray(Iq)
_Iqabcf = _Iqabc

if USE_NUMBA:
    # Override simple numpy solution with numba if available
    def _Iqabc_py(weight, x, y, z, qa, qb, qc):
        #print("calling numba")
        Iq = np.empty_like(qa)
        for j in prange(len(Iq)):
            #total = 0. + 0j
            #for k in range(len(weight)):
            #    total += weight[k]*np.exp(1j*(qa[j]*x[k] + qb[j]*y[k] + qc[j]*z[k]))
            total = np.sum(weight * np.exp(1j*(qa[j]*x + qb[j]*y + qc[j]*z)))
            Iq[j] = abs(total)**2
        return Iq
    sig = "f8[:](f8[:],f8[:],f8[:],f8[:],f8[:],f8[:],f8[:])"
    _Iqabc = njit(sig, parallel=True, fastmath=True)(_Iqabc_py)
    _Iqabcf = njit(sig.replace("f8", "f4"), parallel=True, fastmath=True)(_Iqabc_py)

if USE_CUDA:
    # delayed loading of cuda
    _IQABC_CUDA_KERNELS = {}
    def _get_Iqabc_kernel(dtype):
        #print("calling cuda")
        if not _IQABC_CUDA_KERNELS:
            from numba import cuda
            if not cuda.list_devices():
                raise RuntimeError("no cuda devices found")
            def _kernel_py(weight, x, y, z, qa, qb, qc, Iq):
                j = cuda.grid(1)
                if j < qa.size:
                    total = 0. + 0j
                    for k in range(x.size):
                        total += weight[k]*cmath.exp(1j*(qa[j]*x[k] + qb[j]*y[k] + qc[j]*z[k]))
                    Iq[j] = abs(total)**2
            sig_d = "void(f8[:],f8[:],f8[:],f8[:],f8[:],f8[:],f8[:],f8[:])"
            sig_f = sig_d.replace("f8", "f4")
            kernel_f = cuda.jit(sig_f, parallel=True, fastmath=True)(_kernel_py)
            kernel_d = cuda.jit(sig_d, parallel=True, fastmath=True)(_kernel_py)
            _IQABC_CUDA_KERNELS['f'] = kernel_f
            _IQABC_CUDA_KERNELS['d'] = kernel_d
        kernel = _IQABC_CUDA_KERNELS[dtype.char]
        return kernel

    def _Iqabc(weight, x, y, z, qa, qb, qc):
        Iq = np.empty_like(qa)
        # Apparently numba deals with all the necessary padding of vectors to nice boundaries
        # before transfering to the GPU, so we don't need to do so by hand here.
        threadsperblock = 32
        blockspergrid = (Iq.size + (threadsperblock - 1)) // threadsperblock
        kernel = _get_Iqabc_kernel(qa.dtype)
        kernel[blockspergrid, threadsperblock](weight, x, y, z, qa, qb, qc, Iq)
        return Iq
    _Iqabcf = _Iqabc


if 0 and USE_CUDA:
    ### *** DEPRECATED ***
    ### Variant on the kernel with padding of vectors that is no faster and doesn't appear
    ### to be more correct.  Leave it around for now in case we decide we don't trust numba.
    _IQABC_CUDA_KERNELS = {}
    def _get_Iqabc_kernel(dtype):
        #print("calling cuda")
        if not _IQABC_CUDA_KERNELS:
            from numba import cuda
            if not cuda.list_devices():
                raise RuntimeError("no cuda devices found")
            def _kernel_py(nx, nq, weight, x, y, z, qa, qb, qc, Iq):
                j = cuda.grid(1)
                if j < nq:
                    total = 0. + 0j
                    for k in range(nx):
                        total += weight[k]*cmath.exp(1j*(qa[j]*x[k] + qb[j]*y[k] + qc[j]*z[k]))
                    Iq[j] = abs(total)**2
            sig_d = "void(i4,i4,f8[:],f8[:],f8[:],f8[:],f8[:],f8[:],f8[:],f8[:])"
            sig_f = sig_d.replace("f8", "f4")
            kernel_f = cuda.jit(sig_f, parallel=True, fastmath=True)(_kernel_py)
            kernel_d = cuda.jit(sig_d, parallel=True, fastmath=True)(_kernel_py)
            _IQABC_CUDA_KERNELS['f'] = kernel_f
            _IQABC_CUDA_KERNELS['d'] = kernel_d
        kernel = _IQABC_CUDA_KERNELS[dtype.char]
        return kernel

    def _Iqabc(weight, x, y, z, qa, qb, qc):
        kernel = _get_Iqabc_kernel(qa.dtype)
        nx, nq = len(x), len(qa)
        threadsperblock = 32
        blockspergrid = (nq + (threadsperblock - 1)) // threadsperblock
        weight, x, y, z = pad_vectors(4, weight, x, y, z)
        qa, qb, qc = pad_vectors(threadsperblock, qa, qb, qc)
        Iq = np.empty_like(qa)
        kernel[blockspergrid, threadsperblock](nx, nq, weight, x, y, z, qa, qb, qc, Iq)
        return Iq[:nq]
    _Iqabcf = _Iqabc

    def pad_vectors(boundary, *vectors):
        """
        Yields a list of vectors padded with NaN to a multiple of *boundary*.
        Yields the original vector if the size is already a mulitple of *boundary*.
        """
        for old in vectors:
            old_size = len(old)
            new_size = ((old_size + boundary-1)//boundary)*boundary
            if new_size > old_size:
                new = np.empty(new_size, dtype=old.dtype)
                new[:old_size] = old
                new[old_size:] = np.nan
                yield new
            else:
                yield old

def calc_Iqxy(qx, qy, rho, points, volume=1.0, view=(0, 0, 0), dtype='f'):
    """
    *qx*, *qy* correspond to the detector pixels at which to calculate the
    scattering, relative to the beam along the negative z axis.
    *points* are three columns (x, y, z), one for each sample in the shape.
    *rho* (1e-6/Ang) is the scattering length density of each point.
    *volume* should be 1/number_density.  That is, each of n particles in the
    total value represents volume/n contribution to the scattering.
    *view* rotates the points about the axes using Euler angles for pitch
    yaw and roll for a beam travelling along the negative z axis.
    *dtype* is the numerical precision of the calculation.
    """
    # TODO: maybe slightly faster to rotate points, and drop qc*z
    qx, qy = np.broadcast_arrays(qx, qy)
    qa, qb, qc = invert_view(qx, qy, view)
    rho, volume = np.broadcast_arrays(rho, volume)
    weight = rho*volume
    x, y, z = points.T

    # I(q) = |sum V(r) rho(r) e^(1j q.r)|^2 / sum V(r)
    if np.dtype(dtype) == np.float64:
        weight, x, y, z, qa, qb, qc = [np.asarray(v, 'd') for v in (weight, x, y, z, qa, qb, qc)]
        Iq = _Iqabc(weight, x, y, z, qa.flatten(), qb.flatten(), qc.flatten())
    else:  # float32
        weight, x, y, z, qa, qb, qc = [np.asarray(v, 'f') for v in (weight, x, y, z, qa, qb, qc)]
        Iq = _Iqabcf(weight, x, y, z, qa.flatten(), qb.flatten(), qc.flatten())
    # The scale factor 1e-4 is due to the conversion from rho = 1e-6 squared
    # times the conversion of 1e-8 from inverse angstroms to inverse cm.
    return np.asarray(Iq).reshape(qx.shape) * (1e-4 / np.sum(volume))

def spin_weights(in_spin, out_spin):
    """
    Compute spin cross sections given in_spin and out_spin
    To convert spin cross sections to sld b:
        uu * (sld - m_sigma_x);
        dd * (sld + m_sigma_x);
        ud * (m_sigma_y - 1j*m_sigma_z);
        du * (m_sigma_y + 1j*m_sigma_z);
    weights for spin crosssections: dd du real, ud real, uu, du imag, ud imag
    """

    in_spin = clip(in_spin, 0.0, 1.0)
    out_spin = clip(out_spin, 0.0, 1.0)
    # Previous version of this function took the square root of the weights,
    # under the assumption that
    #
    #     w*I(q, rho1, rho2, ...) = I(q, sqrt(w)*rho1, sqrt(w)*rho2, ...)
    #
    # However, since the weights are applied to the final intensity and
    # are not interned inside the I(q) function, we want the full
    # weight and not the square root.  Anyway no function will ever use
    # set_spin_weights as part of calculating an amplitude, as the weights are
    # related to polarisation efficiency of the instrument. The weights serve to
    # construct various magnet scattering cross sections, which are linear combinations
    # of the spin-resolved cross sections. The polarisation efficiency e_in and e_out
    # are parameters ranging from 0.5 (unpolarised) beam to 1 (perfect optics).
    # For in_spin or out_spin <0.5 one assumes a CS, where the spin is reversed/flipped
    # with respect to the initial supermirror polariser. The actual polarisation efficiency
    # in this case is however e_in/out = 1-in/out_spin.

    norm = 1 - out_spin if out_spin < 0.5 else out_spin

    # The norm is needed to make sure that the scattering cross sections are
    # correctly weighted, such that the sum of spin-resolved measurements adds up to
    # the unpolarised or half-polarised scattering cross section. No intensity weighting
    # needed on the incoming polariser side (assuming that a user), has normalised
    # to the incoming flux with polariser in for SANSPOl and unpolarised beam, respectively.

    weight = (
        (1.0 - in_spin) * (1.0 - out_spin) / norm, # dd
        (1.0 - in_spin) * out_spin / norm,       # du
        in_spin * (1.0 - out_spin) / norm,       # ud
        in_spin * out_spin / norm,             # uu
    )
    return weight

def orth(A, b_hat): # A = 3 x n, and b_hat unit vector
    #return A - np.sum(A*b_hat[:, None], axis=0)[None, :]*b_hat[:, None]
    return A - np.outer(b_hat, b_hat)@A


def magnetic_sld(qx, qy, up_theta, up_phi, rho, rho_m):
    """
    Compute the complex sld for the magnetic spin states.
    Returns effective rho for spin states [dd, du, ud, uu].
    """
    # For q=0 one would see the demagnetising field of the sample, equivalent
    # to direction q_hat = [sqrt(1/2), sqrt(1/2), 0] for a disc shaped sample
    # that is very thin along the beam.
    # Note: This is different from kernel_iq.c, which sets I(0, 0) to zero.
    q_norm = sqrt(qx**2 + qy**2)
    if abs(q_norm) < 1.e-16:
        q_hat = np.array([1., 1., 0.]) / np.sqrt(2)
    else:
        q_hat = np.array([qx, qy, 0]) / q_norm
    M_perp = orth(rho_m, q_hat)  # M = rho_m

    # perpy_hat and perpz_hat are unit vectors spanning up the plane
    # perpendicular to polarisation for SF scattering (avoiding
    # repetitive computation of orthogonal vectors)
    cos_theta, sin_theta = cos(radians(up_theta)), sin(radians(up_theta))
    cos_phi, sin_phi = cos(radians(up_phi)), sin(radians(up_phi))
    p_hat = np.array([sin_theta * cos_phi, sin_theta * sin_phi, cos_theta])
    perpy_hat = np.array([-sin_phi, cos_phi, 0])
    perpz_hat = np.array([-cos_theta * cos_phi, -cos_theta * sin_phi, sin_theta])

    perpx = p_hat @ M_perp
    perpy = perpy_hat @ M_perp
    perpz = perpz_hat @ M_perp

    return (
        rho - perpx,   # dd => sld - D M_perpx
        perpy - 1j * perpz, # du => -D (M_perpy + j M_perpz)
        perpy + 1j * perpz, # ud => -D (M_perpy - j M_perpz)
        rho + perpx,   # uu => sld + D M_perpx
    )

# TODO: provide numba and cuda version
def calc_Iq_magnetic(qx, qy, rho, rho_m, points, volume=1.0, view=(0, 0, 0),
                     up_frac_i=0.5, up_frac_f=0.5, up_theta=0., up_phi=0.):
    """
    *qx*, *qy* correspond to the detector pixels at which to calculate the
    scattering, relative to the beam along the negative z axis.
    *points* are three columns (x, y, z), one for each sample in the shape.
    *rho* (1e-6/Ang) is the scattering length density of each point.
    *rho_m* (1e-6/Ang) are the (mx, my, mz) components of the magnetic
    scattering length density for each point.
    *volume* should be 1/number_density.  That is, each of n particles in the
    total value represents volume/n contribution to the scattering.
    *view* rotates the points about the axes using Euler angles for pitch
    yaw and roll for a beam travelling along the negative z axis.
    *up_frac_i* is the portion of polarizer neutrons which are spin up.
    *up_frac_f* is the portion of analyzer neutrons which are spin up.
    *up_theta* is the inclination from the beam direction (z-axis).
    *up_phi* is the rotation in the detector plane.
    *dtype* is the numerical precision of the calculation. [not implemented]
    """
    # TODO: maybe slightly faster to rotate points and rho_m, and drop qc*z
    qx, qy = np.broadcast_arrays(qx, qy)
    qa, qb, qc = invert_view(qx, qy, view)
    rho, volume = np.broadcast_arrays(rho, volume)
    weights = spin_weights(up_frac_i, up_frac_f)

    # I(q) = |sum V(r) rho(r) e^(1j q.r)|^2 / sum V(r)
    shape = qx.shape
    Iq = np.zeros(qx.size, 'd')
    x, y, z = points.T
    qx, qy = (v.flatten() for v in (qx, qy))
    for k in range(qx.size):
        ephase = volume*np.exp(1j*(qa[k]*x + qb[k]*y + qc[k]*z))
        dd, du, ud, uu = magnetic_sld(qx[k], qy[k], up_theta, up_phi, rho, rho_m)
        for w, xs in zip(weights, (dd, du, ud, uu)):
            if w == 0.0:
                continue
            Iq[k] += w * abs(np.sum(xs*ephase))**2
    # The scale factor 1e-4 is due to the conversion from rho = 1e-6 squared
    # times the conversion of 1e-8 from inverse angstroms to inverse cm.
    return np.asarray(Iq).reshape(shape) * (1e-4 / np.sum(volume))

def _calc_Pr_nonuniform(r, rho, points, volume):
    # Make Pr a little be bigger than necessary so that only distances
    # min < d < max end up in Pr
    n_max = len(r)+1
    extended_Pr = np.zeros(n_max+1, 'd')
    # r refers to bin centers; find corresponding bin edges
    bins = bin_edges(r)
    t_next = timer() + 3
    for k, rho_k in enumerate(rho[:-1]):
        distance = np.linalg.norm(points[k] - points[k+1:], axis=1)
        weights = (rho_k * volume[k]) * (rho[k+1:] * volume[k+1:])
        #weights = (rho_k * volume[k]) * rho[k+1:]
        index = np.searchsorted(bins, distance)
        # Note: indices may be duplicated, so "Pr[index] += w" will not work!!
        extended_Pr += np.bincount(index, weights, n_max+1)
        t = timer()
        if t > t_next:
            t_next = t + 3
            print("processing %d of %d"%(k, len(rho)-1))
    Pr = extended_Pr[1:-1]
    return Pr

def _calc_Pr_uniform(r, rho, points, volume):
    # Make Pr a little be bigger than necessary so that only distances
    # min < d < max end up in Pr
    dr, n_max = r[0], len(r)
    extended_Pr = np.zeros(n_max+1, 'd')
    t0 = timer()
    t_next = t0 + 3
    for k, rho_k in enumerate(rho[:-1]):
        distance = np.linalg.norm(points[k] - points[k+1:], axis=1)
        #weights = (rho_k * volume[k]) * (rho[k+1:] * volume[k+1:])
        weights = rho_k * rho[k+1:] * (volume[k] + volume[k+1:])
        index = np.minimum(np.asarray(distance/dr, 'i'), n_max)
        # Note: indices may be duplicated, so "Pr[index] += w" will not work!!
        extended_Pr += np.bincount(index, weights, n_max+1)
        t = timer()
        if t > t_next:
            t_next = t + 3
            print("processing %d of %d"%(k, len(rho)-1))
    #print("time py:", timer() - t0)
    Pr = extended_Pr[:-1]
    #print("vol", np.sum(volume))
    return Pr*1e-4

    # Can get an additional 2x by going to C. Cuda/OpenCL will allow even
    # more speedup, though still bounded by the O(n^2) cost.
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
    #@njit("f8[:](f8[::1], f8[::1], f8[::1,:], f8[:])", parallel=True, fastmath=True)
    @njit(parallel=True, fastmath=True)
    def _calc_Pr_uniform(r, rho, points, volume):
        dr = r[0]
        n_max = len(r)
        Pr = np.zeros_like(r)
        for j in prange(len(rho) - 1):
            x, y, z = points[j, 0], points[j, 1], points[j, 2]
            rho_j, volume_j = rho[j], volume[j]
            for k in range(j+1, len(rho)):
                distance = sqrt((x - points[k, 0])**2
                                + (y - points[k, 1])**2
                                + (z - points[k, 2])**2)
                index = int(distance/dr)
                if index < n_max:
                    Pr[index] += rho_j*rho[k]*(volume_j + volume[k])
        return Pr


def calc_Pr(r, rho, points, volume):
    # P(r) with uniform steps in r is 3x faster; check if we are uniform
    # before continuing
    r, points = [np.asarray(v, 'd') for v in (r, points)]
    npoints = points.shape[0]
    rho = np.broadcast_to(np.asarray(rho, 'd'), npoints)
    volume = np.broadcast_to(np.asarray(volume, 'd'), npoints)
    if np.max(np.abs(np.diff(r) - r[0])) > r[0]*0.01:
        Pr = _calc_Pr_nonuniform(r, rho, points, volume)
    else:
        Pr = _calc_Pr_uniform(r, rho, points, volume)
    # Note: 1e-4 because (1e-6 rho)^2 = 1e-12 rho^2 time 1e-8 for 1/A to 1/cm
    return Pr * 1e-4


def r_bins(q, r_max=None, r_step=None, over_sampling=1):
    if r_max is None:
        r_max = 2 * pi / q[0]
    if r_step is None:
        r_step = 2 * pi / q[-1] / over_sampling
    return np.arange(r_step, r_max, r_step)


def j0(x):
    # use q/pi since np.sinc = sin(pi x)/(pi x)
    return np.sinc(x/np.pi)


def calc_Iq_from_Pr(q, r, Pr):
    Iq = np.array([simps(Pr * j0(qk*r), r) for qk in q])
    #Iq = np.array([np.trapz(Pr * j0(qk*r), r) for qk in q])
    #Iq /= Iq[0]
    return Iq


def _calc_Iq_avg(Iq, q, r, sld, volume):
    weight = sld * volume
    for i, qi in enumerate(q):
        Fq = np.sum(weight * np.sinc((qi/np.pi)*r))
        Iq[i] = Fq**2
if USE_NUMBA:
    #sig = njit('(f8[:], f8[:], f8[:], f8[:], f8[:])', parallel=True, fastmath=True)
    sig = njit(parallel=True, fastmath=True)
    _calc_Iq_avg = sig(_calc_Iq_avg)


def calc_Iq_avg(q, rho, points, volume=1.0):
    # Centralize the data
    center = 0.5*(np.min(points, axis=0, keepdims=True)
                  + np.max(points, axis=0, keepdims=True))
    points = points - center
    # Find distance from center
    r = np.linalg.norm(points, axis=1)
    # Call calculator
    Iq = np.empty_like(q)
    rho = np.broadcast_to(np.asarray(rho, 'd'), points.shape[:1])
    volume = np.broadcast_to(np.asarray(volume, 'd'), points.shape[:1])
    _calc_Iq_avg(Iq, q, r, rho, volume)
    return Iq * (1e-4/np.sum(volume))

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
def plot_calc(r, Pr, q, Iq, theory=None, title=None, Iq_avg=None):
    import matplotlib.pyplot as plt
    plt.subplot(211)
    plt.plot(r, Pr, '-', label="Pr")
    plt.xlabel('r (A)')
    plt.ylabel('Pr (1/A^2)')
    if title is not None:
        plt.title(title)
    plt.grid(True)
    plt.subplot(212)
    plt.loglog(q, Iq, '-', label='from Pr')
    #plt.loglog(q, Iq/theory[1], '-', label='Pr/theory')
    if Iq_avg is not None:
        plt.loglog(q, Iq_avg, '-', label='from Iq_avg')
    plt.xlabel('q (1/A)')
    plt.ylabel('Iq')
    plt.grid(True)
    if theory is not None:
        #plt.loglog(theory[0], theory[1]/theory[1][0], '-', label='analytic')
        plt.loglog(theory[0], theory[1], '-', label='analytic')
        plt.legend()

def plot_calc_2d(qx, qy, Iqxy, theory=None, title=None):
    import matplotlib.pyplot as plt
    qx, qy = bin_edges(qx), bin_edges(qy)
    #qx, qy = np.meshgrid(qx, qy)
    if theory is not None:
        plt.subplot(131)
    #plt.pcolor(qx, qy, np.log10(Iqxy))
    extent = [qx[0], qx[-1], qy[0], qy[-1]]
    plt.imshow(np.log10(Iqxy), extent=extent, interpolation="nearest",
               origin='lower')
    plt.colorbar()
    plt.xlabel('qx (1/A)')
    plt.ylabel('qy (1/A)')
    plt.axis('equal')
    plt.axis(extent)
    #plt.grid(True)
    if title is not None:
        plt.title(title)
    if theory is not None:
        plt.subplot(132)
        # Skip bad values in theory
        index = np.isnan(theory)
        theory[index] = Iqxy[index]
        plt.imshow(np.log10(theory), extent=extent, interpolation="nearest",
                   origin='lower')
        plt.title("theory")
        plt.colorbar()
        plt.axis('equal')
        plt.axis(extent)
        plt.xlabel('qx (1/A)')

    if theory is not None:
        plt.subplot(133)
        rel = (theory-Iqxy)/theory
        plt.imshow(rel, extent=extent, interpolation="nearest", origin='lower')
        plt.colorbar()
        plt.axis('equal')
        plt.axis(extent)
        plt.xlabel('qx (1/A)')
        plt.title('max rel. err=%g' % np.max(abs(rel)))

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

def _sasmodels_Iq(kernel, q, pars):
    from sasmodels.data import empty_data1D
    from sasmodels.direct_model import DirectModel
    data = empty_data1D(q)
    calculator = DirectModel(data, kernel)
    Iq = calculator(**pars)
    return Iq

def _sasmodels_Iqxy(kernel, qx, qy, pars, view):
    from sasmodels.data import Data2D
    from sasmodels.direct_model import DirectModel
    Iq = np.full_like(qx, 100)
    data = Data2D(x=qx, y=qy, z=Iq, dx=None, dy=None, dz=np.sqrt(Iq))
    data.x_bins = qx[0, :]
    data.y_bins = qy[:, 0]
    data.filename = "fake data"

    calculator = DirectModel(data, kernel)
    pars_plus_view = pars.copy()
    pars_plus_view.update(theta=view[0], phi=view[1], psi=view[2])
    Iqxy = calculator(**pars_plus_view)
    # calculator avoids masked values; instead set masked values to NaN
    result = np.empty_like(qx)
    result[calculator.index] = Iqxy
    result[~calculator.index] = np.nan
    return result

def wrap_sasmodel(name, **pars):
    from sasmodels.core import load_model
    kernel = load_model(name)
    fn = lambda q: _sasmodels_Iq(kernel, q, pars)
    fn_xy = lambda qx, qy, view: _sasmodels_Iqxy(kernel, qx, qy, pars, view)
    return fn, fn_xy


# --------- Test cases -----------

def build_box(a=10, b=20, c=30, rho=2.):
    shape = Box(a, b, c, rho)
    fn = lambda q: box_Iq(q, a, b, c)*rho**2
    fn_xy = lambda qx, qy, view: box_Iqxy(qx, qy, a, b, c, view=view)*rho**2
    return shape, fn, fn_xy

def build_superball(a=10, p=3, rho=2.):
    shape = Superball(a, p, rho)
    fn, fn_xy = wrap_sasmodel(
        'superball',
        scale=1,
        background=0,
        length_a=a,
        exponent_p=p,
        sld=rho,
        sld_solvent=0,
    )
    return shape, fn, fn_xy

def build_csbox(a=10, b=20, c=30, da=1, db=2, dc=3, slda=1, sldb=2, sldc=3, sld_core=4):
    shape = csbox(a, b, c, da, db, dc, slda, sldb, sldc, sld_core)
    fn = lambda q: csbox_Iq(q, a, b, c, da, db, dc, slda, sldb, sldc, sld_core)
    fn_xy = lambda qx, qy, view: csbox_Iqxy(qx, qy, a, b, c, da, db, dc,
                                            slda, sldb, sldc, sld_core, view=view)
    return shape, fn, fn_xy

def build_sphere(radius=125, rho=2,
                 rho_m=0, theta_m=0, phi_m=0, up_i=0, up_f=0, up_theta=0, up_phi=0):
    magnetism = pol2rec(rho_m, theta_m, phi_m) if rho_m != 0.0 else None
    shape = TriaxialEllipsoid(radius, radius, radius, rho, magnetism=magnetism)
    shape.spin = (up_i, up_f, up_theta, up_phi)
    fn, fn_xy = wrap_sasmodel(
        'sphere',
        scale=1,
        background=0,
        radius=radius,
        sld=rho,
        sld_solvent=0,
        sld_M0=rho_m,
        sld_mtheta=theta_m,
        sld_mphi=phi_m,
        up_frac_i=up_i,
        up_frac_f=up_f,
        up_theta=up_theta,
        up_phi=up_phi,
    )
    return shape, fn, fn_xy

def build_ellip(rab=125, rc=50, rho=2,
                rho_m=0, theta_m=0, phi_m=0, up_i=0, up_f=0, up_theta=0, up_phi=0):
    magnetism = pol2rec(rho_m, theta_m, phi_m) if rho_m != 0.0 else None
    shape = TriaxialEllipsoid(rab, rab, rc, rho, magnetism=magnetism)
    # TODO: polarization spec doesn't belong in shape
    # Put spin state info into the shape since we have it available.
    shape.spin = (up_i, up_f, up_theta, up_phi)
    fn, fn_xy = wrap_sasmodel(
        'ellipsoid',
        scale=1,
        background=0,
        radius_equatorial=rab,
        radius_polar=rc,
        sld=rho,
        sld_solvent=0,
        sld_M0=rho_m,
        sld_mtheta=theta_m,
        sld_mphi=phi_m,
        up_frac_i=up_i,
        up_frac_f=up_f,
        up_theta=up_theta,
        up_phi=up_phi,
    )
    return shape, fn, fn_xy

def build_triell(ra=125, rb=200, rc=50, rho=2,
                 rho_m=0, theta_m=0, phi_m=0, up_i=0, up_f=0, up_theta=0, up_phi=0):
    magnetism = pol2rec(rho_m, theta_m, phi_m) if rho_m != 0.0 else None
    shape = TriaxialEllipsoid(ra, rb, rc, rho, magnetism=magnetism)
    shape.spin = (up_i, up_f, up_theta, up_phi)
    fn, fn_xy = wrap_sasmodel(
        'triaxial_ellipsoid',
        scale=1,
        background=0,
        radius_equat_minor=ra,
        radius_equat_major=rb,
        radius_polar=rc,
        sld=rho,
        sld_solvent=0,
        sld_M0=rho_m,
        sld_mtheta=theta_m,
        sld_mphi=phi_m,
        up_frac_i=up_i,
        up_frac_f=up_f,
        up_theta=up_theta,
        up_phi=up_phi,
    )
    return shape, fn, fn_xy

def build_cylinder(radius=25, length=125, rho=2.):
    shape = EllipticalCylinder(radius, radius, length, rho)
    fn = lambda q: cylinder_Iq(q, radius, length)*rho**2
    fn_xy = lambda qx, qy, view: cylinder_Iqxy(qx, qy, radius, length, view=view)*rho**2
    return shape, fn, fn_xy

def build_truncated_sphere(radius=25, h=0.5, rho=2.):
    shape = TruncatedSphere(radius, h=radius*h, value=rho)
    return shape, None, None

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

def build_barbell(r=20, rbell=50, length=20, rho=2.):
    shape = barbell(r, rbell, length, rho)
    fn, fn_xy = wrap_sasmodel(
        'barbell',
        scale=1,
        background=0,
        radius=r,
        radius_bell=rbell,
        length=length,
        sld=rho,
        sld_solvent=0,
    )
    return shape, fn, fn_xy

def build_capcyl(r=20, rcap=50, length=20, rho=2.):
    shape = capped_cylinder(r, rcap, length, rho)
    fn, fn_xy = wrap_sasmodel(
        'capped_cylinder',
        scale=1,
        background=0,
        radius=r,
        radius_cap=rcap,
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

def build_cubic_lattice(shape, nx=1, ny=1, nz=1, dx=2, dy=2, dz=2,
                        shuffle=0, rotate=0):
    a, b, c = shape.dims
    shapes = [copy(shape)
              .shift((ix+(randn() if shuffle < 0.3 else rand())*shuffle)*dx*a,
                     (iy+(randn() if shuffle < 0.3 else rand())*shuffle)*dy*b,
                     (iz+(randn() if shuffle < 0.3 else rand())*shuffle)*dz*c)
              .rotate(*((randn(3) if rotate < 30 else rand(3))*rotate))
              for ix in range(nx)
              for iy in range(ny)
              for iz in range(nz)]
    lattice = Composite(shapes)
    return lattice

SHAPE_FUNCTIONS = OrderedDict([
    ("cyl", build_cylinder),
    ("ellcyl", build_ellcyl),
    ("ellip", build_ellip),
    ("triell", build_triell),
    ("barbell", build_barbell),
    ("capcyl", build_capcyl),
    ("sphere", build_sphere),
    ("tsphere", build_truncated_sphere),
    ("box", build_box),
    ("csbox", build_csbox),
    ("cscyl", build_cscyl),
    ("superball", build_superball),
])
SHAPES = list(SHAPE_FUNCTIONS.keys())

def check_shape(title, shape, fn=None, show_points=False,
                mesh=100, qmax=1.0, r_step=0.01, samples=5000):
    rho_solvent = 0
    qmin = qmax/100.
    q = np.logspace(np.log10(qmin), np.log10(qmax), mesh)
    r = shape.r_bins(q, r_step=r_step)
    sampling_density = samples / shape.volume
    rho, points = shape.sample(sampling_density)
    volume = shape.volume / len(points)
    t0 = timer()
    Pr = calc_Pr(r, rho-rho_solvent, points, volume)
    print("calc Pr time", timer() - t0)
    Iq = calc_Iq_from_Pr(q, r, Pr)
    t0 = timer()
    Iq_avg = calc_Iq_avg(q, rho-rho_solvent, points, volume)
    print("calc Iq_avg time", timer() - t0)
    theory = (q, fn(q)) if fn is not None else None

    import pylab
    if show_points:
        plot_points(rho, points)
        pylab.figure()
    plot_calc(r, Pr, q, Iq, theory=theory, title=title, Iq_avg=Iq_avg)
    pylab.gcf().canvas.manager.set_window_title(title)
    pylab.show()

def check_shape_2d(title, shape, fn=None, view=(0, 0, 0), show_points=False,
                   mesh=100, qmax=1.0, samples=5000):
    rho_solvent = 0
    #qx = np.linspace(0.0, qmax, mesh)
    #qy = np.linspace(0.0, qmax, mesh)
    qx = np.linspace(-qmax, qmax, mesh)
    qy = np.linspace(-qmax, qmax, mesh)
    Qx, Qy = np.meshgrid(qx, qy)
    t0 = timer()
    theory = fn(Qx, Qy, view) if fn is not None else None
    print("calc theory time", timer() - t0)

    t0 = timer()
    sampling_density = samples / shape.volume
    if False: # point orientation test: rotate shape rather than view
        shape.rotate(*view)
        view = (0, 0, 0)
    rho, points = shape.sample(sampling_density)
    # The volume of each sample is approximately 1/sampling_density, except
    # that the number of points actually sampled may be slightly more or
    # less due to nature of how points are sampled from the shape.
    volume = shape.volume / len(rho)
    print("point generation time", timer() - t0)
    #print("bounding box", points.min(axis=0), points.max(axis=0))
    # Call calculator with a trivial problem to hide compile times.
    calc_Iqxy(Qx[:1,:1], Qy[:1,:1], rho[:1], points[:1])
    t0 = timer()
    Iqxy = calc_Iqxy(Qx, Qy, rho, points, volume=volume, view=view)
    print("calc Iqxy time", timer() - t0)

    # Add floor to limit colorbar range.
    Iqxy += 0.001 * Iqxy.max()
    if theory is not None:
        theory += 0.001 * theory[~np.isnan(theory)].max()

    import pylab
    if show_points:
        plot_points(rho, apply_view(points, view))
        pylab.figure()
    plot_calc_2d(qx, qy, Iqxy, theory=theory, title=title)
    pylab.gcf().canvas.manager.set_window_title(title)

    ## Histogram of point density in the z direction.
    #pylab.figure()
    ##bins = 100  # fixed number of bins
    #limit=25; bins = np.arange(-limit, limit+0.001, 1) # particular limits
    #pylab.hist(points[:,2], bins=bins)

    pylab.show()

def check_shape_mag(title, shape, fn=None, view=(0, 0, 0), show_points=False,
                    mesh=100, qmax=1.0, samples=5000,
                    up_frac_i=0, up_frac_f=0, up_theta=0, up_phi=0):
    rho_solvent = 0
    #qx = np.linspace(0.0, qmax, mesh)
    #qy = np.linspace(0.0, qmax, mesh)
    qx = np.linspace(-qmax, qmax, mesh)
    qy = np.linspace(-qmax, qmax, mesh)
    Qx, Qy = np.meshgrid(qx, qy)
    sampling_density = samples / shape.volume
    t0 = timer()
    rho, rho_m, points = shape.sample_magnetic(sampling_density)
    # The volume of each sample is approximately 1/sampling_density, except
    # that the number of points actually sampled may be slightly more or
    # less due to nature of how points are sampled from the shape.
    volume = shape.volume / len(rho)
    print("point generation time", timer() - t0)
    #print("bounding box", points.min(axis=0), points.max(axis=0))
    t0 = timer()
    Iqxy = calc_Iq_magnetic(Qx, Qy, rho, rho_m, points, volume=volume, view=view,
                            up_frac_i=up_frac_i, up_frac_f=up_frac_f,
                            up_theta=up_theta, up_phi=up_phi)
    print("calc I_mag time", timer() - t0)
    t0 = timer()
    theory = fn(Qx, Qy, view) if fn is not None else None
    print("calc theory time", timer() - t0)

    # Add floor to limit colorbar range.
    Iqxy += 0.001 * Iqxy.max()
    if theory is not None:
        theory += 0.001 * theory.max()

    import pylab
    if show_points:
        plot_points(rho, points)
        pylab.figure()
    plot_calc_2d(qx, qy, Iqxy, theory=theory, title=title)
    pylab.gcf().canvas.manager.set_window_title(title)

    ## Histogram of point density in the z direction.
    #pylab.figure()
    ##bins = 100  # fixed number of bins
    #limit=25; bins = np.arange(-limit, limit+0.001, 1) # particular limits
    #pylab.hist(points[:,2], bins=bins)

    pylab.show()


def check_pars(function, pars, name):
    spec = getfullargspec(function)
    incorrect = [p for p in pars if p not in spec.args]
    if any(incorrect):
        print("realspace.py: error: invalid parameters to %s: [%s].  Available paramters are:"
              % (name, ", ".join(incorrect)))
        print("   ", "\n    ".join("%s: %g"%(k,v) for k,v in zip(spec.args, spec.defaults)))
        return False
    return True

def main():
    parser = argparse.ArgumentParser(
        description="Compute scattering from realspace sampling",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
    parser.add_argument('-d', '--dim', type=int, default=1,
                        help='dimension 1 or 2')
    parser.add_argument('-m', '--mesh', type=int, default=100,
                        help='number of mesh points in the computed scattering')
    parser.add_argument('-s', '--samples', type=int, default=5000,
                        help="number of sample points in the Monte Carlo estimate")
    parser.add_argument('-q', '--qmax', type=float, default=0.5,
                        help='max q')
    parser.add_argument('-v', '--view', type=str, default='0,0,0',
                        help='theta,phi,psi angles')
    parser.add_argument('-n', '--lattice', type=str, default='1,1,1',
                        help='lattice size')
    parser.add_argument('-z', '--spacing', type=str, default='2,2,2',
                        help='lattice spacing')
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
    shuffle, rotate = opts.shuffle, opts.rotate
    shape_generator = SHAPE_FUNCTIONS[opts.shape]
    if not check_pars(shape_generator, pars, name=opts.shape):
        return
    shape, fn, fn_xy = shape_generator(**pars)
    if nx > 1 or ny > 1 or nz > 1:
        shape = build_cubic_lattice(shape, nx, ny, nz, dx, dy, dz, shuffle, rotate)
    title = "%s(%s)" % (opts.shape, " ".join(opts.pars))
    if shape.is_magnetic:
        view = tuple(float(v) for v in opts.view.split(','))
        up_frac_i, up_frac_f, up_theta, up_phi = shape.spin
        check_shape_mag(title, shape, fn_xy, view=view, show_points=opts.plot,
                       mesh=opts.mesh, qmax=opts.qmax, samples=opts.samples,
                       up_frac_i=up_frac_i, up_frac_f=up_frac_f, up_theta=up_theta, up_phi=up_phi,
                       )
    elif opts.dim == 1:
        check_shape(title, shape, fn, show_points=opts.plot,
                    mesh=opts.mesh, qmax=opts.qmax, samples=opts.samples)
    else:
        view = tuple(float(v) for v in opts.view.split(','))
        check_shape_2d(title, shape, fn_xy, view=view, show_points=opts.plot,
                       mesh=opts.mesh, qmax=opts.qmax, samples=opts.samples)


if __name__ == "__main__":
    # Make sure sasmodels in on the path
    try:
        import sasmodels
    except ImportError:
        import sys
        from os.path import realpath, dirname, join as joinpath
        sys.path.insert(0, dirname(dirname(realpath(__file__))))
    main()
