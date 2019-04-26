#!/usr/bin/env python
"""
Asymmetric shape integration

Usage:

    explore/asymint.py [MODEL] [q-value]

Computes the numerical integral over theta and phi of the given model at a
single point q using different algorithms or the same algorithm with different
precision.  It also displays a 2-D image of the theta-phi surface that is
being integrated.

The available models are:

    triaxial_ellipsoid, parallelpiped, paracrystal, cylinder, sphere

Cylinder and sphere are included as simple checks on the integration
algorithms. Cylinder is better investigated using 1-D integration methods in
explore/symint.py.  Sphere has an easily computed analytic value which is
identical for all theta-phi for a given q, so it is useful for checking
that the normalization constants are correct for the different algorithms.
"""

from __future__ import print_function, division

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import warnings

import numpy as np
import mpmath as mp
from numpy import pi, sin, cos, sqrt, exp, expm1, degrees, log10, arccos
from numpy.polynomial.legendre import leggauss
from scipy.integrate import dblquad, simps, romb, romberg
import pylab

import sasmodels.special as sp

DTYPE = 'd'

class MPenv:
    sqrt = staticmethod(mp.sqrt)
    exp = staticmethod(mp.exp)
    expm1 = staticmethod(mp.expm1)
    cos = staticmethod(mp.cos)
    sin = staticmethod(mp.sin)
    tan = staticmethod(mp.tan)
    @staticmethod
    def sas_3j1x_x(x):
        return 3*(mp.sin(x)/x - mp.cos(x))/(x*x)
    @staticmethod
    def sas_2J1x_x(x):
        return 2*mp.j1(x)/x
    @staticmethod
    def sas_sinx_x(x):
        return mp.sin(x)/x
    pi = mp.pi
    mpf = staticmethod(mp.mpf)

class NPenv:
    sqrt = staticmethod(np.sqrt)
    exp = staticmethod(np.exp)
    expm1 = staticmethod(np.expm1)
    cos = staticmethod(np.cos)
    sin = staticmethod(np.sin)
    tan = staticmethod(np.tan)
    sas_3j1x_x = staticmethod(sp.sas_3j1x_x)
    sas_2J1x_x = staticmethod(sp.sas_2J1x_x)
    sas_sinx_x = staticmethod(sp.sas_sinx_x)
    pi = np.pi
    #mpf = staticmethod(float)
    mpf = staticmethod(lambda x: np.array(x, DTYPE))

SLD = 3
SLD_SOLVENT = 6
CONTRAST = SLD - SLD_SOLVENT

# Carefully code models so that mpmath will use full precision.  That means:
#    * wrap inputs in env.mpf
#    * don't use floating point constants, only integers
#    * for division, make sure the numerator or denominator is env.mpf
#    * use env.pi, env.sas_sinx_x, etc. for functions
def make_parallelepiped(a, b, c, env=NPenv):
    a, b, c = env.mpf(a), env.mpf(b), env.mpf(c)
    def Fq(qa, qb, qc):
        siA = env.sas_sinx_x(a*qa/2)
        siB = env.sas_sinx_x(b*qb/2)
        siC = env.sas_sinx_x(c*qc/2)
        return siA * siB * siC
    Fq.__doc__ = "parallelepiped a=%g, b=%g c=%g"%(a, b, c)
    volume = a*b*c
    norm = CONTRAST**2*volume/10000
    return norm, Fq

def make_core_shell_parallelepiped(a, b, c, da, db, dc, slda, sldb, sldc, env=NPenv):
    overlapping = False
    a, b, c = env.mpf(a), env.mpf(b), env.mpf(c)
    da, db, dc = env.mpf(da), env.mpf(db), env.mpf(dc)
    slda, sldb, sldc = env.mpf(slda), env.mpf(sldb), env.mpf(sldc)
    dr0 = CONTRAST
    drA, drB, drC = slda-SLD_SOLVENT, sldb-SLD_SOLVENT, sldc-SLD_SOLVENT
    tA, tB, tC = a + 2*da, b + 2*db, c + 2*dc
    def Fq(qa, qb, qc):
        siA = a*env.sas_sinx_x(a*qa/2)
        siB = b*env.sas_sinx_x(b*qb/2)
        siC = c*env.sas_sinx_x(c*qc/2)
        siAt = tA*env.sas_sinx_x(tA*qa/2)
        siBt = tB*env.sas_sinx_x(tB*qb/2)
        siCt = tC*env.sas_sinx_x(tC*qc/2)
        if overlapping:
            return (dr0*siA*siB*siC
                    + drA*(siAt-siA)*siB*siC
                    + drB*siAt*(siBt-siB)*siC
                    + drC*siAt*siBt*(siCt-siC))
        else:
            return (dr0*siA*siB*siC
                    + drA*(siAt-siA)*siB*siC
                    + drB*siA*(siBt-siB)*siC
                    + drC*siA*siB*(siCt-siC))
    Fq.__doc__ = "core-shell parallelepiped a=%g, b=%g c=%g"%(a, b, c)
    if overlapping:
        volume = a*b*c + 2*da*b*c + 2*tA*db*c + 2*tA*tB*dc
    else:
        volume = a*b*c + 2*da*b*c + 2*a*db*c + 2*a*b*dc
    norm = 1/(volume*10000)
    return norm, Fq

def make_triaxial_ellipsoid(a, b, c, env=NPenv):
    a, b, c = env.mpf(a), env.mpf(b), env.mpf(c)
    def Fq(qa, qb, qc):
        qr = env.sqrt((a*qa)**2 + (b*qb)**2 + (c*qc)**2)
        return env.sas_3j1x_x(qr)
    Fq.__doc__ = "triaxial ellipsoid minor=%g, major=%g polar=%g"%(a, b, c)
    volume = 4*env.pi*a*b*c/3
    norm = CONTRAST**2*volume/10000
    return norm, Fq

def make_cylinder(radius, length, env=NPenv):
    radius, length = env.mpf(radius), env.mpf(length)
    def Fq(qa, qb, qc):
        qab = env.sqrt(qa**2 + qb**2)
        return env.sas_2J1x_x(qab*radius) * env.sas_sinx_x((qc*length)/2)
    Fq.__doc__ = "cylinder radius=%g, length=%g"%(radius, length)
    volume = env.pi*radius**2*length
    norm = CONTRAST**2*volume/10000
    return norm, Fq

def make_sphere(radius, env=NPenv):
    radius = env.mpf(radius)
    def Fq(qa, qb, qc):
        q = env.sqrt(qa**2 + qb**2 + qc**2)
        return env.sas_3j1x_x(q*radius)
    Fq.__doc__ = "sphere radius=%g"%(radius, )
    volume = 4*pi*radius**3
    norm = CONTRAST**2*volume/10000
    return norm, Fq

def make_paracrystal(radius, dnn, d_factor, lattice='bcc', env=NPenv):
    radius, dnn, d_factor = env.mpf(radius), env.mpf(dnn), env.mpf(d_factor)
    def sc(qa, qb, qc):
        return qa, qb, qc
    def bcc(qa, qb, qc):
        a1 = (+qa + qb + qc)/2
        a2 = (-qa - qb + qc)/2
        a3 = (-qa + qb - qc)/2
        return a1, a2, a3
    def fcc(qa, qb, qc):
        a1 = ( 0 + qb + qc)/2
        a2 = (-qa + 0 + qc)/2
        a3 = (-qa + qb + 0)/2
        return a1, a2, a3
    lattice_fn = {'sc': sc, 'bcc': bcc, 'fcc': fcc}[lattice]
    radius, dnn, d_factor = env.mpf(radius), env.mpf(dnn), env.mpf(d_factor)
    def Fq(qa, qb, qc):
        a1, a2, a3 = lattice_fn(qa, qb, qc)
        # Note: paper says that different directions can have different
        # distoration factors.  Easy enough to add to the code.
        arg = -(dnn*d_factor)**2*(a1**2 + a2**2 + a3**2)/2
        exp_arg = env.exp(arg)
        den = [((exp_arg - 2*env.cos(dnn*a))*exp_arg + 1) for a in (a1, a2, a3)]
        Sq = -env.expm1(2*arg)**3/(den[0]*den[1]*den[2])

        q = env.sqrt(qa**2 + qb**2 + qc**2)
        Fq = env.sas_3j1x_x(q*radius)
        # the caller computes F(q)**2, but we need it to compute S(q)*F(q)**2
        return env.sqrt(Sq)*Fq
    Fq.__doc__ = "%s paracrystal a=%g da=%g r=%g"%(lattice, dnn, d_factor, radius)
    def sphere_volume(r): return 4*env.pi*r**3/3
    Vf = {
        'sc': sphere_volume(radius/dnn),
        'bcc': 2*sphere_volume(env.sqrt(3)/2*radius/dnn),
        'fcc': 4*sphere_volume(1/env.sqrt(2)*radius/dnn),
    }[lattice]
    volume = sphere_volume(radius)
    norm = CONTRAST**2*volume/10000*Vf
    return norm, Fq

NORM = 1.0  # type: float
KERNEL = None  # type: CALLABLE[[ndarray, ndarray, ndarray], ndarray]
NORM_MP = 1  # type: mpf
KERNEL = None  # type: CALLABLE[[mpf, mpf, mpf], mpf]

SHAPES = [
    'sphere',
    'cylinder',
    'triaxial_ellipsoid',
    'parallelepiped',
    'core_shell_parallelepiped',
    'fcc_paracrystal',
    'bcc_paracrystal',
    'sc_paracrystal',
]
def build_shape(shape, **pars):
    global NORM, KERNEL
    global NORM_MP, KERNEL_MP

    # Note: using integer or string defaults for the sake of mpf
    if shape == 'sphere':
        RADIUS = pars.get('radius', 50)
        NORM, KERNEL = make_sphere(radius=RADIUS)
        NORM_MP, KERNEL_MP = make_sphere(radius=RADIUS, env=MPenv)
    elif shape == 'cylinder':
        #RADIUS, LENGTH = 10, 100000
        RADIUS = pars.get('radius', 10)
        LENGTH = pars.get('radius', 300)
        NORM, KERNEL = make_cylinder(radius=RADIUS, length=LENGTH)
        NORM_MP, KERNEL_MP = make_cylinder(radius=RADIUS, length=LENGTH, env=MPenv)
    elif shape == 'triaxial_ellipsoid':
        #A, B, C = 4450, 14000, 47
        A = pars.get('a', 445)
        B = pars.get('b', 140)
        C = pars.get('c', 47)
        NORM, KERNEL = make_triaxial_ellipsoid(A, B, C)
        NORM_MP, KERNEL_MP = make_triaxial_ellipsoid(A, B, C, env=MPenv)
    elif shape == 'parallelepiped':
        #A, B, C = 4450, 14000, 47
        A = pars.get('a', 445)
        B = pars.get('b', 140)
        C = pars.get('c', 47)
        NORM, KERNEL = make_parallelepiped(A, B, C)
        NORM_MP, KERNEL_MP = make_parallelepiped(A, B, C, env=MPenv)
    elif shape == 'core_shell_parallelepiped':
        #A, B, C = 4450, 14000, 47
        #A, B, C = 445, 140, 47  # integer for the sake of mpf
        A = pars.get('a', 114)
        B = pars.get('b', 1380)
        C = pars.get('c', 6800)
        DA = pars.get('da', 21)
        DB = pars.get('db', 58)
        DC = pars.get('dc', 2300)
        SLDA = pars.get('slda', "5")
        SLDB = pars.get('sldb', "-0.3")
        SLDC = pars.get('sldc', "11.5")
        ## default parameters from sasmodels
        #A,B,C,DA,DB,DC,SLDA,SLDB,SLDC = 400,75,35,10,10,10,2,4,2
        ## swap A-B-C to C-B-A
        #A, B, C, DA, DB, DC, SLDA, SLDB, SLDC = C, B, A, DC, DB, DA, SLDC, SLDB, SLDA
        #A,B,C,DA,DB,DC,SLDA,SLDB,SLDC = 10,20,30,100,200,300,1,2,3
        #SLD_SOLVENT,CONTRAST = 0, 4
        if 1: # C shortest
            B, C = C, B
            DB, DC = DC, DB
            SLDB, SLDC = SLDC, SLDB
        elif 0: # C longest
            A, C = C, A
            DA, DC = DC, DA
            SLDA, SLDC = SLDC, SLDA
        #NORM, KERNEL = make_core_shell_parallelepiped(A, B, C, DA, DB, DC, SLDA, SLDB, SLDC)
        NORM, KERNEL = make_core_shell_parallelepiped(A, B, C, DA, DB, DC, SLDA, SLDB, SLDC)
        NORM_MP, KERNEL_MP = make_core_shell_parallelepiped(A, B, C, DA, DB, DC, SLDA, SLDB, SLDC, env=MPenv)
    elif shape.endswith('paracrystal'):
        LATTICE, _ = shape.split('_')
        DNN = pars.get('dnn', 220)
        D_FACTOR = pars.get('d_factor', '0.06')
        RADIUS = pars.get('radius', 40)
        NORM, KERNEL = make_paracrystal(
            radius=RADIUS, dnn=DNN, d_factor=D_FACTOR, lattice=LATTICE)
        NORM_MP, KERNEL_MP = make_paracrystal(
            radius=RADIUS, dnn=DNN, d_factor=D_FACTOR, lattice=LATTICE, env=MPenv)
    else:
        raise ValueError("Unknown shape %r"%shape)

# Note: hardcoded in mp_quad
THETA_LOW, THETA_HIGH = 0, pi
PHI_LOW, PHI_HIGH = 0, 2*pi
SCALE = 1

# mathematica code for triaxial_ellipsoid (untested)
_ = """
R[theta_, phi_, a_, b_, c_] := Sqrt[(a Sin[theta]Cos[phi])^2 + (b Sin[theta]Sin[phi])^2 + (c Cos[theta])^2]
Sphere[q_, r_] := 3 SphericalBesselJ[q r]/(q r)
V[a_, b_, c_] := 4/3 pi a b c
Norm[sld_, solvent_, a_, b_, c_] := V[a, b, c] (solvent - sld)^2
F[q_, theta_, phi_, a_, b_, c_] := Sphere[q, R[theta, phi, a, b, c]]
I[q_, sld_, solvent_, a_, b_, c_] := Norm[sld, solvent, a, b, c]/(4 pi) Integrate[F[q, theta, phi, a, b, c]^2 Sin[theta], {phi, 0, 2 pi}, {theta, 0, pi}]
I[6/10^3, 63/10, 3, 445, 140, 47]
"""

# 2D integration functions
def mp_quad_2d(q):
    evals = [0]
    def integrand(theta, phi):
        evals[0] += 1
        qab = q*mp.sin(theta)
        qa = qab*mp.cos(phi)
        qb = qab*mp.sin(phi)
        qc = q*mp.cos(theta)
        Zq = KERNEL_MP(qa, qb, qc)**2
        return Zq*mp.sin(theta)
    ans = mp.quad(integrand, (0, mp.pi), (0, 2*mp.pi))
    Iq = NORM_MP*ans/(4*mp.pi)
    return evals[0], Iq

def kernel_2d(q, theta, phi):
    """
    S(q) kernel for paracrystal forms.
    """
    qab = q*sin(theta)
    qa = qab*cos(phi)
    qb = qab*sin(phi)
    qc = q*cos(theta)
    return NORM*KERNEL(qa, qb, qc)**2

def scipy_dblquad_2d(q):
    """
    Compute the integral using scipy dblquad.  This gets the correct answer
    eventually, but it is slow.
    """
    evals = [0]
    def integrand(phi, theta):
        evals[0] += 1
        Zq = kernel_2d(q, theta=theta, phi=phi)
        return Zq*sin(theta)
    ans = dblquad(integrand, THETA_LOW, THETA_HIGH, lambda x: PHI_LOW, lambda x: PHI_HIGH)[0]
    return evals[0], ans*SCALE/(4*pi)


def scipy_romberg_2d(q):
    """
    Compute the integral using romberg integration.  This function does not
    complete in a reasonable time.  No idea if it is accurate.
    """
    evals = [0]
    def inner(phi, theta):
        evals[0] += 1
        return kernel_2d(q, theta=theta, phi=phi)
    def outer(theta):
        Zq = romberg(inner, PHI_LOW, PHI_HIGH, divmax=100, args=(theta,))
        return Zq*sin(theta)
    ans = romberg(outer, THETA_LOW, THETA_HIGH, divmax=100)
    return evals[0], ans*SCALE/(4*pi)


def semi_romberg_2d(q, n=100):
    """
    Use 1D romberg integration in phi and regular simpsons rule in theta.
    """
    evals = [0]
    def inner(phi, theta):
        evals[0] += 1
        return kernel_2d(q, theta=theta, phi=phi)
    theta = np.linspace(THETA_LOW, THETA_HIGH, n)
    Zq = [romberg(inner, PHI_LOW, PHI_HIGH, divmax=100, args=(t,)) for t in theta]
    ans = simps(np.array(Zq)*sin(theta), dx=theta[1]-theta[0])
    return evals[0], ans*SCALE/(4*pi)

def gauss_quad_2d(q, n=150):
    """
    Compute the integral using gaussian quadrature for n = 20, 76 or 150.
    """
    z, w = leggauss(n)
    theta = (THETA_HIGH-THETA_LOW)*(z + 1)/2 + THETA_LOW
    phi = (PHI_HIGH-PHI_LOW)*(z + 1)/2 + PHI_LOW
    Atheta, Aphi = np.meshgrid(theta, phi)
    Aw = w[None, :] * w[:, None]
    sin_theta = abs(sin(Atheta))
    Zq = kernel_2d(q=q, theta=Atheta, phi=Aphi)
    # change from [-1,1] x [-1,1] range to [0, pi] x [0, 2 pi] range
    dxdy_stretch = (THETA_HIGH-THETA_LOW)/2 * (PHI_HIGH-PHI_LOW)/2
    Iq = np.sum(Zq*Aw*sin_theta)*SCALE/(4*pi) * dxdy_stretch
    return n**2, Iq

def gauss_quad_usub(q, n=150, dtype=DTYPE):
    """
    Compute the integral using gaussian quadrature for n = 20, 76 or 150.

    Use *u = sin theta* substitution, and restrict integration over a single
    quadrant for shapes that are mirror symmetric about AB, AC and BC planes.

    Note that this doesn't work for fcc/bcc paracrystals, which instead step
    over the entire 4 pi surface uniformly in theta-phi.
    """
    z, w = leggauss(n)
    cos_theta = 0.5 * (z + 1)
    theta = arccos(cos_theta)
    phi = pi/2*(0.5 * (z + 1))
    Atheta, Aphi = np.meshgrid(theta, phi)
    Aw = w[None, :] * w[:, None]
    q, Atheta, Aphi, Aw = [np.asarray(v, dtype=dtype) for v in (q, Atheta, Aphi, Aw)]
    Zq = kernel_2d(q=q, theta=Atheta, phi=Aphi)
    Iq = np.sum(Zq*Aw)*0.25
    return n**2, Iq

def gridded_2d(q, n=300):
    """
    Compute the integral on a regular grid using rectangular, trapezoidal,
    simpsons, and romberg integration.  Romberg integration requires that
    the grid be of size n = 2**k + 1.
    """
    theta = np.linspace(THETA_LOW, THETA_HIGH, n)
    phi = np.linspace(PHI_LOW, PHI_HIGH, n)
    Atheta, Aphi = np.meshgrid(theta, phi)
    Zq = kernel_2d(q=q, theta=Atheta, phi=Aphi)
    Zq *= abs(sin(Atheta))
    dx, dy = theta[1]-theta[0], phi[1]-phi[0]
    print("rect-%d"%n, n**2, np.sum(Zq)*dx*dy*SCALE/(4*pi))
    print("trapz-%d"%n, n**2, np.trapz(np.trapz(Zq, dx=dx), dx=dy)*SCALE/(4*pi))
    print("simpson-%d"%n, n**2, simps(simps(Zq, dx=dx), dx=dy)*SCALE/(4*pi))
    print("romb-%d"%n, n**2, romb(romb(Zq, dx=dx), dx=dy)*SCALE/(4*pi))

def quadpy_method(q, rule):
    """
    Use *rule*="name:index" where name and index are chosen from below.

    Available rule names and the corresponding indices::

        AlbrechtCollatz: [1-5]
        BazantOh: 9, 11, 13
        HeoXu: 13, 15, 17, 19-[1-2], 21-[1-6], 23-[1-3], 25-[1-2], 27-[1-3],
            29, 31, 33, 35, 37, 39-[1-2]
        FliegeMaier: 4, 9, 16, 25
        Lebedev: 3[a-c], 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35,
            41, 47, 53, 59, 65, 71, 77 83, 89, 95, 101, 107, 113, 119, 125, 131
        McLaren: [1-10]
        Stroud: U3 3-1, U3 5-[1-5], U3 7-[1-2], U3 8-1, U3 9-[1-3],
            U3 11-[1-3], U3 14-1
    """
    try:
        import quadpy
    except ImportError:
        warnings.warn("use 'pip install quadpy' to enable quadpy.sphere tests")
        return

    from quadpy.sphere import (AlbrechtCollatz, BazantOh, HeoXu,
        FliegeMaier, Lebedev, McLaren, Stroud, integrate_spherical)
    RULES = {
        'AlbrechtCollatz': AlbrechtCollatz,
        'BazantOh': BazantOh,
        'HeoXu': HeoXu,
        'FliegeMaier': FliegeMaier,
        'Lebedev': Lebedev,
        'McLaren': McLaren,
        'Stroud': Stroud,
    }
    int_index = 'AlbrechtCollatz', 'McLaren'

    rule_name, rule_index = rule.split(':')
    index = int(rule_index) if rule_name in int_index else rule_index
    rule_obj = RULES[rule_name](index)
    fn = lambda azimuthal, polar: kernel_2d(q=q, theta=polar, phi=azimuthal)
    Iq = integrate_spherical(fn, rule=rule_obj)/(4*pi)
    print("%s degree=%d points=%s => %.15g"
          % (rule, rule_obj.degree, len(rule_obj.points), Iq))

def plot_2d(q, n=300):
    """
    Plot the 2D surface that needs to be integrated in order to compute
    the BCC S(q) at a particular q, dnn and d_factor.  *n* is the number
    of points in the grid.
    """
    theta = np.linspace(THETA_LOW, THETA_HIGH, n)
    phi = np.linspace(PHI_LOW, PHI_HIGH, n)
    Atheta, Aphi = np.meshgrid(theta, phi)
    Zq = kernel_2d(q=q, theta=Atheta, phi=Aphi)
    #Zq *= abs(sin(Atheta))
    pylab.pcolor(degrees(theta), degrees(phi), log10(np.fmax(Zq, 1.e-6)))
    pylab.axis('tight')
    pylab.title("%s I(q,t) sin(t) for q=%g" % (KERNEL.__doc__, q))
    pylab.xlabel("theta (degrees)")
    pylab.ylabel("phi (degrees)")
    cbar = pylab.colorbar()
    cbar.set_label('log10 S(q)')
    pylab.show()

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="asymmetric integration explorer",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
    parser.add_argument('-s', '--shape', choices=SHAPES,
                        default='parallelepiped',
                        help='oriented shape')
    parser.add_argument('-q', '--q_value', type=str, default='0.005',
                        help='Q value to evaluate')
    parser.add_argument('pars', type=str, nargs='*', default=[],
                        help='p=val for p in shape parameters')
    opts = parser.parse_args()
    pars = {k: v for par in opts.pars for k, v in [par.split('=')]}
    build_shape(opts.shape, **pars)

    Q = float(opts.q_value)
    if opts.shape == 'sphere':
        print("exact", NORM*sp.sas_3j1x_x(Q*RADIUS)**2)

    # Methods from quadpy, if quadpy is available
    #  AlbrechtCollatz: [1-5]
    #  BazantOh: 9, 11, 13
    #  HeoXu: 13, 15, 17, 19-[1-2], 21-[1-6], 23-[1-3], 25-[1-2], 27-[1-3],
    #     29, 31, 33, 35, 37, 39-[1-2]
    #  FliegeMaier: 4, 9, 16, 25
    #  Lebedev: 3[a-c], 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35,
    #     41, 47, 53, 59, 65, 71, 77 83, 89, 95, 101, 107, 113, 119, 125, 131
    #  McLaren: [1-10]
    #  Stroud: U3 3-1, U3 5-[1-5], U3 7-[1-2], U3 8-1, U3 9-[1-3],
    #     U3 11-[1-3], U3 14-1
    quadpy_method(Q, "AlbrechtCollatz:5")
    quadpy_method(Q, "HeoXu:39-2")
    quadpy_method(Q, "FliegeMaier:25")
    quadpy_method(Q, "Lebedev:19")
    quadpy_method(Q, "Lebedev:131")
    quadpy_method(Q, "McLaren:10")
    quadpy_method(Q, "Stroud:U3 14-1")

    print("gauss-20 points=%d => %.15g" % gauss_quad_2d(Q, n=20))
    print("gauss-76 points=%d => %.15g" % gauss_quad_2d(Q, n=76))
    print("gauss-150 points=%d => %.15g" % gauss_quad_2d(Q, n=150))
    print("gauss-500 points=%d => %.15g" % gauss_quad_2d(Q, n=500))
    print("gauss-1025 points=%d => %.15g" % gauss_quad_2d(Q, n=1025))
    print("gauss-2049 points=%d => %.15g" % gauss_quad_2d(Q, n=2049))
    print("gauss-20 usub points=%d => %.15g" % gauss_quad_usub(Q, n=20))
    print("gauss-76 usub points=%d => %.15g" % gauss_quad_usub(Q, n=76))
    print("gauss-150 usub points=%d => %.15g" % gauss_quad_usub(Q, n=150))

    #gridded_2d(Q, n=2**8+1)
    gridded_2d(Q, n=2**10+1)
    #gridded_2d(Q, n=2**12+1)
    #gridded_2d(Q, n=2**15+1)
    # adaptive forms on models for which the calculations are fast enough
    SLOW_SHAPES = {
        'fcc_paracrystal', 'bcc_paracrystal', 'sc_paracrystal',
        'core_shell_parallelepiped',
    }
    if opts.shape not in SLOW_SHAPES:
        print("dblquad", *scipy_dblquad_2d(Q))
        print("semi-romberg-100", *semi_romberg_2d(Q, n=100))
        print("romberg", *scipy_romberg_2d(Q))
        with mp.workprec(100):
            print("mpmath", *mp_quad_2d(mp.mpf(opts.q_value)))
    plot_2d(Q, n=200)

if __name__ == "__main__":
    main()
