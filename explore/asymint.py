"""
Asymmetric shape integration
"""

from __future__ import print_function, division

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np
import mpmath as mp
from numpy import pi, sin, cos, sqrt, exp, expm1, degrees, log10
from numpy.polynomial.legendre import leggauss
from scipy.integrate import dblquad, simps, romb, romberg
import pylab

from sasmodels.special import square
from sasmodels.special import Gauss20Wt, Gauss20Z
from sasmodels.special import Gauss76Wt, Gauss76Z
from sasmodels.special import Gauss150Wt, Gauss150Z
from sasmodels.special import sas_2J1x_x, sas_sinx_x, sas_3j1x_x

def mp_3j1x_x(x):
    return 3*(mp.sin(x)/x - mp.cos(x))/(x*x)
def mp_2J1x_x(x):
    return 2*mp.j1(x)/x
def mp_sinx_x(x):
    return mp.sin(x)/x

SLD = 3
SLD_SOLVENT = 6
CONTRAST = SLD - SLD_SOLVENT

def make_parallelepiped(a, b, c):
    def Fq(qa, qb, qc):
        siA = sas_sinx_x(0.5*a*qa)
        siB = sas_sinx_x(0.5*b*qb)
        siC = sas_sinx_x(0.5*c*qc)
        return siA * siB * siC
    volume = a*b*c
    norm = volume*CONTRAST**2/10**4
    return norm, Fq

def make_parallelepiped_mp(a, b, c):
    a, b, c = mp.mpf(a), mp.mpf(b), mp.mpf(c)
    def Fq(qa, qb, qc):
        siA = mp_sinx_x(a*qa/2)
        siB = mp_sinx_x(b*qb/2)
        siC = mp_sinx_x(c*qc/2)
        return siA * siB * siC
    volume = a*b*c
    norm = (volume*CONTRAST**2)/10000 # mpf since volume=a*b*c is mpf
    return norm, Fq

def make_triellip(a, b, c):
    def Fq(qa, qb, qc):
        qr = sqrt((a*qa)**2 + (b*qb)**2 + (c*qc)**2)
        return sas_3j1x_x(qr)
    volume = 4*pi*a*b*c/3
    norm = volume*CONTRAST**2/10**4
    return norm, Fq

def make_triellip_mp(a, b, c):
    a, b, c = mp.mpf(a), mp.mpf(b), mp.mpf(c)
    def Fq(qa, qb, qc):
        qr = mp.sqrt((a*qa)**2 + (b*qb)**2 + (c*qc)**2)
        return mp_3j1x_x(qr)
    volume = (4*mp.pi*a*b*c)/3
    norm = (volume*CONTRAST**2)/10000  # mpf since mp.pi is mpf
    return norm, Fq

def make_cylinder(radius, length):
    def Fq(qa, qb, qc):
        qab = sqrt(qa**2 + qb**2)
        return sas_2J1x_x(qab*radius) * sas_sinx_x(qc*0.5*length)
    volume = pi*radius**2*length
    norm = volume*CONTRAST**2/10**4
    return norm, Fq

def make_cylinder_mp(radius, length):
    radius, length = mp.mpf(radius), mp.mpf(length)
    def Fq(qa, qb, qc):
        qab = mp.sqrt(qa**2 + qb**2)
        return mp_2J1x_x(qab*radius) * mp_sinx_x((qc*length)/2)
    volume = mp.pi*radius**2*length
    norm = (volume*CONTRAST**2)/10000  # mpf since mp.pi is mpf
    return norm, Fq

def make_sphere(radius):
    def Fq(qa, qb, qc):
        q = sqrt(qa**2 + qb**2 + qc**2)
        return sas_3j1x_x(q*radius)
    volume = 4*pi*radius**3/3
    norm = volume*CONTRAST**2/10**4
    return norm, Fq

def make_sphere_mp(radius):
    radius = mp.mpf(radius)
    def Fq(qa, qb, qc):
        q = mp.sqrt(qa**2 + qb**2 + qc**2)
        return mp_3j1x_x(q*radius)
    volume = (4*mp.pi*radius**3)/3
    norm = (volume*CONTRAST**2)/10000  # mpf since mp.pi is mpf
    return norm, Fq

shape = 'parallelepiped'
#shape = 'triellip'
#shape = 'sphere'
#shape = 'cylinder'
if shape == 'cylinder':
    #RADIUS, LENGTH = 10, 100000
    RADIUS, LENGTH = 10, 300  # integer for the sake of mpf
    NORM, KERNEL = make_cylinder(radius=RADIUS, length=LENGTH)
    NORM_MP, KERNEL_MP = make_cylinder_mp(radius=RADIUS, length=LENGTH)
elif shape == 'triellip':
    #A, B, C = 4450, 14000, 47
    A, B, C = 445, 140, 47  # integer for the sake of mpf
    NORM, KERNEL = make_triellip(A, B, C)
    NORM_MP, KERNEL_MP = make_triellip_mp(A, B, C)
elif shape == 'parallelepiped':
    #A, B, C = 4450, 14000, 47
    A, B, C = 445, 140, 47  # integer for the sake of mpf
    NORM, KERNEL = make_parallelepiped(A, B, C)
    NORM_MP, KERNEL_MP = make_parallelepiped_mp(A, B, C)
elif shape == 'sphere':
    RADIUS = 50  # integer for the sake of mpf
    NORM, KERNEL = make_sphere(radius=RADIUS)
    NORM_MP, KERNEL_MP = make_sphere_mp(radius=RADIUS)
else:
    raise ValueError("Unknown shape %r"%shape)

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


def mp_quad(q, shape):
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

def kernel(q, theta, phi):
    """
    S(q) kernel for paracrystal forms.
    """
    qab = q*sin(theta)
    qa = qab*cos(phi)
    qb = qab*sin(phi)
    qc = q*cos(theta)
    return NORM*KERNEL(qa, qb, qc)**2

def scipy_dblquad(q):
    """
    Compute the integral using scipy dblquad.  This gets the correct answer
    eventually, but it is slow.
    """
    evals = [0]
    def integrand(phi, theta):
        evals[0] += 1
        Zq = kernel(q, theta=theta, phi=phi)
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
        return kernel(q, theta=theta, phi=phi)
    def outer(theta):
        Zq = romberg(inner, PHI_LOW, PHI_HIGH, divmax=100, args=(theta,))
        return Zq*sin(theta)
    ans = romberg(outer, THETA_LOW, THETA_HIGH, divmax=100)
    return evals[0], ans*SCALE/(4*pi)


def semi_romberg(q, n=100):
    """
    Use 1D romberg integration in phi and regular simpsons rule in theta.
    """
    evals = [0]
    def inner(phi, theta):
        evals[0] += 1
        return kernel(q, theta=theta, phi=phi)
    theta = np.linspace(THETA_LOW, THETA_HIGH, n)
    Zq = [romberg(inner, PHI_LOW, PHI_HIGH, divmax=100, args=(t,)) for t in theta]
    ans = simps(np.array(Zq)*sin(theta), dx=theta[1]-theta[0])
    return evals[0], ans*SCALE/(4*pi)

def gauss_quad(q, n=150):
    """
    Compute the integral using gaussian quadrature for n = 20, 76 or 150.
    """
    if n == 20:
        z, w = Gauss20Z, Gauss20Wt
    elif n == 76:
        z, w = Gauss76Z, Gauss76Wt
    elif n == 150:
        z, w = Gauss150Z, Gauss150Wt
    else:
        z, w = leggauss(n)
    theta = (THETA_HIGH-THETA_LOW)*(z + 1)/2 + THETA_LOW
    phi = (PHI_HIGH-PHI_LOW)*(z + 1)/2 + PHI_LOW
    Atheta, Aphi = np.meshgrid(theta, phi)
    Aw = w[None, :] * w[:, None]
    sin_theta = abs(sin(Atheta))
    Zq = kernel(q=q, theta=Atheta, phi=Aphi)
    # change from [-1,1] x [-1,1] range to [0, pi] x [0, 2 pi] range
    dxdy_stretch = (THETA_HIGH-THETA_LOW)/2 * (PHI_HIGH-PHI_LOW)/2
    Iq = np.sum(Zq*Aw*sin_theta)*SCALE/(4*pi) * dxdy_stretch
    return n**2, Iq


def gridded_integrals(q, n=300):
    """
    Compute the integral on a regular grid using rectangular, trapezoidal,
    simpsons, and romberg integration.  Romberg integration requires that
    the grid be of size n = 2**k + 1.
    """
    theta = np.linspace(THETA_LOW, THETA_HIGH, n)
    phi = np.linspace(PHI_LOW, PHI_HIGH, n)
    Atheta, Aphi = np.meshgrid(theta, phi)
    Zq = kernel(q=q, theta=Atheta, phi=Aphi)
    Zq *= abs(sin(Atheta))
    dx, dy = theta[1]-theta[0], phi[1]-phi[0]
    print("rect-%d"%n, n**2, np.sum(Zq)*dx*dy*SCALE/(4*pi))
    print("trapz-%d"%n, n**2, np.trapz(np.trapz(Zq, dx=dx), dx=dy)*SCALE/(4*pi))
    print("simpson-%d"%n, n**2, simps(simps(Zq, dx=dx), dx=dy)*SCALE/(4*pi))
    print("romb-%d"%n, n**2, romb(romb(Zq, dx=dx), dx=dy)*SCALE/(4*pi))

def plot(q, n=300):
    """
    Plot the 2D surface that needs to be integrated in order to compute
    the BCC S(q) at a particular q, dnn and d_factor.  *n* is the number
    of points in the grid.
    """
    theta = np.linspace(THETA_LOW, THETA_HIGH, n)
    phi = np.linspace(PHI_LOW, PHI_HIGH, n)
    Atheta, Aphi = np.meshgrid(theta, phi)
    Zq = kernel(q=q, theta=Atheta, phi=Aphi)
    #Zq *= abs(sin(Atheta))
    pylab.pcolor(degrees(theta), degrees(phi), log10(np.fmax(Zq, 1.e-6)))
    pylab.axis('tight')
    pylab.title("%s Z(q) for q=%g" % (KERNEL.__name__, q))
    pylab.xlabel("theta (degrees)")
    pylab.ylabel("phi (degrees)")
    cbar = pylab.colorbar()
    cbar.set_label('log10 S(q)')
    pylab.show()

if __name__ == "__main__":
    Qstr = '0.005'
    #Qstr = '0.8'
    #Qstr = '0.0003'
    Q = float(Qstr)
    if shape == 'sphere':
        print("exact", NORM*sas_3j1x_x(Q*RADIUS)**2)
    print("gauss-20", *gauss_quad(Q, n=20))
    print("gauss-76", *gauss_quad(Q, n=76))
    print("gauss-150", *gauss_quad(Q, n=150))
    print("gauss-500", *gauss_quad(Q, n=500))
    print("dblquad", *scipy_dblquad(Q))
    print("semi-romberg-100", *semi_romberg(Q, n=100))
    print("romberg", *scipy_romberg_2d(Q))
    #gridded_integrals(Q, n=2**8+1)
    gridded_integrals(Q, n=2**10+1)
    #gridded_integrals(Q, n=2**13+1)
    #gridded_integrals(Q, n=2**15+1)
    with mp.workprec(100):
        print("mpmath", *mp_quad(mp.mpf(Qstr), shape))
    #plot(Q, n=200)
