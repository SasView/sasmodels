"""
Explore integration of rotationally symmetric shapes
"""

from __future__ import print_function, division

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np
from numpy import pi, sin, cos, sqrt, exp, expm1, degrees, log10
from numpy.polynomial.legendre import leggauss
from scipy.integrate import dblquad, simps, romb, romberg
import pylab

from sasmodels.special import square
from sasmodels.special import Gauss20Wt, Gauss20Z
from sasmodels.special import Gauss76Wt, Gauss76Z
from sasmodels.special import Gauss150Wt, Gauss150Z
from sasmodels.special import sas_2J1x_x, sas_sinx_x, sas_3j1x_x

SLD = 3.0
SLD_SOLVENT = 6
CONTRAST = SLD - SLD_SOLVENT

def make_cylinder(radius, length):
    def cylinder(qab, qc):
        return sas_2J1x_x(qab*radius) * sas_sinx_x(qc*0.5*length)
    cylinder.__doc__ = "cylinder radius=%g, length=%g"%(radius, length)
    volume = pi*radius**2*length
    norm = CONTRAST**2*volume/10000
    return norm, cylinder

def make_long_cylinder(radius, length):
    def long_cylinder(q):
        return norm/q * sas_2J1x_x(q*radius)**2
    long_cylinder.__doc__ = "long cylinder radius=%g, length=%g"%(radius, length)
    volume = pi*radius**2*length
    norm = CONTRAST**2*volume/10000*pi/length
    return long_cylinder

def make_sphere(radius):
    def sphere(qab, qc):
        q = sqrt(qab**2 + qc**2)
        return sas_3j1x_x(q*radius)
    sphere.__doc__ = "sphere radius=%g"%(radius,)
    volume = 4*pi*radius**3/3
    norm = CONTRAST**2*volume/10000
    return norm, sphere

THETA_LOW, THETA_HIGH = 0, pi/2
SCALE = 1


def kernel_1d(q, theta):
    """
    S(q) kernel for paracrystal forms.
    """
    qab = q*sin(theta)
    qc = q*cos(theta)
    return NORM*KERNEL(qab, qc)**2

def gauss_quad_1d(q, n=150):
    """
    Compute the integral using gaussian quadrature for n = 20, 76 or 150.
    """
    z, w = leggauss(n)
    theta = (THETA_HIGH-THETA_LOW)*(z + 1)/2 + THETA_LOW
    sin_theta = abs(sin(theta))
    Zq = kernel_1d(q=q, theta=theta)
    return np.sum(Zq*w*sin_theta)*(THETA_HIGH-THETA_LOW)/2

def gridded_1d(q, n=300):
    """
    Compute the integral on a regular grid using rectangular, trapezoidal,
    simpsons, and romberg integration.  Romberg integration requires that
    the grid be of size n = 2**k + 1.
    """
    theta = np.linspace(THETA_LOW, THETA_HIGH, n)
    Zq = kernel_1d(q=q, theta=theta)
    Zq *= abs(sin(theta))
    dx = theta[1]-theta[0]
    print("rect-%d"%n, np.sum(Zq)*dx*SCALE)
    print("trapz-%d"%n, np.trapz(Zq, dx=dx)*SCALE)
    print("simpson-%d"%n, simps(Zq, dx=dx)*SCALE)
    print("romb-%d"%n, romb(Zq, dx=dx)*SCALE)

def scipy_romberg_1d(q):
    """
    Compute the integral using romberg integration.  This function does not
    complete in a reasonable time.  No idea if it is accurate.
    """
    evals = [0]
    def outer(theta):
        evals[0] += 1
        return kernel_1d(q, theta=theta)*abs(sin(theta))
    result = romberg(outer, THETA_LOW, THETA_HIGH, divmax=100)*SCALE
    print("scipy romberg", evals[0], result)

def plot_1d(q, n=300):
    """
    Plot the function that needs to be integrated in order to compute
    the I(q) at a particular q.  *n* is the number of points in the grid.
    """
    theta = np.linspace(THETA_LOW, THETA_HIGH, n)
    Zq = kernel_1d(q=q, theta=theta)
    Zq *= abs(sin(theta))
    pylab.semilogy(degrees(theta), np.fmax(Zq, 1.e-6), label="Q=%g"%q)
    pylab.title("%s I(q, theta) sin(theta)" % (KERNEL.__doc__,))
    pylab.xlabel("theta (degrees)")
    pylab.ylabel("Iq 1/cm")

def Iq_trapz(q, n):
    theta = np.linspace(THETA_LOW, THETA_HIGH, n)
    Zq = kernel_1d(q=q, theta=theta)
    Zq *= abs(sin(theta))
    dx = theta[1]-theta[0]
    return np.trapz(Zq, dx=dx)*SCALE

def plot_Iq(q, n, form="trapz"):
    if form == "trapz":
        Iq = np.array([Iq_trapz(qk, n) for qk in q])
    elif form == "gauss":
        Iq = np.array([gauss_quad_1d(qk, n) for qk in q])
    pylab.loglog(q, Iq, label="%s, n=%d"%(form, n))
    pylab.xlabel("q (1/A)")
    pylab.ylabel("Iq (1/cm)")
    pylab.title(KERNEL.__doc__ + " I(q) circular average")
    return Iq

radius = 10.
length = 1e5
NORM, KERNEL = make_cylinder(radius=radius, length=length)
long_cyl = make_long_cylinder(radius=radius, length=length)
#NORM, KERNEL = make_sphere(radius=50.)


if __name__ == "__main__":
    Q = 0.386
    for n in (20, 76, 150, 300, 1000): #, 10000, 30000):
        print("gauss-%d"%n, gauss_quad_1d(Q, n=n))
    for k in (8, 10, 13, 16, 19):
        gridded_1d(Q, n=2**k+1)
    #print("inf cyl", 0, long_cyl(Q))
    #scipy_romberg(Q)

    plot_1d(0.386, n=2000)
    plot_1d(0.5, n=2000)
    plot_1d(0.8, n=2000)
    pylab.legend()
    pylab.figure()

    q = np.logspace(-3, 0, 400)
    I1 = long_cyl(q)
    I2 = plot_Iq(q, n=2**19+1, form="trapz")
    #plot_Iq(q, n=2**16+1, form="trapz")
    #plot_Iq(q, n=2**10+1, form="trapz")
    plot_Iq(q, n=1024, form="gauss")
    #plot_Iq(q, n=300, form="gauss")
    #plot_Iq(q, n=150, form="gauss")
    #plot_Iq(q, n=76, form="gauss")
    pylab.loglog(q, long_cyl(q), label="limit")
    pylab.legend()

    pylab.figure()
    pylab.semilogx(q, (I2 - I1)/I1)

    pylab.show()
