"""
Explore integration of rotationally symmetric shapes
"""

from __future__ import print_function, division

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np
from numpy import pi, sin, cos, sqrt, exp, expm1, degrees, log10
from scipy.integrate import dblquad, simps, romb, romberg
import pylab

from sasmodels.special import square
from sasmodels.special import Gauss20Wt, Gauss20Z
from sasmodels.special import Gauss76Wt, Gauss76Z
from sasmodels.special import Gauss150Wt, Gauss150Z
from sasmodels.special import sas_2J1x_x, sas_sinx_x, sas_3j1x_x

SLD = 3.0
SLD_SOLVENT = 6.3
CONTRAST = SLD - SLD_SOLVENT

def make_cylinder(radius, length):
    def cylinder(qab, qc):
        return sas_2J1x_x(qab*radius) * sas_sinx_x(qc*0.5*length)
    volume = pi*radius**2*length
    norm = 1e-4*volume*CONTRAST**2
    return norm, cylinder

def make_sphere(radius):
    def sphere(qab, qc):
        q = sqrt(qab**2 + qc**2)
        return sas_3j1x_x(q*radius)
    volume = 4*pi*radius**3/3
    norm = 1e-4*volume*CONTRAST**2
    return norm, sphere 
    
THETA_LOW, THETA_HIGH = 0, pi
SCALE = 1


def kernel(q, theta):
    """
    S(q) kernel for paracrystal forms.
    """
    qab = q*sin(theta)
    qc = q*cos(theta)
    return NORM*KERNEL(qab, qc)**2


def gauss_quad(q, n=150):
    """
    Compute the integral using gaussian quadrature for n = 20, 76 or 150.
    """
    if n == 20:
        z, w = Gauss20Z, Gauss20Wt
    elif n == 76:
        z, w = Gauss76Z, Gauss76Wt
    else:
        z, w = Gauss150Z, Gauss150Wt
    theta = (THETA_HIGH-THETA_LOW)*(z + 1)/2 + THETA_LOW
    sin_theta = abs(sin(theta))
    Zq = kernel(q=q, theta=theta)
    return np.sum(Zq*w*sin_theta)*SCALE/2


def gridded_integrals(q, n=300):
    """
    Compute the integral on a regular grid using rectangular, trapezoidal,
    simpsons, and romberg integration.  Romberg integration requires that
    the grid be of size n = 2**k + 1.
    """
    theta = np.linspace(THETA_LOW, THETA_HIGH, n)
    Zq = kernel(q=q, theta=theta)
    Zq *= abs(sin(theta))
    dx = theta[1]-theta[0]
    print("rect", n, np.sum(Zq)*dx*SCALE/pi)
    print("trapz", n, np.trapz(Zq, dx=dx)*SCALE/pi)
    print("simpson", n, simps(Zq, dx=dx)*SCALE/pi)
    print("romb", n, romb(Zq, dx=dx)*SCALE/pi)

def scipy_romberg(q):
    """
    Compute the integral using romberg integration.  This function does not
    complete in a reasonable time.  No idea if it is accurate.
    """
    evals = [0]
    def outer(theta):
        evals[0] += 1
        return kernel(q, theta=theta)*abs(sin(theta))
    result = romberg(outer, THETA_LOW, THETA_HIGH, divmax=100)*SCALE/pi
    print("scipy romberg", evals[0], result)

def plot(q, n=300):
    """
    Plot the 2D surface that needs to be integrated in order to compute
    the BCC S(q) at a particular q, dnn and d_factor.  *n* is the number
    of points in the grid.
    """
    theta = np.linspace(THETA_LOW, THETA_HIGH, n)
    Zq = kernel(q=q, theta=theta)
    Zq *= abs(sin(theta))
    pylab.semilogy(degrees(theta), np.fmax(Zq, 1.e-6), label="Q=%g"%q)
    pylab.title("%s I(q, theta) sin(theta)" % (KERNEL.__name__,))
    pylab.xlabel("theta (degrees)")
    pylab.ylabel("Iq 1/cm")

def Iq_trapz(q, n):
    theta = np.linspace(THETA_LOW, THETA_HIGH, n)
    Zq = kernel(q=q, theta=theta)
    Zq *= abs(sin(theta))
    dx = theta[1]-theta[0]
    return np.trapz(Zq, dx=dx)*SCALE/pi

def plot_Iq(q, n, form="trapz"):
    if form == "trapz":
        I = np.array([Iq_trapz(qk, n) for qk in q])
    elif form == "gauss":
        I = np.array([gauss_quad(qk, n) for qk in q])
    pylab.loglog(q, I, label="%s, n=%d"%(form, n))
    pylab.xlabel("q (1/A)")
    pylab.ylabel("Iq (1/cm)")

NORM, KERNEL = make_cylinder(radius=10., length=100000.)
#NORM, KERNEL = make_cylinder(radius=10., length=10000.)
#NORM, KERNEL = make_cylinder(radius=10., length=30.)
#NORM, KERNEL = make_sphere(radius=50.)

if __name__ == "__main__":
    Q = 0.8
    print("gauss", 20, gauss_quad(Q, n=20))
    print("gauss", 76, gauss_quad(Q, n=76))
    print("gauss", 150, gauss_quad(Q, n=150))
    gridded_integrals(Q, n=2**8+1)
    gridded_integrals(Q, n=2**10+1)
    gridded_integrals(Q, n=2**13+1)
    gridded_integrals(Q, n=2**16+1)
    gridded_integrals(Q, n=2**19+1)
    #scipy_romberg(Q)
    plot(0.5, n=2000)
    plot(0.6, n=2000)
    plot(0.8, n=2000)
    pylab.legend()
    pylab.figure()
    plot_Iq(np.logspace(-3,0,200), n=2**16+1, form="trapz")
    plot_Iq(np.logspace(-3,0,200), n=2**10+1, form="trapz")
    plot_Iq(np.logspace(-3,0,200), n=150, form="gauss")
    pylab.legend()
    pylab.show()
