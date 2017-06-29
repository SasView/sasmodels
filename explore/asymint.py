"""
Asymmetric shape integration
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
    def cylinder(qa, qb, qc):
        qab = sqrt(qa**2 + qb**2)
        return sas_2J1x_x(qab*radius) * sas_sinx_x(qc*0.5*length)
    volume = pi*radius**2*length
    norm = 1e-4*volume*CONTRAST**2
    return norm, cylinder

def make_sphere(radius):
    def sphere(qa, qb, qc):
        qab = sqrt(qa**2 + qb**2)
        q = sqrt(qab**2 + qc**2)
        return sas_3j1x_x(q*radius)
    volume = 4*pi*radius**3/3
    norm = 1e-4*volume*CONTRAST**2
    return norm, sphere 
    
#NORM, KERNEL = make_cylinder(radius=10., length=100000.)
NORM, KERNEL = make_cylinder(radius=10., length=30.)
#NORM, KERNEL = make_sphere(radius=50.)

THETA_LOW, THETA_HIGH = 0, pi
PHI_LOW, PHI_HIGH = 0, 2*pi
SCALE = 1


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
    def integrand(theta, phi):
        evals[0] += 1
        Zq = kernel(q, theta=theta, phi=phi)
        return Zq*sin(theta)
    ans = dblquad(integrand, THETA_LOW, THETA_HIGH, lambda x: PHI_LOW, lambda x: PHI_HIGH)[0]*SCALE/(4*pi)
    print("dblquad evals =", evals[0])
    return ans


def scipy_romberg_2d(q):
    """
    Compute the integral using romberg integration.  This function does not
    complete in a reasonable time.  No idea if it is accurate.
    """
    def inner(phi, theta):
        return kernel(q, theta=theta, phi=phi)
    def outer(theta):
        return romberg(inner, PHI_LOW, PHI_HIGH, divmax=100, args=(theta,))*sin(theta)
    return romberg(outer, THETA_LOW, THETA_HIGH, divmax=100)*SCALE/(4*pi)


def semi_romberg(q, n=100):
    """
    Use 1D romberg integration in phi and regular simpsons rule in theta.
    """
    evals = [0]
    def inner(phi, theta):
        evals[0] += 1
        return kernel(q, theta=theta, phi=phi)
    theta = np.linspace(THETA_LOW, THETA_HIGH, n)
    f_phi = [romberg(inner, PHI_LOW, PHI_HIGH, divmax=100, args=(t,))
             for t in theta]
    ans = simps(sin(theta)*np.array(f_phi), dx=theta[1]-theta[0])
    print("semi romberg evals =", evals[0])
    return ans*SCALE/(4*pi)

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
    phi = (PHI_HIGH-PHI_LOW)*(z + 1)/2 + PHI_LOW
    Atheta, Aphi = np.meshgrid(theta, phi)
    Aw = w[None, :] * w[:, None]
    sin_theta = np.fmax(abs(sin(Atheta)), 1e-6)
    Zq = kernel(q=q, theta=Atheta, phi=Aphi)
    print("gauss %d evals ="%n, n**2)
    return np.sum(Zq*Aw*sin_theta)*SCALE/(4*pi)


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
    print("rect", n, np.sum(Zq)*dx*dy*SCALE/(4*pi))
    print("trapz", n, np.trapz(np.trapz(Zq, dx=dx), dx=dy)*SCALE/(4*pi))
    print("simpson", n, simps(simps(Zq, dx=dx), dx=dy)*SCALE/(4*pi))
    print("romb", n, romb(romb(Zq, dx=dx), dx=dy)*SCALE/(4*pi))
    print("gridded %d evals ="%n, n**2)

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
    Q = 0.8
    print("gauss", 20, gauss_quad(Q, n=20))
    print("gauss", 76, gauss_quad(Q, n=76))
    print("gauss", 150, gauss_quad(Q, n=150))
    #print("dblquad", scipy_dblquad(Q))
    #print("semi romberg", semi_romberg(Q))
    #gridded_integrals(Q, n=2**8+1)
    #gridded_integrals(Q, n=2**10+1)
    #gridded_integrals(Q, n=2**13+1)
    #print("romberg", scipy_romberg(Q))
    plot(Q, n=400)
