"""
The current 1D calculations for BCC paracrystal are very wrong at low q, orders
of magnitude wrong.  The integration fails to capture a very narrow,
very steep ridge.

Uncomment the plot() line at the bottom of the code to show an image of the
set of S(q, theta, phi) values that must be summed to form the 1D S(q=0.001)
for the BCC paracrystal form.  Note particularly that this is a log scale
image spanning 10 orders of magnitude.  This pattern repeats itself 8 times
over the entire 4 pi surface integral.

You can explore various integration options by uncommenting more lines.
Adaptive integration using scpy.integrate.dbsquad is very slow.  Romberg
didn't even complete in the time I gave it. Accurate brute force calculation
requires a 4000x4000 grid to get enough precision.

We may need a specialized integrator for low q which can identify and integrate
the ridges properly.

This program need sasmodels on the path so it is inserted automatically,
assuming that the explore/bccpy.py is beside sasmodels/special.py in the
source tree.  Run from the sasmodels directory using:

    python explore/bccpy.py
"""

from __future__ import print_function, division

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np
from numpy import pi, sin, cos, exp, expm1, degrees, log10
from scipy.integrate import dblquad, simps, romb, romberg
import pylab

from sasmodels.special import square
from sasmodels.special import Gauss20Wt, Gauss20Z
from sasmodels.special import Gauss76Wt, Gauss76Z
from sasmodels.special import Gauss150Wt, Gauss150Z

Q = 0.001
DNN = 220.
D_FACTOR = 0.06
RADIUS = 40.0
SLD = 3.0
SLD_SOLVENT = 6.3

def kernel(q, dnn, d_factor, theta, phi):
    """
    S(q) kernel for paracrystal forms.
    """
    qab = q*sin(theta)
    qa = qab*cos(phi)
    qb = qab*sin(phi)
    qc = q*cos(theta)

    if 0: # sc
        a1, a2, a3 = qa, qb, qc
        dcos = dnn
    if 1: # bcc
        a1 = +qa - qc + qb
        a2 = +qa + qc - qb
        a3 = -qa + qc + qb
        dcos = dnn/2
    if 0: # fcc
        a1 = qb + qa
        a2 = qa + qc
        a3 = qb + qc
        dcos = dnn/2

    arg = 0.5*square(dnn*d_factor)*(a1**2 + a2**2 + a3**2)
    exp_arg = exp(-arg)
    den = [((exp_arg - 2*cos(dcos*a))*exp_arg + 1.0) for a in (a1, a2, a3)]
    Sq = -expm1(-2*arg)**3/np.prod(den, axis=0)
    return Sq


def scipy_dblquad(q=Q, dnn=DNN, d_factor=D_FACTOR):
    """
    Compute the integral using scipy dblquad.  This gets the correct answer
    eventually, but it is slow.
    """
    evals = [0]
    def integrand(theta, phi):
        evals[0] += 1
        Sq = kernel(q=q, dnn=dnn, d_factor=d_factor, theta=theta, phi=phi)
        return Sq*sin(theta)
    ans = dblquad(integrand, 0, pi/2, lambda x: 0, lambda x: pi/2)[0]*8/(4*pi)
    print("dblquad evals =", evals[0])
    return ans


def scipy_romberg_2d(q=Q, dnn=DNN, d_factor=D_FACTOR):
    """
    Compute the integral using romberg integration.  This function does not
    complete in a reasonable time.  No idea if it is accurate.
    """
    def inner(phi, theta):
        return kernel(q=q, dnn=dnn, d_factor=d_factor, theta=theta, phi=phi)
    def outer(theta):
        return romberg(inner, 0, pi/2, divmax=100, args=(theta,))*sin(theta)
    return romberg(outer, 0, pi/2, divmax=100)*8/(4*pi)


def semi_romberg(q=Q, dnn=DNN, d_factor=D_FACTOR, n=100):
    """
    Use 1D romberg integration in phi and regular simpsons rule in theta.
    """
    evals = [0]
    def inner(phi, theta):
        evals[0] += 1
        return kernel(q=q, dnn=dnn, d_factor=d_factor, theta=theta, phi=phi)
    theta = np.linspace(0, pi/2, n)
    f_phi = [romberg(inner, 0, pi/2, divmax=100, args=(t,))
             for t in theta]
    ans = simps(sin(theta)*np.array(f_phi), dx=theta[1]-theta[0])
    print("semi romberg evals =", evals[0])
    return ans*8/(4*pi)

def gauss_quad(q=Q, dnn=DNN, d_factor=D_FACTOR, n=150):
    """
    Compute the integral using gaussian quadrature for n = 20, 76 or 150.
    """
    if n == 20:
        z, w = Gauss20Z, Gauss20Wt
    elif n == 76:
        z, w = Gauss76Z, Gauss76Wt
    else:
        z, w = Gauss150Z, Gauss150Wt
    theta = pi/4*(z + 1)
    phi = pi/4*(z + 1)
    Atheta, Aphi = np.meshgrid(theta, phi)
    Aw = w[None, :] * w[:, None]
    sin_theta = np.fmax(abs(sin(Atheta)), 1e-6)
    Sq = kernel(q=q, dnn=dnn, d_factor=d_factor, theta=Atheta, phi=Aphi)
    print("gauss %d evals ="%n, n**2)
    return np.sum(Sq*Aw*sin_theta)*8/(4*pi)


def gridded_integrals(q=0.001, dnn=DNN, d_factor=D_FACTOR, n=300):
    """
    Compute the integral on a regular grid using rectangular, trapezoidal,
    simpsons, and romberg integration.  Romberg integration requires that
    the grid be of size n = 2**k + 1.
    """
    theta = np.linspace(0, pi/2, n)
    phi = np.linspace(0, pi/2, n)
    Atheta, Aphi = np.meshgrid(theta, phi)
    Sq = kernel(q=Q, dnn=dnn, d_factor=d_factor, theta=Atheta, phi=Aphi)
    Sq *= abs(sin(Atheta))
    dx, dy = theta[1]-theta[0], phi[1]-phi[0]
    print("rect", n, np.sum(Sq)*dx*dy*8/(4*pi))
    print("trapz", n, np.trapz(np.trapz(Sq, dx=dx), dx=dy)*8/(4*pi))
    print("simpson", n, simps(simps(Sq, dx=dx), dx=dy)*8/(4*pi))
    print("romb", n, romb(romb(Sq, dx=dx), dx=dy)*8/(4*pi))
    print("gridded %d evals ="%n, n**2)

def plot(q=0.001, dnn=DNN, d_factor=D_FACTOR, n=300):
    """
    Plot the 2D surface that needs to be integrated in order to compute
    the BCC S(q) at a particular q, dnn and d_factor.  *n* is the number
    of points in the grid.
    """
    theta = np.linspace(0, pi, n)
    phi = np.linspace(0, 2*pi, n)
    #theta = np.linspace(0, pi/2, n)
    #phi = np.linspace(0, pi/2, n)
    Atheta, Aphi = np.meshgrid(theta, phi)
    Sq = kernel(q=Q, dnn=dnn, d_factor=d_factor, theta=Atheta, phi=Aphi)
    Sq *= abs(sin(Atheta))
    pylab.pcolor(degrees(theta), degrees(phi), log10(np.fmax(Sq, 1.e-6)))
    pylab.axis('tight')
    pylab.title("BCC S(q) for q=%g, dnn=%g d_factor=%g" % (q, dnn, d_factor))
    pylab.xlabel("theta (degrees)")
    pylab.ylabel("phi (degrees)")
    cbar = pylab.colorbar()
    cbar.set_label('log10 S(q)')
    pylab.show()

if __name__ == "__main__":
    #print("gauss", 20, gauss_quad(n=20))
    #print("gauss", 76, gauss_quad(n=76))
    #print("gauss", 150, gauss_quad(n=150))
    #print("dblquad", scipy_dblquad())
    #print("semi romberg", semi_romberg())
    #gridded_integrals(n=2**8+1)
    #gridded_integrals(n=2**10+1)
    #gridded_integrals(n=2**13+1)
    #print("romberg", scipy_romberg())
    plot(n=400)
