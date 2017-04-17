from __future__ import print_function, division

import numpy as np
from numpy import pi, sin, cos, exp, expm1, degrees, log10
from scipy.integrate import dblquad, simps, romb, romberg
import pylab

from sasmodels.special import square
from sasmodels.special import Gauss20Wt, Gauss20Z
from sasmodels.special import Gauss76Wt, Gauss76Z
from sasmodels.special import Gauss150Wt, Gauss150Z

def square(x):
    return x**2

Q = 0.001
DNN = 220.
D_FACTOR = 0.06
RADIUS = 40.0
SLD = 3.0
SLD_SOLVENT = 6.3

def kernel(q, dnn, d_factor, theta, phi):
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

    arg = 0.5*square(dnn*d_factor)*(a1**2 + a2**2 + a3**2)
    exp_arg = exp(-arg)
    den = [((exp_arg - 2*cos(dcos*a))*exp_arg + 1.0) for a in (a1, a2, a3)]
    Sq = -expm1(-2*arg)**3/np.prod(den, axis=0)
    return Sq


def scipy_dblquad(q=Q, dnn=DNN, d_factor=D_FACTOR):
    def integrand(theta, phi):
        Sq = kernel(q=q, dnn=dnn, d_factor=d_factor, theta=theta, phi=phi)
        return Sq*sin(theta)
    return dblquad(integrand, 0, pi/2, lambda x: 0, lambda x: pi/2)[0]*8/(4*pi)


def scipy_romberg(q=Q, dnn=DNN, d_factor=D_FACTOR):
    def inner(phi, theta):
        return kernel(q=q, dnn=dnn, d_factor=d_factor, theta=theta, phi=phi)
    def outer(theta):
        return romberg(inner, 0, pi/2, divmax=100, args=(theta,))*sin(theta)
    return romberg(outer, 0, pi/2, divmax=100)*8/(4*pi)


def gauss_quad(q=Q, dnn=DNN, d_factor=D_FACTOR, n=150):
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
    return np.sum(Sq*Aw*sin_theta)*8/(4*pi)


def gridded_integrals(q=0.001, dnn=DNN, d_factor=D_FACTOR, n=300):
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

def plot(q=0.001, dnn=DNN, d_factor=D_FACTOR, n=300):
    #theta = np.linspace(0, pi, n)
    #phi = np.linspace(0, 2*pi, n)
    theta = np.linspace(0, pi/2, n)
    phi = np.linspace(0, pi/2, n)
    Atheta, Aphi = np.meshgrid(theta, phi)
    Sq = kernel(q=Q, dnn=dnn, d_factor=d_factor, theta=Atheta, phi=Aphi)
    Sq *= abs(sin(Atheta))
    pylab.pcolor(degrees(theta), degrees(phi), log10(np.fmax(Sq, 1.e-6)))
    pylab.title("I(q) for q=%g"%q)
    pylab.xlabel("theta (degrees)")
    pylab.ylabel("phi (degrees)")
    cbar = pylab.colorbar()
    cbar.set_label('log10 I(q)')
    pylab.show()

if __name__ == "__main__":
    print("gauss", 20, gauss_quad(n=20))
    print("gauss", 76, gauss_quad(n=76))
    print("gauss", 150, gauss_quad(n=150))
    print("dblquad", scipy_dblquad())
    gridded_integrals(n=2**8+1)
    gridded_integrals(n=2**10+1)
    gridded_integrals(n=2**13+1)
    #print("romberg", scipy_romberg())
    #plot(n=300)