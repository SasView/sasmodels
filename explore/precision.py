#!/usr/bin/env python
r"""
Show numerical precision of various expressions.

Evaluates the same function(s) in single and double precision and compares
the results to 500 digit mpmath evaluation of the same function.

Note: a quick way to generation C and python code for taylor series
expansions from sympy:

    import sympy as sp
    x = sp.var("x")
    f = sp.sin(x)/x
    t = sp.series(f, n=12).removeO()  # taylor series with no O(x^n) term
    p = sp.horner(t)   # Horner representation
    p = p.replace(x**2, sp.var("xsq")  # simplify if alternate terms are zero
    p = p.n(15)  # evaluate coefficients to 15 digits (optional)
    c_code = sp.ccode(p, assign_to=sp.var("p"))  # convert to c code
    py_code = c[:-1]  # strip semicolon to convert c to python

    # mpmath has pade() rational function approximation, which might work
    # better than the taylor series for some functions:
    P, Q = mp.pade(sp.Poly(t.n(15),x).coeffs(), L, M)
    P = sum(a*x**n for n,a in enumerate(reversed(P)))
    Q = sum(a*x**n for n,a in enumerate(reversed(Q)))
    c_code = sp.ccode(sp.horner(P)/sp.horner(Q), assign_to=sp.var("p"))

    # There are richardson and shanks series accelerators in both sympy
    # and mpmath that may be helpful.
"""
from __future__ import division, print_function

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
# Note: mpmath.pi and numpy.pi are not interchangeable; don't import pi from
# numpy as we usually do, but instead be explicit about which one we want each
# time it we use it.
from numpy import inf
import scipy.special
try:
    from mpmath import mp
except ImportError:
    # CRUFT: mpmath used to be a package in sympy
    from sympy.mpmath import mp
#import matplotlib; matplotlib.use('TkAgg')
import pylab

from sasmodels import core, data, direct_model, modelinfo

class Comparator(object):
    def __init__(self, name, mp_function, np_function, ocl_function, xaxis, limits):
        self.name = name
        self.mp_function = mp_function
        self.np_function = np_function
        self.ocl_function = ocl_function
        self.xaxis = xaxis
        self.limits = limits

    def __repr__(self):
        return "Comparator(%s)"%self.name

    def call_mpmath(self, vec, bits=500):
        """
        Direct calculation using mpmath extended precision library.
        """
        with mp.workprec(bits):
            return [self.mp_function(mp.mpf(x)) for x in vec]

    def call_numpy(self, x, dtype):
        """
        Direct calculation using numpy/scipy.
        """
        x = np.asarray(x, dtype)
        return self.np_function(x)

    def call_ocl(self, x, dtype, platform='ocl'):
        """
        Calculation using sasmodels ocl libraries.
        """
        x = np.asarray(x, dtype)
        model = core.build_model(self.ocl_function, dtype=dtype)
        calculator = direct_model.DirectModel(data.empty_data1D(x), model)
        return calculator(background=0)

    def run(self, xrange="log", diff="relative"):
        r"""
        Compare accuracy of different methods for computing f.

        *xrange* is::

            log:    [10^-3,10^5]
            logq:   [10^-4, 10^1]
            linear: [1,1000]
            zoom:   [1000,1010]
            neg:    [-100,100]

        For arbitrary range use "start:stop:steps:scale" where scale is
        one of log, lin, or linear.

        *diff* is "relative", "absolute" or "none"

        *x_bits* is the precision with which the x values are specified.  The
        default 23 should reproduce the equivalent of a single precisio
        """
        linear = not xrange.startswith("log")
        if xrange == "zoom":
            start, stop, steps = 1000, 1010, 2000
        elif xrange == "neg":
            start, stop, steps = -100.1, 100.1, 2000
        elif xrange == "linear":
            start, stop, steps = 1, 1000, 2000
            start, stop, steps = 0.001, 2, 2000
        elif xrange == "log":
            start, stop, steps = -3, 5, 400
        elif xrange == "logq":
            start, stop, steps = -4, 1, 400
        elif ':' in xrange:
            parts = xrange.split(':')
            linear = parts[3] != "log" if len(parts) == 4 else True
            steps = int(parts[2]) if len(parts) > 2 else 400
            start = float(parts[0])
            stop = float(parts[1])

        else:
            raise ValueError("unknown range "+xrange)
        with mp.workprec(500):
            # Note: we make sure that we are comparing apples to apples...
            # The x points are set using single precision so that we are
            # examining the accuracy of the transformation from x to f(x)
            # rather than x to f(nearest(x)) where nearest(x) is the nearest
            # value to x in the given precision.
            if linear:
                start = max(start, self.limits[0])
                stop = min(stop, self.limits[1])
                qrf = np.linspace(start, stop, steps, dtype='single')
                #qrf = np.linspace(start, stop, steps, dtype='double')
                qr = [mp.mpf(float(v)) for v in qrf]
                #qr = mp.linspace(start, stop, steps)
            else:
                start = np.log10(max(10**start, self.limits[0]))
                stop = np.log10(min(10**stop, self.limits[1]))
                qrf = np.logspace(start, stop, steps, dtype='single')
                #qrf = np.logspace(start, stop, steps, dtype='double')
                qr = [mp.mpf(float(v)) for v in qrf]
                #qr = [10**v for v in mp.linspace(start, stop, steps)]

        target = self.call_mpmath(qr, bits=500)
        pylab.subplot(121)
        self.compare(qr, 'single', target, linear, diff)
        pylab.legend(loc='best')
        pylab.subplot(122)
        self.compare(qr, 'double', target, linear, diff)
        pylab.legend(loc='best')
        pylab.suptitle(self.name + " compared to 500-bit mpmath")

    def compare(self, x, precision, target, linear=False, diff="relative"):
        r"""
        Compare the different computation methods using the given precision.
        """
        if precision == 'single':
            #n=11; plotdiff(x, target, self.call_mpmath(x, n), 'mp %d bits'%n, diff=diff)
            #n=23; plotdiff(x, target, self.call_mpmath(x, n), 'mp %d bits'%n, diff=diff)
            pass
        elif precision == 'double':
            #n=53; plotdiff(x, target, self.call_mpmath(x, n), 'mp %d bits'%n, diff=diff)
            #n=83; plotdiff(x, target, self.call_mpmath(x, n), 'mp %d bits'%n, diff=diff)
            pass
        plotdiff(x, target, self.call_numpy(x, precision), 'numpy '+precision, diff=diff)
        plotdiff(x, target, self.call_ocl(x, precision, 0), 'OpenCL '+precision, diff=diff)
        pylab.xlabel(self.xaxis)
        if diff == "relative":
            pylab.ylabel("relative error")
        elif diff == "absolute":
            pylab.ylabel("absolute error")
        else:
            pylab.ylabel(self.name)
            pylab.semilogx(x, target, '-', label="true value")
        if linear:
            pylab.xscale('linear')

def plotdiff(x, target, actual, label, diff):
    """
    Plot the computed value.

    Use relative error if SHOW_DIFF, otherwise just plot the value directly.
    """
    if diff == "relative":
        err = np.array([(abs((t-a)/t) if t != 0 else a) for t, a in zip(target, actual)], 'd')
        #err = np.clip(err, 0, 1)
        pylab.loglog(x, err, '-', label=label, alpha=0.7)
    elif diff == "absolute":
        err = np.array([abs((t-a)) for t, a in zip(target, actual)], 'd')
        pylab.loglog(x, err, '-', label=label, alpha=0.7)
    else:
        limits = np.min(target), np.max(target)
        pylab.semilogx(x, np.clip(actual, *limits), '-', label=label, alpha=0.7)

def make_ocl(function, name, source=[]):
    class Kernel(object):
        pass
    Kernel.__file__ = name+".py"
    Kernel.name = name
    Kernel.parameters = []
    Kernel.source = source
    Kernel.Iq = function
    model_info = modelinfo.make_model_info(Kernel)
    return model_info

# Hack to allow second parameter A in the gammainc and gammaincc functions.
# Create a 2-D variant of the precision test if we need to handle other two
# parameter functions.
A = 1
def parse_extra_pars():
    """
    Parse the command line looking for the second parameter "A=..." for the
    gammainc/gammaincc functions.
    """
    global A

    A_str = str(A)
    pop = []
    for k, v in enumerate(sys.argv[1:]):
        if v.startswith("A="):
            A_str = v[2:]
            pop.append(k+1)
    if pop:
        sys.argv = [v for k, v in enumerate(sys.argv) if k not in pop]
        A = float(A_str)

parse_extra_pars()


# =============== FUNCTION DEFINITIONS ================

FUNCTIONS = {}
def add_function(name, mp_function, np_function, ocl_function,
                 shortname=None, xaxis="x", limits=(-inf, inf)):
    if shortname is None:
        shortname = name.replace('(x)', '').replace(' ', '')
    FUNCTIONS[shortname] = Comparator(name, mp_function, np_function, ocl_function, xaxis, limits)

add_function(
    name="J0(x)",
    mp_function=mp.j0,
    np_function=scipy.special.j0,
    ocl_function=make_ocl("return sas_J0(q);", "sas_J0", ["lib/polevl.c", "lib/sas_J0.c"]),
)
add_function(
    name="J1(x)",
    mp_function=mp.j1,
    np_function=scipy.special.j1,
    ocl_function=make_ocl("return sas_J1(q);", "sas_J1", ["lib/polevl.c", "lib/sas_J1.c"]),
)
add_function(
    name="JN(-3, x)",
    mp_function=lambda x: mp.besselj(-3, x),
    np_function=lambda x: scipy.special.jn(-3, x),
    ocl_function=make_ocl("return sas_JN(-3, q);", "sas_JN",
                          ["lib/polevl.c", "lib/sas_J0.c", "lib/sas_J1.c", "lib/sas_JN.c"]),
    shortname="J-3",
)
add_function(
    name="JN(3, x)",
    mp_function=lambda x: mp.besselj(3, x),
    np_function=lambda x: scipy.special.jn(3, x),
    ocl_function=make_ocl("return sas_JN(3, q);", "sas_JN",
                          ["lib/polevl.c", "lib/sas_J0.c", "lib/sas_J1.c", "lib/sas_JN.c"]),
    shortname="J3",
)
add_function(
    name="JN(2, x)",
    mp_function=lambda x: mp.besselj(2, x),
    np_function=lambda x: scipy.special.jn(2, x),
    ocl_function=make_ocl("return sas_JN(2, q);", "sas_JN",
                          ["lib/polevl.c", "lib/sas_J0.c", "lib/sas_J1.c", "lib/sas_JN.c"]),
    shortname="J2",
)
add_function(
    name="2 J1(x)/x",
    mp_function=lambda x: 2*mp.j1(x)/x,
    np_function=lambda x: 2*scipy.special.j1(x)/x,
    ocl_function=make_ocl("return sas_2J1x_x(q);", "sas_2J1x_x", ["lib/polevl.c", "lib/sas_J1.c"]),
)
add_function(
    name="J1(x)",
    mp_function=mp.j1,
    np_function=scipy.special.j1,
    ocl_function=make_ocl("return sas_J1(q);", "sas_J1", ["lib/polevl.c", "lib/sas_J1.c"]),
)
add_function(
    name="Si(x)",
    mp_function=mp.si,
    np_function=lambda x: scipy.special.sici(x)[0],
    ocl_function=make_ocl("return sas_Si(q);", "sas_Si", ["lib/sas_Si.c"]),
)
#import fnlib
#add_function(
#    name="fnlibJ1",
#    mp_function=mp.j1,
#    np_function=fnlib.J1,
#    ocl_function=make_ocl("return sas_J1(q);", "sas_J1", ["lib/polevl.c", "lib/sas_J1.c"]),
#)
add_function(
    name="sin(x)",
    mp_function=mp.sin,
    np_function=np.sin,
    #ocl_function=make_ocl("double sn, cn; SINCOS(q,sn,cn); return sn;", "sas_sin"),
    ocl_function=make_ocl("return sin(q);", "sas_sin"),
)
add_function(
    name="sin(x)/x",
    mp_function=lambda x: mp.sin(x)/x if x != 0 else 1,
    ## scipy sinc function is inaccurate and has an implied pi*x term
    #np_function=lambda x: scipy.special.sinc(x/pi),
    ## numpy sin(x)/x needs to check for x=0
    np_function=lambda x: np.sin(x)/x,
    ocl_function=make_ocl("return sas_sinx_x(q);", "sas_sinc"),
)
add_function(
    name="cos(x)",
    mp_function=mp.cos,
    np_function=np.cos,
    #ocl_function=make_ocl("double sn, cn; SINCOS(q,sn,cn); return cn;", "sas_cos"),
    ocl_function=make_ocl("return cos(q);", "sas_cos"),
)
add_function(
    name="gamma(x)",
    mp_function=mp.gamma,
    np_function=scipy.special.gamma,
    ocl_function=make_ocl("return sas_gamma(q);", "sas_gamma", ["lib/sas_gamma.c"]),
    limits=(-3.1, 10),
)
add_function(
    name="gammaln(x)",
    mp_function=mp.loggamma,
    np_function=scipy.special.gammaln,
    ocl_function=make_ocl("return sas_gammaln(q);", "sas_gammaln", ["lib/sas_gammainc.c"]),
    #ocl_function=make_ocl("return lgamma(q);", "sas_gammaln"),
)
add_function(
    # Note: "a" is given as A=... on the command line via parse_extra_pars
    name="gammainc(x)",
    mp_function=lambda x, a=A: mp.gammainc(a, a=0, b=x)/mp.gamma(a),
    np_function=lambda x, a=A: scipy.special.gammainc(a, x),
    ocl_function=make_ocl("return sas_gammainc(%.15g,q);"%A, "sas_gammainc", ["lib/sas_gammainc.c"]),
)
add_function(
    # Note: "a" is given as A=... on the command line via parse_extra_pars
    name="gammaincc(x)",
    mp_function=lambda x, a=A: mp.gammainc(a, a=x, b=mp.inf)/mp.gamma(a),
    np_function=lambda x, a=A: scipy.special.gammaincc(a, x),
    ocl_function=make_ocl("return sas_gammaincc(%.15g,q);"%A, "sas_gammaincc", ["lib/sas_gammainc.c"]),
)
add_function(
    name="erf(x)",
    mp_function=mp.erf,
    np_function=scipy.special.erf,
    ocl_function=make_ocl("return sas_erf(q);", "sas_erf", ["lib/polevl.c", "lib/sas_erf.c"]),
    limits=(-5., 5.),
)
add_function(
    name="erfc(x)",
    mp_function=mp.erfc,
    np_function=scipy.special.erfc,
    ocl_function=make_ocl("return sas_erfc(q);", "sas_erfc", ["lib/polevl.c", "lib/sas_erf.c"]),
    limits=(-5., 5.),
)
add_function(
    name="expm1(x)",
    mp_function=mp.expm1,
    np_function=np.expm1,
    ocl_function=make_ocl("return expm1(q);", "sas_expm1"),
    limits=(-5., 5.),
)
add_function(
    name="arctan(x)",
    mp_function=mp.atan,
    np_function=np.arctan,
    ocl_function=make_ocl("return atan(q);", "sas_arctan"),
)
add_function(
    name="3 j1(x)/x",
    mp_function=lambda x: 3*(mp.sin(x)/x - mp.cos(x))/(x*x),
    # Note: no taylor expansion near 0
    np_function=lambda x: 3*(np.sin(x)/x - np.cos(x))/(x*x),
    ocl_function=make_ocl("return sas_3j1x_x(q);", "sas_j1c", ["lib/sas_3j1x_x.c"]),
)
add_function(
    name="(1-cos(x))/x^2",
    mp_function=lambda x: (1 - mp.cos(x))/(x*x),
    np_function=lambda x: (1 - np.cos(x))/(x*x),
    ocl_function=make_ocl("return (1-cos(q))/q/q;", "sas_1mcosx_x2"),
)
add_function(
    name="(1-sin(x)/x)/x",
    mp_function=lambda x: 1/x - mp.sin(x)/(x*x),
    np_function=lambda x: 1/x - np.sin(x)/(x*x),
    ocl_function=make_ocl("return (1-sas_sinx_x(q))/q;", "sas_1msinx_x_x"),
)
add_function(
    name="(1/2-sin(x)/x+(1-cos(x))/x^2)/x",
    mp_function=lambda x: (0.5 - mp.sin(x)/x + (1-mp.cos(x))/(x*x))/x,
    np_function=lambda x: (0.5 - np.sin(x)/x + (1-np.cos(x))/(x*x))/x,
    ocl_function=make_ocl("return (0.5-sin(q)/q + (1-cos(q))/q/q)/q;", "sas_T2"),
)
add_function(
    name="fmod_2pi",
    mp_function=lambda x: mp.fmod(x, 2*mp.pi),
    np_function=lambda x: np.fmod(x, 2*np.pi),
    ocl_function=make_ocl("return fmod(q, 2*M_PI);", "sas_fmod"),
)
add_function(
    name="expm1(x)/x",
    # Note: should be 1 when x = 0
    mp_function=lambda x: mp.expm1(x)/x,
    np_function=lambda x: np.expm1(x)/x,
    ocl_function=make_ocl("return (q==0.) ? 1. : expm1(q)/q;", "sas_exp1_x"),
)
add_function(
    name="sq_expm1(x)/x",
    # Note: should be 1 when x = 0
    mp_function=lambda x: (mp.expm1(x)/x)**2,
    np_function=lambda x: (np.expm1(x)/x)**2,
    ocl_function=make_ocl("return (q==0.) ? 1. : square(expm1(q)/q);", "sas_square_exp1_x"),
)

# TODO: move to sas_special
def sas_langevin(x):
    scalar = np.isscalar(x)
    if scalar:
        x = np.array([x]) # should inherit dtype for single if given single
    f = np.empty_like(x)
    cutoff = 0.1 if f.dtype == np.float64 else 1.0
    #cutoff *= 10
    index = x < cutoff
    xp = x[index]
    xpsq = xp*xp
    f[index] = xp / (3. + xpsq / (5. + xpsq/(7. + xpsq/(9.))))
    # 4 terms gets to 1e-7 single, 1e-14 double. Can get to 1e-15 double by adding
    # another 4 terms and setting cutoff at 1.0. Not worthwhile. Instead we would
    # need an expansion about x somewhere between 1 and 10 for the interval [0.1, 100.]
    #f[index] = xp / (3. + xpsq / (5. + xpsq/(7. + xpsq/(9. + xpsq/(11.0 + xpsq/(13. + xpsq/(15. + xpsq/17.)))))))
    xp = x[~index]
    f[~index] = 1/np.tanh(xp) - 1/xp
    return f[0] if scalar else f

def sas_langevin_x(x):
    scalar = np.isscalar(x)
    if scalar:
        x = np.array([x]) # should inherit dtype for single if given single
    f = np.empty_like(x)
    cutoff = 0.1 if f.dtype == np.float64 else 1.0
    index = x < cutoff
    xp = x[index]
    xpsq = xp*xp
    f[index] = 1. / (3. + xpsq / (5. + xpsq/(7. + xpsq/(9.))))
    xp = x[~index]
    f[~index] = (1/np.tanh(xp) - 1/xp)/xp
    return f[0] if scalar else f

add_function(
    name="langevin(x)",
    mp_function=lambda x: (1/mp.tanh(x) - 1/x),
    np_function=sas_langevin,
    #ocl_function=make_ocl("return q < 0.7 ? q*(1./3. + q*q*(-1./45. + q*q*(2./945. + q*q*(-1./4725.) + q*q*(2./93555.)))) : 1/tanh(q) - 1/q;", "sas_langevin"),
    #ocl_function=make_ocl("return q < 1e-5 ? q/3. : 1/tanh(q) - 1/q;", "sas_langevin"),
    ocl_function=make_ocl("""
#if FLOAT_SIZE>4  // DOUBLE_PRECISION
#  define LANGEVIN_CUTOFF 0.1
#else
#  define LANGEVIN_CUTOFF 1.0
#endif
    const double qsq = q*q;
    return (q < LANGEVIN_CUTOFF) ? q / (3. + qsq / (5. + qsq/(7. + qsq/(9.)))) : 1/tanh(q) - 1/q;
    """, "sas_langevin"),
)
add_function(
    name="langevin(x)/x",
    mp_function=lambda x: (1/mp.tanh(x) - 1/x)/x,
    np_function=sas_langevin_x,
    ocl_function=make_ocl("""
#if FLOAT_SIZE>4  // DOUBLE_PRECISION
#  define LANGEVIN_CUTOFF 0.1
#else
#  define LANGEVIN_CUTOFF 1.0
#endif
    const double qsq = q*q;
    return (q < LANGEVIN_CUTOFF) ? 1. / (3. + qsq / (5. + qsq/(7. + qsq/(9.)))) : (1/tanh(q) - 1/q)/q;
    """, "sas_langevin_x"),
)
add_function(
    name="gauss_coil",
    mp_function=lambda x: 2*(mp.exp(-x**2) + x**2 - 1)/x**4,
    np_function=lambda x: 2*(np.expm1(-x**2) + x**2)/x**4,
    ocl_function=make_ocl("""
    const double qsq = q*q;
    // For double: use O(5) Pade with 0.5 cutoff (10 mad + 1 divide)
    // For single: use O(7) Taylor with 0.8 cutoff (7 mad)
    if (qsq < 0.0) {
        const double x = qsq;
        if (0) { // 0.36 single
            // PadeApproximant[2*Exp[-x^2] + x^2-1)/x^4, {x, 0, 4}]
            return (x*x/180. + 1.)/((1./30.*x + 1./3.)*x + 1);
        } else if (0) { // 1.0 for single
            // padeapproximant[2*exp[-x^2] + x^2-1)/x^4, {x, 0, 6}]
            const double A1=1./24., A2=1./84, A3=-1./3360;
            const double B1=3./8., B2=3./56., B3=1./336.;
            return (((A3*x + A2)*x + A1)*x + 1.)/(((B3*x + B2)*x + B1)*x + 1.);
        } else if (0) { // 1.0 for single, 0.25 for double
            // PadeApproximant[2*Exp[-x^2] + x^2-1)/x^4, {x, 0, 8}]
            const double A1=1./15., A2=1./60, A3=0., A4=1./75600.;
            const double B1=2./5., B2=1./15., B3=1./180., B4=1./5040.;
            return ((((A4*x + A3)*x + A2)*x + A1)*x + 1.)
                  /((((B4*x + B3)*x + B2)*x + B1)*x + 1.);
        } else { // 1.0 for single, 0.5 for double
            // PadeApproximant[2*Exp[-x^2] + x^2-1)/x^4, {x, 0, 8}]
            const double A1=1./12., A2=2./99., A3=1./2640., A4=1./23760., A5=-1./1995840.;
            const double B1=5./12., B2=5./66., B3=1./132., B4=1./2376., B5=1./95040.;
            return (((((A5*x + A4)*x + A3)*x + A2)*x + A1)*x + 1.)
                  /(((((B5*x + B4)*x + B3)*x + B2)*x + B1)*x + 1.);
        }
    } else if (qsq < 0.8) {
        const double x = qsq;
        const double C0 = +1.;
        const double C1 = -1./3.;
        const double C2 = +1./12.;
        const double C3 = -1./60.;
        const double C4 = +1./360.;
        const double C5 = -1./2520.;
        const double C6 = +1./20160.;
        const double C7 = -1./181440.;
        //return ((((C5*x + C4)*x + C3)*x + C2)*x + C1)*x + C0;
        //return (((((C6*x + C5)*x + C4)*x + C3)*x + C2)*x + C1)*x + C0;
        return ((((((C7*x + C6)*x + C5)*x + C4)*x + C3)*x + C2)*x + C1)*x + C0;
    } else {
        return 2.*(expm1(-qsq) + qsq)/(qsq*qsq);
    }
    """, "sas_debye"),
)

def mp_star_polymer(x, arms=3): # x = q*Rg
    from mpmath import expm1
    v = x * arms / (3*arms - 2)
    term1 = v + expm1(-v)
    term2 = (arms- 1)/2 * expm1(-v)**2
    return 2 * (term1 + term2) / (arms * v**2) if v > 0 else 1

def np_star_polymer(x, arms=3):
    from numpy import expm1, polyval
    scalar = np.isscalar(x)
    if scalar:
        x = np.array([x]) # should inherit dtype for single if given single
    #T1 = [1, -1/3, 1/12, -1/60, 1/360, -1/2520, 1/20160, -1/181440][::-1]
    T1 = [1, -1/3, 1/12, -1/60, 1/360, -1/2520][::-1]
    #T2 = [1, -1, 7/12, -1/4, 31/360, -1/40][::-1]
    f = np.empty_like(x)
    cutoff = 0.03 if f.dtype == np.float64 else 1.0
    index = (x == 0.)
    f[index] = 1.0
    index = ~index & (x < cutoff)
    v = x * arms / (3*arms - 2)
    vi = v[index]
    #f[index] = polyval(T1, vi)/arms + (1 - 1/arms)*polyval(T2, vi)
    #f[index] = 2/(arms*vi)*(1 + expm1(-vi)/vi) + (1 - 1/arms)*polyval(T2, vi)
    f[index] = polyval(T1, vi)/arms + (1 - 1/arms)*(expm1(-vi)/vi)**2
    #f[index] = 2/(arms*vi)*(1 + expm1(-vi)/vi) + (1 - 1/arms)*(expm1(-vi)/vi)**2

    term1 = v[~index] + expm1(-v[~index])
    term2 = (arms - 1)/2 * expm1(-v[~index])**2
    f[~index] = 2 * (term1 + term2) / (arms * v[~index]**2)
    return f

ocl_star_polymer = """
    const int arms = 3;
    const double v = q * arms / (3.0 * arms - 2.0);
// Note: cutoff values are for (Q Rg)^2.
#if FLOAT_SIZE>4
#define STAR_POLYMER_CUTOFF 0.03
#else
#define STAR_POLYMER_CUTOFF 1.0
#endif

    if (q == 0.) {
        return 1.;
    } else if (q <= STAR_POLYMER_CUTOFF) {
        double P1 = 1. + v*(-1./3. + v*(1./12. + v*(-1./60. + v*(1./360. + v*(-1./2520)))));
        //double P2 = 1. + v*(-1. + v*(7./12. - v*(-1./4.)));
        return P1/arms + (1. - 1./arms)*square(expm1(-v)/v);
    } else {
        double term1 = v + expm1(-v);
        double term2 = ((arms - 1.0)/2.0) * square(expm1(-v));
        return (2.0 * (term1 + term2)) / (arms * v * v);
    }
"""

add_function(
    name="star_polymer(arms=3)",
    mp_function=mp_star_polymer,
    np_function=np_star_polymer,
    ocl_function=make_ocl(ocl_star_polymer, "star_polymer", []),
    shortname="star_polymer",
    xaxis="$(Q R_g)^2$ (unitless)",
)


RADIUS=3000
LENGTH=30
THETA=45
def mp_cyl(x):
    f = mp.mpf
    theta = f(THETA)*mp.pi/f(180)
    qr = x * f(RADIUS)*mp.sin(theta)
    qh = x * f(LENGTH)/f(2)*mp.cos(theta)
    be = f(2)*mp.j1(qr)/qr
    si = mp.sin(qh)/qh
    background = f(0)
    #background = f(1)/f(1000)
    volume = mp.pi*f(RADIUS)**f(2)*f(LENGTH)
    contrast = f(5)
    units = f(1)/f(10000)
    #return be
    #return si
    return units*(volume*contrast*be*si)**f(2)/volume + background
def np_cyl(x):
    f = np.float64 if x.dtype == np.float64 else np.float32
    theta = f(THETA)*f(np.pi)/f(180)
    qr = x * f(RADIUS)*np.sin(theta)
    qh = x * f(LENGTH)/f(2)*np.cos(theta)
    be = f(2)*scipy.special.j1(qr)/qr
    si = np.sin(qh)/qh
    background = f(0)
    #background = f(1)/f(1000)
    volume = f(np.pi)*f(RADIUS)**2*f(LENGTH)
    contrast = f(5)
    units = f(1)/f(10000)
    #return be
    #return si
    return units*(volume*contrast*be*si)**f(2)/volume + background
ocl_cyl = """\
    double THETA = %(THETA).15e*M_PI_180;
    double qr = q*%(RADIUS).15e*sin(THETA);
    double qh = q*0.5*%(LENGTH).15e*cos(THETA);
    double be = sas_2J1x_x(qr);
    double si = sas_sinx_x(qh);
    double background = 0;
    //double background = 0.001;
    double volume = M_PI*square(%(RADIUS).15e)*%(LENGTH).15e;
    double contrast = 5.0;
    double units = 1e-4;
    //return be;
    //return si;
    return units*square(volume*contrast*be*si)/volume + background;
"""%{"LENGTH":LENGTH, "RADIUS": RADIUS, "THETA": THETA}
add_function(
    name="cylinder(r=%g, l=%g, theta=%g)"%(RADIUS, LENGTH, THETA),
    mp_function=mp_cyl,
    np_function=np_cyl,
    ocl_function=make_ocl(ocl_cyl, "ocl_cyl", ["lib/polevl.c", "lib/sas_J1.c"]),
    shortname="cylinder",
    xaxis="$q/A^{-1}$",
)

lanczos_gamma = """\
    const double coeff[] = {
            76.18009172947146, -86.50532032941677,
            24.01409824083091, -1.231739572450155,
            0.1208650973866179e-2,-0.5395239384953e-5
            };
    const double x = q;
    double tmp  = x + 5.5;
    tmp -= (x + 0.5)*log(tmp);
    double ser = 1.000000000190015;
    for (int k=0; k < 6; k++) ser += coeff[k]/(x + k+1);
    return -tmp + log(2.5066282746310005*ser/x);
"""
add_function(
    name="loggamma(x)",
    mp_function=mp.loggamma,
    np_function=scipy.special.gammaln,
    ocl_function=make_ocl(lanczos_gamma, "lgamma"),
)

replacement_expm1 = """\
      double x = (double)q;  // go back to float for single precision kernels
      // Adapted from the cephes math library.
      // Copyright 1984 - 1992 by Stephen L. Moshier
      if (x != x || x == 0.0) {
         return x; // NaN and +/- 0
      } else if (x < -0.5 || x > 0.5) {
         return exp(x) - 1.0;
      } else {
         const double xsq = x*x;
         const double p = (((
            +1.2617719307481059087798E-4)*xsq
            +3.0299440770744196129956E-2)*xsq
            +9.9999999999999999991025E-1);
         const double q = ((((
            +3.0019850513866445504159E-6)*xsq
            +2.5244834034968410419224E-3)*xsq
            +2.2726554820815502876593E-1)*xsq
            +2.0000000000000000000897E0);
         double r = x * p;
         r =  r / (q - r);
         return r+r;
       }
"""
add_function(
    name="sas_expm1(x)",
    mp_function=mp.expm1,
    np_function=np.expm1,
    ocl_function=make_ocl(replacement_expm1, "sas_expm1"),
)

# Alternate versions of 3 j1(x)/x, for posterity
def taylor_3j1x_x(x):
    """
    Calculation using taylor series.
    """
    # Generate coefficients using the precision of the target value.
    n = 5
    cinv = [3991680, -45360, 840, -30, 3]
    three = x.dtype.type(3)
    p = three/np.array(cinv, x.dtype)
    return np.polyval(p[-n:], x*x)
add_function(
    name="3 j1(x)/x: taylor",
    mp_function=lambda x: 3*(mp.sin(x)/x - mp.cos(x))/(x*x),
    np_function=taylor_3j1x_x,
    ocl_function=make_ocl("return sas_3j1x_x(q);", "sas_j1c", ["lib/sas_3j1x_x.c"]),
)
def trig_3j1x_x(x):
    r"""
    Direct calculation using linear combination of sin/cos.

    Use the following trig identity:

    .. math::

        a \sin(x) + b \cos(x) = c \sin(x + \phi)

    where $c = \surd(a^2+b^2)$ and $\phi = \tan^{-1}(b/a) to calculate the
    numerator $\sin(x) - x\cos(x)$.
    """
    one = x.dtype.type(1)
    three = x.dtype.type(3)
    c = np.sqrt(one + x*x)
    phi = np.arctan2(-x, one)
    return three*(c*np.sin(x+phi))/(x*x*x)
add_function(
    name="3 j1(x)/x: trig",
    mp_function=lambda x: 3*(mp.sin(x)/x - mp.cos(x))/(x*x),
    np_function=trig_3j1x_x,
    ocl_function=make_ocl("return sas_3j1x_x(q);", "sas_j1c", ["lib/sas_3j1x_x.c"]),
)
def np_2J1x_x(x):
    """
    numpy implementation of 2J1(x)/x using single precision algorithm
    """
    # pylint: disable=bad-continuation
    f = x.dtype.type
    ax = abs(x)
    if ax < f(8.0):
        y = x*x
        ans1 = f(2)*(f(72362614232.0)
                  + y*(f(-7895059235.0)
                  + y*(f(242396853.1)
                  + y*(f(-2972611.439)
                  + y*(f(15704.48260)
                  + y*(f(-30.16036606)))))))
        ans2 = (f(144725228442.0)
                  + y*(f(2300535178.0)
                  + y*(f(18583304.74)
                  + y*(f(99447.43394)
                  + y*(f(376.9991397)
                  + y)))))
        return ans1/ans2
    else:
        y = f(64.0)/(ax*ax)
        xx = ax - f(2.356194491)
        ans1 = (f(1.0)
                  + y*(f(0.183105e-2)
                  + y*(f(-0.3516396496e-4)
                  + y*(f(0.2457520174e-5)
                  + y*f(-0.240337019e-6)))))
        ans2 = (f(0.04687499995)
                  + y*(f(-0.2002690873e-3)
                  + y*(f(0.8449199096e-5)
                  + y*(f(-0.88228987e-6)
                  + y*f(0.105787412e-6)))))
        sn, cn = np.sin(xx), np.cos(xx)
        ans = np.sqrt(f(0.636619772)/ax) * (cn*ans1 - (f(8.0)/ax)*sn*ans2) * f(2)/x
        return -ans if (x < f(0.0)) else ans
add_function(
    name="2 J1(x)/x:alt",
    mp_function=lambda x: 2*mp.j1(x)/x,
    np_function=lambda x: np.asarray([np_2J1x_x(v) for v in x], x.dtype),
    ocl_function=make_ocl("return sas_2J1x_x(q);", "sas_2J1x_x", ["lib/polevl.c", "lib/sas_J1.c"]),
)

ALL_FUNCTIONS = set(FUNCTIONS.keys())
ALL_FUNCTIONS.discard("loggamma")  # use cephes-based gammaln instead
ALL_FUNCTIONS.discard("3j1/x:taylor")
ALL_FUNCTIONS.discard("3j1/x:trig")
ALL_FUNCTIONS.discard("2J1/x:alt")

# =============== MAIN PROGRAM ================

def usage():
    names = ", ".join(sorted(ALL_FUNCTIONS))
    print("""\
usage: precision.py [-f/a/r] [-x<range>] "name" ...
where
    -f indicates that the function value should be plotted,
    -a indicates that the absolute error should be plotted,
    -r indicates that the relative error should be plotted (default),
    -x<range> indicates the steps in x, where <range> is one of the following
        log indicates log stepping in [10^-3, 10^5] (default)
        logq indicates log stepping in [10^-4, 10^1]
        linear indicates linear stepping in [1, 1000]
        zoom indicates linear stepping in [1000, 1010]
        neg indicates linear stepping in [-100.1, 100.1]
        start:stop:n[:stepping] indicates an n-step plot in [start, stop]
            or [10^start, 10^stop] if stepping is "log" (default n=400)
Some functions (notably gammainc/gammaincc) have an additional parameter A
which can be set from the command line as A=value.  Default is A=1.

Name is one of:
    """+names)
    sys.exit(1)

def main():
    import sys
    diff = "relative"
    xrange = "log"
    options = [v for v in sys.argv[1:] if v.startswith('-')]
    for opt in options:
        if opt == '-f':
            diff = "none"
        elif opt == '-r':
            diff = "relative"
        elif opt == '-a':
            diff = "absolute"
        elif opt.startswith('-x'):
            xrange = opt[2:]
        else:
            usage()

    names = [v for v in sys.argv[1:] if not v.startswith('-')]
    if not names:
        usage()

    if names[0] == "all":
        cutoff = names[1] if len(names) > 1 else ""
        names = list(sorted(ALL_FUNCTIONS))
        names = [k for k in names if k >= cutoff]
    if any(k not in FUNCTIONS for k in names):
        usage()
    multiple = len(names) > 1
    pylab.interactive(multiple)
    for k in names:
        pylab.clf()
        comparator = FUNCTIONS[k]
        comparator.run(xrange=xrange, diff=diff)
        if multiple:
            input()
    if not multiple:
        pylab.show()

if __name__ == "__main__":
    main()
