r"""
Show numerical precision of cylinder form.

Using::

    qr = q r sin(t)
    qh = q h/2 cos(t)
    F = 2 J_1(qr)/qr sin(qh)/qh
"""

import numpy as np
from sympy.mpmath import mp
#import matplotlib; matplotlib.use('TkAgg')
import pylab

SHOW_DIFF = True  # True if show diff rather than function value
LINEAR_X = False  # True if q is linearly spaced instead of log spaced

RADIUS = 20
LENGTH = 300
CONTRAST = 5
THETA = 45

def mp_form(vec, bits=500):
    """
    Direct calculation using sympy multiprecision library.
    """
    with mp.workprec(bits):
        return [_mp_f(mp.mpf(x)) for x in vec]

def _mp_f(x):
    """
    Helper function for mp_j1c
    """
    f = mp.mpf
    theta = f(THETA)*mp.pi/f(180)
    qr = x * f(RADIUS)*mp.sin(theta)
    qh = x * f(LENGTH)/f(2)*mp.cos(theta)
    return (f(2)*mp.j1(qr)/qr * mp.sin(qh)/qh)**f(2)

def np_form(x, dtype):
    """
    Direct calculation using scipy.
    """
    from scipy.special import j1
    f = np.float64 if np.dtype(dtype) == np.float64 else np.float32
    x = np.asarray(x, dtype)
    theta = f(THETA)*f(np.pi)/f(180)
    qr = x * f(RADIUS)*np.sin(theta)
    qh = x * f(LENGTH)/f(2)*np.cos(theta)
    return (f(2)*j1(qr)/qr*np.sin(qh)/qh)**f(2)

def sasmodels_form(x, dtype):
    f = np.float64 if np.dtype(dtype) == np.float64 else np.float32
    x = np.asarray(x, dtype)
    theta = f(THETA)*f(np.pi)/f(180)
    qr = x * f(RADIUS)*np.sin(theta)
    qh = x * f(LENGTH)/f(2)*np.cos(theta)
    return (J1c(qr, dtype)*np.sin(qh)/qh)**f(2)

def J1c(x, dtype):
    x = np.asarray(x, dtype)
    f = np.float64 if np.dtype(dtype) == np.float64 else np.float32
    return np.asarray([_J1c(xi, f) for xi in x], dtype)

def _J1c(x, f):
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

def plotdiff(x, target, actual, label):
    """
    Plot the computed value.

    Use relative error if SHOW_DIFF, otherwise just plot the value directly.
    """
    if SHOW_DIFF:
        err = np.clip(abs((target-actual)/target), 0, 1)
        pylab.loglog(x, err, '-', label=label)
    else:
        limits = np.min(target), np.max(target)
        pylab.semilogx(x, np.clip(actual,*limits),  '-', label=label)

def compare(x, precision):
    r"""
    Compare the different computation methods using the given precision.
    """
    target = np.asarray(mp_form(x), 'double')
    plotdiff(x, target, mp_form(x, bits=11), '11-bit')
    plotdiff(x, target, np_form(x, precision), 'direct '+precision)
    plotdiff(x, target, sasmodels_form(x, precision), 'sasmodels '+precision)
    pylab.xlabel("qr (1/Ang)")
    if SHOW_DIFF:
        pylab.ylabel("relative error")
    else:
        pylab.ylabel("2 J1(x)/x")
        pylab.semilogx(x, target,  '-', label="true value")
    if LINEAR_X:
        pylab.xscale('linear')

def main():
    r"""
    Compare accuracy of different methods for computing $3 j_1(x)/x$.
    :return:
    """
    if LINEAR_X:
        qr = np.linspace(1,1000,2000)
    else:
        qr = np.logspace(-3,5,400)
    pylab.subplot(121)
    compare(qr, 'single')
    pylab.legend(loc='best')
    pylab.subplot(122)
    compare(qr, 'double')
    pylab.legend(loc='best')
    pylab.suptitle('2 J1(x)/x')

if __name__ == "__main__":
    #print "\n".join(str(x) for x in mp_J1c([1e-6,1e-5,1e-4,1e-3]))
    main()
    pylab.show()
