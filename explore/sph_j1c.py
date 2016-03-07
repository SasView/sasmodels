r"""
Show numerical precision of $3 j_1(x)/x$ vs. its Taylor series expansion.

The choice of the number of terms in the series and the cutoff value for
switching between series and direct calculation depends on the numeric
precision.

Point where direct calculation reaches machine precision::

    single machine precision eps 3e-8 at qr=1.1 (see below)
    double machine precision eps 4e-16 at qr=1.1

Point where Taylor series reaches machine precision (eps), where taylor
series matches direct calculation (cross) and the error at that point::

    prec   n eps  cross  error
    single 3 0.28  0.4   6.2e-7
    single 4 0.68  0.7   2.3e-7
    single 5 1.18  1.2   7.5e-8
    double 3 0.01  0.03  2.3e-13
    double 4 0.06  0.1   3.1e-14
    double 5 0.16  0.2   5.0e-15

Note: relative error on single precision starts increase on the direct
method at qr=1.1, rising from 3e-8 to 5e-5 by qr=1e3.  This should be
safe for the sans range, with objects of 100 nm supported to a q of 0.1
while maintaining 5 digits of precision.  For usans/sesans, the objects
are larger but the q is smaller, so again it should be fine.
"""

import numpy as np
from sympy.mpmath import mp
#import matplotlib; matplotlib.use('TkAgg')
import pylab

mp.dps = 150 # number of digits to use in estimating true value

SHOW_DIFF = True  # True if show diff rather than function value
LINEAR_X = False  # True if q is linearly spaced instead of log spaced

def mp_j1c(vec):
    """
    Direct calculation using sympy multiprecision library.
    """
    return [_mp_j1c(mp.mpf(x)) for x in vec]

def _mp_j1c(x):
    """
    Helper funciton for mp_j1c
    """
    return mp.mpf(3)*(mp.sin(x)/x - mp.cos(x))/(x*x)

def np_j1c(x, dtype):
    """
    Direct calculation using numpy.
    """
    x = np.asarray(x, dtype)
    return np.asarray(3, dtype)*(np.sin(x) - x*np.cos(x))/(x*x*x)

def lin_j1c(x, dtype):
    r"""
    Direct calculation using linear combination of sin/cos.

    Use the following trig identity:

    .. math::

        a \sin(x) + b \cos(x) = c \sin(x + \phi)

    where $c = \surd(a^2+b^2)$ and $\phi = \tan^{-1}(b/a) to calculate the
    numerator $\sin(x) - x\cos(x)$.
    """
    x = np.asarray(x, dtype)
    c = np.sqrt(np.asarray(1,dtype) + x*x)
    phi = np.arctan2(-np.asarray(x,dtype),np.asarray(1,dtype))
    return np.asarray(3, dtype)*(c*np.sin(x+phi))/(x*x*x)

def taylor_j1c(x, dtype, n):
    """
    Calculation using taylor series.
    """
    # Generate coefficients using the precision of the target value.
    cinv = [3991680, -45360, 840, -30, 3]
    x = np.asarray(x, dtype)
    p = np.asarray(3, dtype)/np.array(cinv, dtype)
    return np.polyval(p[-n:], x*x)

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
    target = np.asarray(mp_j1c(x), 'double')
    direct = np_j1c(x, precision)
    comb = lin_j1c(x, precision)
    taylor3 = taylor_j1c(x, precision, 3)
    taylor4 = taylor_j1c(x, precision, 4)
    taylor5 = taylor_j1c(x, precision, 5)
    plotdiff(x, target, direct, 'direct '+precision)
    #plotdiff(x, target, comb, 'c sin(x+phi) '+precision)
    plotdiff(x, target, taylor3, 'taylor 3 '+precision)
    plotdiff(x, target, taylor4, 'taylor 4 '+precision)
    plotdiff(x, target, taylor5, 'taylor 5 '+precision)
    pylab.xlabel("qr (1/Ang)")
    if SHOW_DIFF:
        pylab.ylabel("relative error")
    else:
        pylab.ylabel("3 j1(x)/x")
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
    pylab.suptitle('3 j1(x)/x')

if __name__ == "__main__":
    #print "\n".join(str(x) for x in mp_j1c([1e-6,1e-5,1e-4,1e-3]))
    main()
    pylab.show()
