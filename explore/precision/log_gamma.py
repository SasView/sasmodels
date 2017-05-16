r"""
Show numerical precision of $\ln \Gamma(x)$.
"""

import numpy as np
import scipy.special
from sympy.mpmath import mp
#import matplotlib; matplotlib.use('TkAgg')
import pylab

mp.dps = 150 # number of digits to use in estimating true value

SHOW_DIFF = True  # True if show diff rather than function value
LINEAR_X = False  # True if q is linearly spaced instead of log spaced

def mp_gamma(vec):
    """
    Direct calculation using sympy multiprecision library.
    """
    return [_mp_fn(mp.mpf(x)) for x in vec]

def _mp_fn(x):
    """
    Helper funciton for mp_j1c
    """
    #return mp.gamma(x)
    return mp.loggamma(x)

def np_gamma(x, dtype):
    """
    Direct calculation using scipy.
    """
    x = np.asarray(x, dtype)
    return scipy.special.gammaln(x)
    #return scipy.special.gamma(x)

def lanczos_gamma(x, dtype):
    coeff = np.asarray((
         76.18009172947146,     -86.50532032941677,
         24.01409824083091,     -1.231739572450155,
          0.1208650973866179e-2,-0.5395239384953e-5), dtype)

    x = np.asarray(x, dtype)
    tmp  = x + np.asarray(5.5, dtype)
    tmp -= (x + np.asarray(0.5, dtype))*np.log(tmp)
    ser  = np.ones_like(x)*np.asarray(1.000000000190015, dtype)
    for k,c in enumerate(coeff):
        ser += c/(x + np.asarray(k+1, dtype))
    return -tmp+np.log(np.asarray(2.5066282746310005, dtype)*ser/x);

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
        pylab.loglog(x, np.clip(actual,*limits),  '-', label=label)

def compare(x, precision):
    r"""
    Compare the different computation methods using the given precision.
    """
    target = np.asarray(mp_gamma(x), 'double')
    direct = np_gamma(x, precision)
    approx = lanczos_gamma(x, precision)
    plotdiff(x, target, direct, 'scipy '+precision)
    plotdiff(x, target, approx, 'sasmodels '+precision)
    pylab.xlabel("x (arbitrary)")
    if SHOW_DIFF:
        pylab.ylabel("relative error")
    else:
        pylab.ylabel("ln(gamma(x))")
        pylab.loglog(x, target,  '-', label="true value")
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
    pylab.suptitle('ln gamma')

if __name__ == "__main__":
    main()
    pylab.show()
