r"""
Show numerical precision of $sin(x)/x$.
"""

import numpy as np
from sympy.mpmath import mp
#import matplotlib; matplotlib.use('TkAgg')
import pylab

mp.dps = 150 # number of digits to use in estimating true value

SHOW_DIFF = True  # True if show diff rather than function value
LINEAR_X = False  # True if q is linearly spaced instead of log spaced

def mp_sinc(vec):
    """
    Direct calculation using sympy multiprecision library.
    """
    return [_mp_sinc(mp.mpf(x)) for x in vec]

def _mp_sinc(x):
    """
    Helper funciton for mp_j1c
    """
    return mp.sin(x)/x

def np_sinc(x, dtype):
    """
    Direct calculation using scipy.
    """
    x = np.asarray(x, dtype)
    return np.sin(x)/x
    #return np.asarray(np.sin(np.double(x))/np.double(x),dtype)

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
    target = np.asarray(mp_sinc(x), 'double')
    direct = np_sinc(x, precision)
    plotdiff(x, target, direct, 'direct '+precision)
    pylab.xlabel("qr (1/Ang)")
    if SHOW_DIFF:
        pylab.ylabel("relative error")
    else:
        pylab.ylabel("sin(x)/x")
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
    pylab.suptitle('sin(x)/x')

if __name__ == "__main__":
    print "\n".join(str(x) for x in mp_sinc([1e-6,1e-5,1e-4,1e-3]))
    main()
    pylab.show()
