r"""
Show numerical precision of $2 J_1(x)/x$.
"""

import numpy as np
from sympy.mpmath import mp
#import matplotlib; matplotlib.use('TkAgg')
import pylab


SHOW_DIFF = True # True if show diff rather than function value
#SHOW_DIFF = False # True if show diff rather than function value
LINEAR_X = False  # True if q is linearly spaced instead of log spaced
#LINEAR_X = True # True if q is linearly spaced instead of log spaced
FUNCTION = "2*J1(x)/x"

def mp_fn(vec, bits=500):
    """
    Direct calculation using sympy multiprecision library.
    """
    with mp.workprec(bits):
        return [_mp_fn(mp.mpf(x)) for x in vec]

def _mp_fn(x):
    """
    Actual function that gets evaluated.  The caller just vectorizes.
    """
    return mp.mpf(2)*mp.j1(x)/x

def np_fn(x, dtype):
    """
    Direct calculation using scipy.
    """
    from scipy.special import j1 as J1
    x = np.asarray(x, dtype)
    return np.asarray(2, dtype)*J1(x)/x

def sasmodels_fn(x, dtype, platform='ocl'):
    """
    Calculation using pade approximant.
    """
    from sasmodels import core, data, direct_model
    model = core.load_model('bessel', dtype=dtype)
    calculator = direct_model.DirectModel(data.empty_data1D(x), model)
    return calculator(background=0)

def plotdiff(x, target, actual, label):
    """
    Plot the computed value.

    Use relative error if SHOW_DIFF, otherwise just plot the value directly.
    """
    if SHOW_DIFF:
        err = abs((target-actual)/target)
        #err = np.clip(err, 0, 1)
        pylab.loglog(x, err, '-', label=label)
    else:
        limits = np.min(target), np.max(target)
        pylab.semilogx(x, np.clip(actual,*limits),  '-', label=label)

def compare(x, precision, target):
    r"""
    Compare the different computation methods using the given precision.
    """
    #plotdiff(x, target, mp_fn(x, 11), 'mp 11 bits')
    plotdiff(x, target, np_fn(x, precision), 'numpy '+precision)
    plotdiff(x, target, sasmodels_fn(x, precision, 0), 'sasmodels '+precision)
    pylab.xlabel("qr (1/Ang)")
    if SHOW_DIFF:
        pylab.ylabel("relative error")
    else:
        pylab.ylabel(FUNCTION)
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
    target = np.asarray(mp_fn(qr), 'double')
    pylab.subplot(121)
    compare(qr, 'single', target)
    pylab.legend(loc='best')
    pylab.subplot(122)
    compare(qr, 'double', target)
    pylab.legend(loc='best')
    pylab.suptitle(FUNCTION)

if __name__ == "__main__":
    #print "\n".join(str(x) for x in mp_J1c([1e-6,1e-5,1e-4,1e-3]))
    main()
    pylab.show()
