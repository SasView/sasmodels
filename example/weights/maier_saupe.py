from __future__ import print_function
import numpy as np
from numpy import exp, sin, degrees, radians, pi, sqrt

from sasmodels.weights import Dispersion as BaseDispersion

class Dispersion(BaseDispersion):
    r"""
    Maier-Saupe dispersion on orientation.

    .. math:

        w(\theta) = e^{a{\cos^2 \theta}}

    This provides a close match to the gaussian distribution for
    low angles, but the tails are limited to $\pm 90^\circ$.  For $a \ll 1$
    the distribution is approximately uniform.  The usual polar coordinate
    projection applies, with $\theta$ weights scaled by $\cos \theta$
    and $\phi$ weights unscaled.

    This is equivalent to a cyclic gaussian distribution
    $w(\theta) = e^{-sin^2(\theta)/(2\a^2)}.

    Note that the incorrect distribution is used for $a=0$.  With the
    distribution parameter labelled as *width*, the value *width=0* is
    is assumed to be completely oriented and no polydispersity distribution
    is generated.  However a value of $a=0$ should be completely unoriented.
    Any small value (e.g., 0.01 or lower) should suffice.

    The order parameter $P_2$ is defined as

    .. math:

        P(a, \beta) = \frac{e^{a \cos^2 \beta}}{
                       4\pi \int_0^{\pi/2} e^{a \cos^2 \beta} \sin \beta\,d\beta}

        P_2 = 4\pi \int_0^{\pi/2} \frac{3 \cos^2 \beta - 1)}{2}
                        P(a, \beta) \sin \beta\,d\beta

    where $a$ is the distribution width $\sigma$ handed to the weights function.
    There is a somewhat complicated closed form solution

    .. math:

        P_2 = \frac{3e^a}{2\sqrt{\pi a} E_a} - \frac{3}{4a} - \frac{1}{2}

    where $E_a = \mathrm{erfi}(\sqrt{a})$ is the imaginary error function,
    which can be coded in python as::

        from numpy import pi, sqrt, exp
        from scipy.special import erfi

        def P_2(a):
            # protect against overflow with a Taylor series at infinity
            if a <= 700:
                r = exp(a)/sqrt(pi*a)/erfi(sqrt(a))
            else:
                r = 1/((((6.5525/a + 1.875)/a + 0.75)/a + 0.5)/a + 1)
            return 1.5*r - 0.75/a - 0.5

    Given an order parameter $S = P_2(a)$, one can also solve for the
    equivalent $a$:

        from scipy.optimize import fsolve

        def P_2_inv(S):
            return fsolve(lambda x: P_2(x) - S, 1.0)[0]

    References
    ----------

    [1] Hardouin, et al., 1995. SANS study of a semiflexible main chain
    liquid crystalline polyether. *Macromolecules* 28, 5427-5433.
    """
    type = "maier_saupe"
    default = dict(npts=35, width=1, nsigmas=None)

    # Note: center is always zero for orientation distributions
    def _weights(self, center, sigma, lb, ub):
        # use the width parameter as the value for Maier-Saupe "a"
        # and find the equivalent width sigma
        a = sigma
        sigma = 1./sqrt(2.*a)

        # Limit width to +/- 90 degrees
        width = min(self.nsigmas*sigma, pi/2)
        x = np.linspace(-width, width, self.npts)

        # Truncate the distribution in case the parameter value is limited
        x[(x >= radians(lb)) & (x <= radians(ub))]

        # Return orientation in degrees with Maier-Saupe weights
        # Note: weights are normalized to sum to 1, so we can scale
        # by an arbitrary scale factor c = exp(m) to get:
        #     w = exp(m*cos(x)**2)/c = exp(-m*sin(x)**2)
        return degrees(x), exp(-a*sin(x)**2)



def P_2(a):
    from numpy import pi, sqrt, exp
    from scipy.special import erfi

    # Double precision e^x overflows so do a Taylor expansion around infinity
    # for large values of a.  With five terms it is accurate to 12 digits
    # at a = 700, and 7 digits at a = 75.
    if a <= 700:
        r = exp(a)/sqrt(pi*a)/erfi(sqrt(a))
    else:
        r = 1/((((6.5525/a + 1.875)/a + 0.75)/a + 0.5)/a + 1)
    return 1.5*r - 0.75/a - 0.5

def P_2_inv(S):
    from scipy.optimize import fsolve
    return fsolve(lambda x: P_2(x) - S, 1.0)[0]

def P_2_numerical(a):
    from scipy.integrate import romberg
    from numpy import cos, sin, pi, exp
    def num(beta):
        return (1.5 * cos(beta)**2 - 0.5) * exp(a * cos(beta)**2) * sin(beta)
    def denom(beta):
        return exp(a * cos(beta)**2) * sin(beta)
    return romberg(num, 0, pi/2) / romberg(denom, 0, pi/2)

if __name__ == "__main__":
    import sys
    a = float(sys.argv[1])
    sigma = 1/(2*radians(a)**2)
    #print("P_2", P_2(a), "difference from integral", P_2(a) - P_2_numerical(a))
    print("a=%g, sigma=%g, P_2=%g, P_2_inv(P_2(a))-a=%g"
          % (a, sigma, P_2(a), P_2_inv(P_2(a))-a))
