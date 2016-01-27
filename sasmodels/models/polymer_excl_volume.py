r"""
This model describes the scattering from polymer chains subject to excluded
volume effects and has been used as a template for describing mass fractals.

Definition
----------

The form factor was originally presented in the following integral form
(Benoit, 1957)

.. math::

    P(Q)=2\int_0^{1}dx(1-x)exp\left[-\frac{Q^2a^2}{6}n^{2v}x^{2v}\right]

where $\nu$ is the excluded volume parameter
(which is related to the Porod exponent $m$ as $\nu=1/m$ ),
$a$ is the statistical segment length of the polymer chain,
and $n$ is the degree of polymerization.
This integral was later put into an almost analytical form as follows
(Hammouda, 1993)

.. math::

    P(Q)=\frac{1}{\nu U^{1/2\nu}}\gamma\left(\frac{1}{2\nu},U\right) -
    \frac{1}{\nu U^{1/\nu}}\gamma\left(\frac{1}{\nu},U\right)

where $\gamma(x,U)$ is the incomplete gamma function

.. math::

    \gamma(x,U)=\int_0^{U}dt\ exp(-t)t^{x-1}

and the variable $U$ is given in terms of the scattering vector $Q$ as

.. math::

    U=\frac{Q^2a^2n^{2\nu}}{6} = \frac{Q^2R_{g}^2(2\nu+1)(2\nu+2)}{6}

The square of the radius-of-gyration is defined as

.. math::

    R_{g}^2 = \frac{a^2n^{2\nu}}{(2\nu+1)(2\nu+2)}

Note that this model applies only in the mass fractal range (ie, $5/3<=m<=3$ )
and **does not apply** to surface fractals ( $3<m<=4$ ).
It also does not reproduce the rigid rod limit (m=1) because it assumes chain
flexibility from the outset. It may cover a portion of the semi-flexible chain
range ( $1<m<5/3$ ).

A low-Q expansion yields the Guinier form and a high-Q expansion yields the
Porod form which is given by

.. math::

    P(Q\rightarrow \infty) = \frac{1}{\nu U^{1/2\nu}}\Gamma\left(
    \frac{1}{2\nu}\right) - \frac{1}{\nu U^{1/\nu}}\Gamma\left(
    \frac{1}{\nu}\right)

Here $\Gamma(x) = \gamma(x,\infty)$ is the gamma function.

The asymptotic limit is dominated by the first term

.. math::

    P(Q\rightarrow \infty) \sim \frac{1}{\nu U^{1/2\nu}}\Gamma\left(\frac{1}{2\nu}\right) =
    \frac{m}{\left(QR_{g}\right)^m}\left[\frac{6}{(2\nu +1)(2\nu +2)} \right]^{m/2}
    \Gamma (m/2)

The special case when $\nu=0.5$ (or $m=2/\nu=2$ ) corresponds to Gaussian chains for
which the form factor is given by the familiar Debye function.

.. math::

    P(Q) = \frac{2}{Q^4R_{g}^4} \left[exp(-Q^2R_{g}^2) - 1 + Q^2R_{g}^2 \right]

For 2D data: The 2D scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}

This example dataset is produced using 200 data points, $qmin=0.001Ang^{-1}$,
$qmax=0.2Ang^{-1}$ and the default values

.. figure:: img/polymer_excl_volume_1d.jpg

    1D plot using the default values (w/500 data point).


References
----------

H Benoit, *Comptes Rendus*, 245 (1957) 2244-2247

B Hammouda, *SANS from Homogeneous Polymer Mixtures - A Unified Overview,
Advances in Polym. Sci.* 106(1993) 87-133

"""

from math import sqrt
from numpy import inf, power
from scipy.special import gammainc, gamma

name = "polymer_excl_volume"
title = "Polymer Excluded Volume model"
description = """Compute the scattering intensity from polymers with excluded
                volume effects.
                rg:         radius of gyration
                porod_exp:  Porod exponent
              """
category = "shape-independent"

#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["rg",        "Ang", 60.0, [0, inf],    "", "Radius of Gyration"],
              ["porod_exp", "",     3.0, [-inf, inf], "", "Porod exponent"],
              ]


def Iq(q, rg, porod_exp):

    """
    :param q:         Input q-value (float or [float, float])
    :param rg:        Radius of gyration
    :param porod_exp: Porod exponent
    :return:          Calculated intensity
    """
    nu = 1.0/porod_exp
    u = q*q*rg*rg*(2.0*nu+1.0) * (2.0*nu+2.0)/6.0
    o2nu = 1.0/(2.0*nu)

    intensity = ((1.0/(nu*power(u, o2nu))) * (gamma(o2nu)*gammainc(o2nu, u) -
                  1.0/power(u, o2nu) * gamma(porod_exp) *
                  gammainc(porod_exp, u))) * (q > 0) + 1.0*(q <= 0)

    return intensity

Iq.vectorized = True  # Iq accepts an array of q values


def Iqxy(qx, qy, *args):
        iq = Iq(sqrt(qx**2 + qy**2), *args)

        return iq

Iqxy.vectorized = True  # Iqxy accepts an array of qx, qy values


demo = dict(scale=1, background=0.0,
            rg=60.0,
            porod_exp=3.0)

oldname = "PolymerExclVolume"
oldpars = dict(background='background', scale='scale',
               rg='rg',
               porod_exp='m')

tests = [
         # Accuracy tests based on content in test/polyexclvol_default_igor.txt
         [{'rg': 60, 'porod_exp': 3.0}, 0.001, 0.998801],
         [{'rg': 60, 'porod_exp': 3.0}, 0.105363, 0.0162751],
         [{'rg': 60, 'porod_exp': 3.0}, 0.665075, 6.56261e-05],

         # Additional tests with larger range of parameters
         [{'rg': 10, 'porod_exp': 4.0}, 0.1, 0.723436675809],
         [{'rg': 2.2, 'porod_exp': 22.0, 'background': 100.0}, 5.0, 100.0],
         [{'rg': 1.1, 'porod_exp': 1, 'background': 10.0, 'scale': 1.25},
         20000., 10.0000712097]
         ]
