r"""
This model calculates the structure factor, $S(q)$, of one-component liquid
systems interacting with a two-term Yukawa potential with the mean spherical
approximation (MSA) by solving the Ornstein-Zernike (OZ) equation. Both Yukawa
terms can be used to simulate either a repulsive or an attraction potential
feature. For example, if one Yukawa term is chosen to be an attractive potential
and the other a repulsive potential, the two Yukawa potential can simulate a
potential with both a short-range attraction and long-range repulsion
(SALR).[1,2] When combined with an appropriate form factor $P(q)$, this allows
for inclusion of the interparticle interference effects due to a relatively
complex potential.

The algorithm used in this model was originally proposed and developed by Liu et
al. in 2005 and implemented using MatLab.[1] This algorithm used a theory
proposed by Hoye and Blum in 1977,[4] in which the structure factor is solved by
using Baxter's Q-method with the MSA closure.

The interaction potential $V(r)$ is

.. math::

    \frac{V(r)}{k_BT} = \begin{cases}
    \infty & r < 1 \\
    K_1 \frac{e^{-Z_1(r-1)}}{r} + K_2 \frac{e^{-Z_2(r-1)}}{r} & r \geq 1
    \end{cases}

where $r$ is the normalized distance, $\frac{r_{cc}}{2R}$. ($R$ is the effective
radius and $r_{cc}$ is the distance between the center-of-mass of particles.)
$K_1$ and $K_2$ are positive for repulsion and negative for attraction.
$Z_1$ and $Z_2$ are the inversely proportional to the decay rate of the
potential. They should always be positive.

To reduce numerical issues during fitting, the lower boundary of the absolute
values for $K_1$ and $K_2$ are set to 0.0001. That is, when $|K| < 0.0001$ it is
adjusted 0.0001 while preserving the sign. Similarly, when $Z_1$ or $Z_2$ are
less than 0.01, they are adjusted to 0.01. To make the two Yukawa terms
functionally different, $Z_1$ is adjusted to $Z_2+0.01$ if the difference
between them is less than 0.01. The impact of these changes on the accuracy of
$S(Q)$ was negligible for all tested cases. The algorithm works best if
$Z_1>Z_2$, particularly for large $Z$ such as $Z>20$. When $Z_1$ is less than
$Z_2$, the code will automatically swap $Z_1$ and $Z_2$, and $K_1$ and $K_2$.
$S(Q)$ is set to 1.0 for volume fraction $\phi < 10^{-10}$.

.. note::

    The so-called $\beta(q)$ ("beta") approximation [3] can be used
    to calculate the effective structure factor for polydispersity and
    non-sphericity when a system is not too far away from the mono-disperse
    spherical systems.

In general, any closure used to solve the OZ equation is an approximation. Thus,
the structure factor calculated using a closure is an approximation of the real
structure factor. Its accuracy varies for different parameter ranges. One should
not "blindly" trust the results. In addition, when using the O-Z equation to
obtain the structure factor, it assumes that a system is in a liquid state.
However, when the attraction potential is too strong the system may not be in a
simple liquid state, and the OZ equation may not have a solution, or the
solution may be unphysical.

Accuracy for the MSA closure used here to solve the OZ equation for two Yukawa
potential was evaluated by Broccio, et al.[2]  When the overall attraction is
moderate, this algorithm produces reasonably accurate results. However, when the
net attraction is very strong, it has been shown that the fitting algorithm
tends to overestimate the attraction strength. Care is needed when discussing
quantitative fitting results for those systems.

When the interaction distances $Z_1$ and $Z_2$ are large, the method implemented
here may be influenced by the numerical accuracy of a computer. In general, when
$Z<20$, the code should function well for most cases. When $Z>25$, sometimes the
intermediate results of this code may run into the limit of the largest number
that a computer can handle. Therefore, results may potentially become less
reliable. Hence, the check of $g(r)$ becomes very important in those situations.
So far, no limitation of the value of $K$ has been encountered except that they
cannot be zero.

References
----------

#. Y. Liu, W. R. Chen, S-H Chen, *J. Chem. Phys.*, 122 (2005) 044507

#. M. Broccio, D. Costa, Y. Liu, S-H Chen, *J. Chem. Phys.*, 124 (2006) 084501

#. M. Kotlarchyk and S-H Chen, *J. Chem. Phys.*, 79 (1983) 2461-2469

#. J. S. Hoye, L. Blum, J. Stat. Phys., 16(1977) 399-413

Authorship and Verification
---------------------------

* **Author:** Yun Liu **Date:** January 22, 2024
* **Last Modified by:** Yun Liu **Date:** January 22, 2024
* **Last Reviewed by:** Yun Liu **Date:** January 22, 2024
"""
from numpy import inf

# TODO: pep8 says packages and modules should not use camel case
from sasmodels.TwoYukawa.CalTYSk import K_MIN, Z_MIN, Z_MIN_DIFF, CalTYSk

# If you want a customized version of two_yukawa as a plugin (for example,
# because you want to use the high precision polynomial root solver from mpmath)
# you will need to load the customized TwoYukawa library. The following loads it
# from the plugin directory each time the plugin is loaded. You will need copy
# the following into you plugins directory before modifying:
#     https://github.com/SasView/sasmodels/tree/master/sasmodels/models/two_yukawa.py
#     https://github.com/SasView/sasmodels/tree/master/sasmodels/TwoYukawa
if 0:
    import importlib.util
    import sys
    from pathlib import Path

    # Remove existing TwoYukawa from sys.modules to force a reload
    remove = [modname for modname in sys.modules if modname.startswith('TwoYukawa.') or modname == 'TwoYukawa']
    for modname in remove:
        del sys.modules[modname]

    # Load or reload the TwoYukawa package
    path = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location('TwoYukawa', str(path / 'TwoYukawa' / '__init__.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules['TwoYukawa'] = module

    # Override sasmodels library symbols with the local symbols.
    from TwoYukawa.CalTYSk import K_MIN, Z_MIN, Z_MIN_DIFF, CalTYSk

name = "two_yukawa"
title = "User model for two Yukawa structure factor (S(q))"
description = """"""

category = "structure-factor"
structure_factor = True
single = False  # make sure that it has double digit precision
opencl = False  # related with parallel computing

valid = "k1 != 0 && k2 != 0 && z1 != z2"

parameters = [
#   ["name", "units", default, [lower, upper], "type", "description"],
    ["radius_effective", "Ang", 50.0, [0, inf], '', ''],
    ['volfraction', '', 0.1, [0, inf], '', ''],
    ['k1', '', 6.0, [-inf, inf], '', ''],
    ['k2', '', -2.0, [-inf, inf], '', ''],
    ['z1', '', 10.0, [0, inf], '', ''],
    ['z2', '', 2.0, [0, inf], '', ''],
    ]

def Iq(q, radius_effective, volfraction, k1, k2, z1, z2):
    # Low volume fractions lead to singularities when solving for the polynomial roots.
    if volfraction < 1e-10:
        return q*0 + 1

    # Clip parameters to prevent numerical precision problems in CalTYSk
    if abs(k1) < K_MIN:
        k1 = -K_MIN if k1 < 0 else K_MIN
    if abs(k2) < K_MIN:
        k2 = -K_MIN if k2 < 0 else K_MIN
    z1 = max(z1, Z_MIN)
    z2 = max(z2, Z_MIN)
    if abs(z1-z2) < Z_MIN_DIFF:
        z1 = z2 + Z_MIN_DIFF

    # Better numerical stability if Z1 > Z2.
    if z1 < z2:
        z1, z2 = z2, z1
        k1, k2 = k2, k1
    # print(z1,z2,k1,k2,volfraction,Qmax)

    # The form computed by CalTYSk assumes a particle diameter of 1.0 but we can correct
    # for this by scaling our measured Q by the proposed diameter.
    Qeff = 2*radius_effective*q

    # Returns Sk, num_roots, r, Gr, error, cVar
    Sk, *rest = CalTYSk(z1,z2,k1,k2,volfraction,Qeff,warnFlag=False,debugFlag=False)
    return Sk
Iq.vectorized = True  # Iq accepts an array of q value

# The test results generated with MatLab Code (TYSQ22) (Jan. 22, 2024) for
# the potential V(r) = -6*exp(-10*(r-1))/r + 2.0*exp(-2*(r-1))/r
# Note that in the matlab code, the definition of K1 (or K2) is different
# from python : in matlab K1 > 0 means attraction and K1 < 0 means repulsion,
# but in python k1 > 0 means repulsion and k1 < 0 means attraction.
tests = [
   [{'scale': 1.0, 'background': 0.0, 'radius_effective' : 50.0,
      'volfraction': 0.1, 'k1': -6, 'k2': 2.0, 'z1': 10, 'z2': 2.0},
     # Values from matlab code, where K1=-k1 and K2=-k2
     [0.0012513, 0.0311018], [0.242317237668712, 0.869463023177384]],
]
