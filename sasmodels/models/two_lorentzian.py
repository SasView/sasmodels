r"""
This model calculates an empirical functional form for SAS data characterized
by two Lorentzian-type functions.

Definition
----------

The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = \frac{A}{1 +(Q\xi_1)^n} + \frac{C}{1 +(Q\xi_2)^m} + \text{B}

where $A$ = Lorentzian scale factor #1, $C$ = Lorentzian scale #2,
$\xi_1$ and $\xi_2$ are the corresponding correlation lengths, and $n$ and
$m$ are the respective power law exponents (set $n = m = 2$ for
Ornstein-Zernicke behaviour).

For 2D data the scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


.. figure:: img/two_lorentzian.jpg

    1D plot using the default values (w/500 data point).

References
----------

None.

"""

from math import sqrt
from numpy import inf, power

name = "two_lorentzian"
title = "Two Lorentzian type peak"
description = """I(q) = scale_1/(1.0 + pow((q*length_1),exponent_1))
             + scale_2/(1.0 + pow((q*length_2),exponent_2) )+ background

             scale_1    = Lorentzian term scaling #1
             length_1   = Lorentzian screening length #1 [A]
             exponent_1 = Lorentzian exponent #1
             scale_2    = Lorentzian term scaling #2
             length_2   = Lorentzian screening length #2 [A]
             exponent_2 = Lorentzian exponent #2
             background = Incoherent background
        """
category = "shape-independent"

#            ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["lorentz_scale_1",  "",     10.0, [-inf, inf], "",
               "First power law scale factor"],
              ["lorentz_length_1", "Ang", 100.0, [-inf, inf], "",
               "First Lorentzian screening length"],
              ["lorentz_exp_1",    "",      3.0, [-inf, inf], "",
               "First exponent of power law"],
              ["lorentz_scale_2",  "",      1.0, [-inf, inf], "",
               "Second scale factor for broad Lorentzian peak"],
              ["lorentz_length_2", "Ang",  10.0, [-inf, inf], "",
               "Second Lorentzian screening length"],
              ["lorentz_exp_2",    "",      2.0, [-inf, inf], "",
               "Second exponent of power law"],
              ]


def Iq(q,
       lorentz_scale_1,
       lorentz_length_1,
       lorentz_exp_1,
       lorentz_scale_2,
       lorentz_length_2,
       lorentz_exp_2):

    """
    :param q:                   Input q-value (float or [float, float])
    :param lorentz_scale_1:     Second scale factor for broad Lorentzian peak
    :param lorentz_length_1:    First Lorentzian screening length
    :param lorentz_exp_1:       Exponent of the second Lorentz function
    :param lorentz_scale_2:     Second scale factor for broad Lorentzian peak
    :param lorentz_length_2:    Second Lorentzian screening length
    :param lorentz_exp_2:       Exponent of the second Lorentz function
    :return:                    Calculated intensity
    """

    intensity  = lorentz_scale_1/(1.0 +
                                  power(q*lorentz_length_1, lorentz_exp_1))
    intensity += lorentz_scale_2/(1.0 +
                                  power(q*lorentz_length_2, lorentz_exp_2))

    return intensity

Iq.vectorized = True  # Iq accepts an array of q values


def Iqxy(qx, qy, *args):
        iq = Iq(sqrt(qx**2 + qy**2), *args)

        return iq

Iqxy.vectorized = True  # Iqxy accepts an array of qx, qy values


demo = dict(scale=1, background=0.1,
            lorentz_scale_1=10, lorentz_length_1=100.0, lorentz_exp_1=3.0,
            lorentz_scale_2=1,  lorentz_length_2=10,    lorentz_exp_2=2.0)

oldname = "TwoLorentzianModel"
oldpars = dict(background='background',
               lorentz_scale_1='scale_1',   lorentz_scale_2='scale_2',
               lorentz_length_1='length_1', lorentz_length_2='length_2',
               lorentz_exp_1='exponent_1',  lorentz_exp_2='exponent_2')

tests = [
         # Accuracy tests based on content in test/utest_extra_models.py
         [{'lorentz_scale_1':   10.0,
           'lorentz_length_1': 100.0,
           'lorentz_exp_1':      3.0,
           'lorentz_scale_2':    1.0,
           'lorentz_length_2':  10.0,
           'lorentz_exp_2':      2.0,
           'background':         0.1,
           }, 0.001, 11.08991],

         [{'lorentz_scale_1':   10.0,
           'lorentz_length_1': 100.0,
           'lorentz_exp_1':      3.0,
           'lorentz_scale_2':    1.0,
           'lorentz_length_2':  10.0,
           'lorentz_exp_2':      2.0,
           'background':         0.1,
           }, 0.150141, 0.410245],

         [{'lorentz_scale_1':   10.0,
           'lorentz_length_1': 100.0,
           'lorentz_exp_1':      3.0,
           'lorentz_scale_2':    1.0,
           'lorentz_length_2':  10.0,
           'lorentz_exp_2':      2.0,
           'background':         0.1,
           }, 0.442528, 0.148699],

         # Additional tests with larger range of parameters
         [{'lorentz_scale_1':   10.0,
           'lorentz_length_1': 100.0,
           'lorentz_exp_1':      3.0,
           'lorentz_scale_2':    1.0,
           'lorentz_length_2':  10.0,
           'lorentz_exp_2':      2.0,
           }, 0.000332070182643, 10.9996228107],

         [{'lorentz_scale_1':  0.0,
           'lorentz_length_1': 0.0,
           'lorentz_exp_1':    0.0,
           'lorentz_scale_2':  0.0,
           'lorentz_length_2': 0.0,
           'lorentz_exp_2':    0.0,
           'background':     100.0
           }, 5.0, 100.0],

         [{'lorentz_scale_1': 200.0,
           'lorentz_length_1': 10.0,
           'lorentz_exp_1':     0.1,
           'lorentz_scale_2':   0.1,
           'lorentz_length_2':  5.0,
           'lorentz_exp_2':     2.0
           }, 20000., 45.5659201896],
         ]
