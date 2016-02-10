r"""
This model calculates an empirical functional form for SAS data characterized
by two power laws.

Definition
----------

The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = \begin{cases}
    A q^{-m1} + \text{background} & q <= qc \\
    C q^{-m2} + \text{background} & q > qc
    \end{cases}

where $qc$ = the location of the crossover from one slope to the other,
$A$ = the scaling coefficent that sets the overall intensity of the lower Q power law region,
$m1$ = power law exponent at low Q,
$m2$ = power law exponent at high Q

The scaling of the second power law region (coefficent C) is then automatically scaled to
match the first by following formula

.. math::
    C = \frac{A}{qc^{-m1} qc^{-m2}}

.. note::
    Be sure to enter the power law exponents as positive values!

For 2D data the scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


.. figure:: img/two_power_law_1d.jpg

    1D plot using the default values (w/500 data point).

References
----------

None.

"""

from numpy import power
from numpy import sqrt
from numpy import inf

name = "two_power_law"
title = "Two Power Law intensity calculation"
description = """
            I(q) = coef_A*pow(qval,-1.0*power1) + background for q<=qc
            =C*pow(qval,-1.0*power2) + background for q>qc
            where C=coef_A*pow(qc,-1.0*power1)/pow(qc,-1.0*power2).

            coef_A = scaling coefficent
            qc = crossover location [1/A]
            power_1 (=m1) = power law exponent at low Q
            power_2 (=m2) = power law exponent at high Q
            background = Incoherent background [1/cm]
        """
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#            ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["coefficent_1",  "",         1.0, [-inf, inf], "", "coefficent A in low Q region"],
              ["crossover",     "1/Ang",    0.04,[0, inf], "", "crossover location"],
              ["power_1",       "",         1.0, [0, inf], "", "power law exponent at low Q"],
              ["power_2",       "",         4.0, [0, inf], "", "power law exponent at high Q"],
              ]
# pylint: enable=bad-whitespace, line-too-long

def Iq(q,
       coefficent_1=1.0,
       crossover=0.04,
       power_1=1.0,
       power_2=4.0,
       ):

    """
    :param q:                   Input q-value (float or [float, float])
    :param coefA:               Scaling coefficent at low Q
    :param crossover:           Crossover location
    :param power_1:             Exponent of power law function at low Q
    :param power_2:             Exponent of power law function at high Q
    :return:                    Calculated intensity
    """
# pylint: disable=bad-whitespace

    if q<=crossover:
        intensity = coefficent_1*power(q,-1.0*power_1)
    else:
        coefficent_2 = coefficent_1*power(crossover,-1.0*power_1)/power(crossover,-1.0*power_2)
        intensity = coefficent_2*power(q,-1.0*power_2)

    return intensity

Iq.vectorized = False  # Iq accepts an array of q values

def Iqxy(qx, qy, *args):
    """
    :param qx:   Input q_x-value
    :param qy:   Input q_y-value
    :param args: Remaining arguments
    :return:     2D-Intensity
    """

    return Iq(sqrt(qx**2 + qy**2), *args)

Iqxy.vectorized = False  # Iqxy accepts an array of qx, qy values

demo = dict(scale=1, background=0.1,
            coefficent_1=1.0,
            crossover=0.04,
            power_1=1.0,
            power_2=4.0)

oldname = "TwoPowerLawModel"
oldpars = dict(coefficent_1='coef_A',
                crossover='qc',
                power_1='power1',
                power_2='power2',
                background='background')

tests = [

    # Accuracy tests based on content in test/utest_extra_models.py
    [{'coefficent_1':     1.0,
      'crossover':  0.04,
      'power_1':    1.0,
      'power_2':    4.0,
      'background': 0.0,
    }, 0.001, 1000],

    [{'coefficent_1':     1.0,
      'crossover':  0.04,
      'power_1':    1.0,
      'power_2':    4.0,
      'background': 0.0,
    }, 0.150141, 0.125945],

    [{'coeffcent_1':    1.0,
      'crossover':  0.04,
      'power_1':    1.0,
      'power_2':    4.0,
      'background': 0.0,
    }, 0.442528, 0.00166884],
]
