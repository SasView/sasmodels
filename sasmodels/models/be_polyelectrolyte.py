r"""
Definition
----------
This model calculates the structure factor of a polyelectrolyte solution with
the RPA expression derived by Borue and Erukhimovich\ [#Borue]_.  Note however
that the fitting procedure here does not follow the notation in that reference
as 's' and 't' are **not** decoupled. Instead the scattering intensity $I(q)$
is calculated as

.. math::

    I(q) = K\frac{q^2+k^2}{4\pi L_b\alpha ^2}
    \frac{1}{1+r_{0}^2(q^2+k^2)(q^2-12hC_a/b^2)} + background

    k^2 = 4\pi L_b(2C_s + \alpha C_a)

    r_{0}^2 = \frac{1}{\alpha \sqrt{C_a} \left( b/\sqrt{48\pi L_b}\right)}

where

$K$ is the contrast factor for the polymer which is defined differently than in
other models and is given in barns where $1 barn = 10^{-24} cm^2$.  $K$ is
defined as:

.. math::

    K = a^2

    a = b_p - (v_p/v_s) b_s

where $b_p$ and $b_s$ are sum of the scattering lengths of the atoms
constituting the monomer of the polymer and the sum of the scattering lengths
of the atoms constituting the solvent molecules respectively, and $v_p$ and
$v_s$ are the partial molar volume of the polymer and the solvent respectively

$L_b$ is the Bjerrum length(|Ang|) - **Note:** This parameter needs to be
kept constant for a given solvent and temperature!

$h$ is the virial parameter (|Ang^3|/mol) - **Note:** See [#Borue]_ for the
correct interpretation of this parameter.  It incorporates second and third
virial coefficients and can be Negative.

$b$ is the monomer length(|Ang|), $C_s$ is the concentration of monovalent
salt(mol/L), $\alpha$ is the ionization degree (ionization degree : ratio of
charged monomers  to total number of monomers), $C_a$ is the polymer molar
concentration(mol/L), and $background$ is the incoherent background.

For 2D data the scattering intensity is calculated in the same way as 1D,
where the $\vec q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}

References
----------

.. [#Borue] V Y Borue, I Y Erukhimovich, *Macromolecules*, 21 (1988) 3240
.. [#] J F Joanny, L Leibler, *Journal de Physique*, 51 (1990) 545
.. [#] A Moussaid, F Schosseler, J P Munch, S Candau, *J. Journal de Physique
   II France*, 3 (1993) 573
.. [#] E Raphael, J F Joanny, *Europhysics Letters*, 11 (1990) 179

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Paul Kienzle **Date:** July 24, 2016
* **Last Reviewed by:** Paul Butler and Richard Heenan **Date:** October 07, 2016
"""

import numpy as np
from numpy import inf, pi, sqrt

name = "be_polyelectrolyte"
title = "Polyelectrolyte with the RPA expression derived by Borue and Erukhimovich"
description = """
            Evaluate
            F(x) = K 1/(4 pi Lb (alpha)^(2)) (q^(2)+k2)/(1+(r02)^(2))
                 (q^(2)+k2) (q^(2)-(12 h C/b^(2)))

            has 3 internal parameters :
                   The inverse Debye Length: K2 = 4 pi Lb (2 Cs+alpha C)
                   r02 =1/alpha/Ca^(0.5) (B/(48 pi Lb)^(0.5))
                   Ca = 6.022136e-4 C
            """
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["contrast_factor",       "barns",   10.0,  [-inf, inf], "", "Contrast factor of the polymer"],
    ["bjerrum_length",        "Ang",      7.1,  [0, inf],    "", "Bjerrum length"],
    ["virial_param",          "Ang^3/mol", 12.0,  [-inf, inf], "", "Virial parameter"],
    ["monomer_length",        "Ang",     10.0,  [0, inf],    "", "Monomer length"],
    ["salt_concentration",    "mol/L",    0.0,  [-inf, inf], "", "Concentration of monovalent salt"],
    ["ionization_degree",     "",         0.05, [0, inf],    "", "Degree of ionization"],
    ["polymer_concentration", "mol/L",    0.7,  [0, inf],    "", "Polymer molar concentration"],
    ]
# pylint: enable=bad-whitespace, line-too-long


def Iq(q,
       contrast_factor=10.0,
       bjerrum_length=7.1,
       virial_param=12.0,
       monomer_length=10.0,
       salt_concentration=0.0,
       ionization_degree=0.05,
       polymer_concentration=0.7):
    """
    :param q:                     Input q-value
    :param contrast_factor:       Contrast factor of the polymer
    :param bjerrum_length:        Bjerrum length
    :param virial_param:          Virial parameter
    :param monomer_length:        Monomer length
    :param salt_concentration:    Concentration of monovalent salt
    :param ionization_degree:     Degree of ionization
    :param polymer_concentration: Polymer molar concentration
    :return:                      1-D intensity
    """

    concentration = polymer_concentration * 6.022136e-4

    k_square = 4.0 * pi * bjerrum_length * (2*salt_concentration +
                                            ionization_degree * concentration)

    r0_square = 1.0/ionization_degree/sqrt(concentration) * \
                (monomer_length/sqrt((48.0*pi*bjerrum_length)))

    term1 = contrast_factor/(4.0 * pi * bjerrum_length *
                             ionization_degree**2) * (q**2 + k_square)

    term2 = 1.0 + r0_square**2 * (q**2 + k_square) * \
        (q**2 - (12.0 * virial_param * concentration/(monomer_length**2)))

    return term1/term2

Iq.vectorized = True  # Iq accepts an array of q values

def random():
    # TODO: review random be_polyelectrolyte model generation
    pars = dict(
        scale=10000, #background=0,
        #polymer_concentration=0.7,
        polymer_concentration=np.random.beta(5, 3), # around 70%
        #salt_concentration=0.0,
        # keep salt concentration extremely low
        # and use explicit molar to match polymer concentration
        salt_concentration=np.random.beta(1, 100)*6.022136e-4,
        #contrast_factor=10.0,
        contrast_fact=np.random.uniform(1, 100),
        #bjerrum_length=7.1,
        bjerrum_length=np.random.uniform(1, 10),
        #virial_param=12.0,
        virial_param=np.random.uniform(-1000, 30),
        #monomer_length=10.0,
        monomer_length=10.0**(4*np.random.beta(1.5, 3)),
        #ionization_degree=0.05,
        ionization_degree=np.random.beta(1.5, 4),
    )
    return pars

demo = dict(scale=1, background=0.1,
            contrast_factor=10.0,
            bjerrum_length=7.1,
            virial_param=12.0,
            monomer_length=10.0,
            salt_concentration=0.0,
            ionization_degree=0.05,
            polymer_concentration=0.7)

tests = [

    # Accuracy tests based on content in test/utest_other_models.py
    [{'contrast_factor':       10.0,
      'bjerrum_length':         7.1,
      'virial_param':          12.0,
      'monomer_length':        10.0,
      'salt_concentration':     0.0,
      'ionization_degree':      0.05,
      'polymer_concentration':  0.7,
      'background':             0.001,
     }, 0.001, 0.0948379],

    # Additional tests with larger range of parameters
    [{'contrast_factor':       10.0,
      'bjerrum_length':       100.0,
      'virial_param':           3.0,
      'monomer_length':         1.0,
      'salt_concentration':    10.0,
      'ionization_degree':      2.0,
      'polymer_concentration': 10.0,
      'background':             0.0,
     }, 0.1, -3.75693800588],

    [{'contrast_factor':       10.0,
      'bjerrum_length':       100.0,
      'virial_param':           3.0,
      'monomer_length':         1.0,
      'salt_concentration':    10.0,
      'ionization_degree':      2.0,
      'polymer_concentration': 10.0,
      'background':           100.0
     }, 5.0, 100.029142149],

    [{'contrast_factor':     100.0,
      'bjerrum_length':       10.0,
      'virial_param':        180.0,
      'monomer_length':        1.0,
      'salt_concentration':    0.1,
      'ionization_degree':     0.5,
      'polymer_concentration': 0.1,
      'background':             0.0,
     }, 200., 1.80664667511e-06],
    ]
