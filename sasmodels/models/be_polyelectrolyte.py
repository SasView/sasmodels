r"""
This model calculates the structure factor of a polyelectrolyte solution with
the RPA expression derived by Borue and Erukhimovich.

Definition
----------

The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = K\frac{q^2+k^2}{4\pi L\alpha ^2}
    \frac{1}{1+r_{0}^2(q^2+k^2)(q^2-12hC_a/b^2)} + background

    k^2 = 4\pi L(2C_s + \alpha C_a)

    r_{0}^2 = \frac{1}{\alpha \sqrt{C_a} \left( b/\sqrt{48\pi L_b}\right)}

where $K$ is the contrast factor for the polymer, $L_b$ is the Bjerrum length,
$h$ is the virial parameter, $b$ is the monomer length,
$C_s$ is the concentration of monovalent salt, $\alpha$ is the ionization
degree, $C_a$ is the polymer molar concentration, and $background$ is the
incoherent background.

For 2D data the scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


.. figure:: img/be_polyelectrolyte_1d.jpg

    1D plot using the default values (w/500 data point).

NB: $1 barn = 10^{-24} cm^2$

References
----------

V Y Borue, I Y Erukhimovich, *Macromolecules*, 21 (1988) 3240

J F Joanny, L Leibler, *Journal de Physique*, 51 (1990) 545

A Moussaid, F Schosseler, J P Munch, S Candau, *J. Journal de Physique II France*, 3 (1993) 573

E Raphael, J F Joanny, *Europhysics Letters*, 11 (1990) 179

"""

from numpy import inf, power, pi, sqrt

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

#            ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["contrast_factor",       "barns",   10.0,  [-inf, inf], "", "Contrast factor of the polymer"],
              ["bjerrum_length",        "Ang",      7.1,  [0, inf],    "", "Bjerrum length"],
              ["virial_param",          "1/Ang^2", 12.0,  [-inf, inf], "", "Virial parameter"],
              ["monomer_length",        "Ang",     10.0,  [0, inf],    "", "Monomer length"],
              ["salt_concentration",    "mol/L",    0.0,  [-inf, inf], "", "Concentration of monovalent salt"],
              ["ionization_degree",     "",         0.05, [0, inf],    "", "Degree of ionization"],
              ["polymer_concentration", "mol/L",    0.7,  [0, inf],    "", "Polymer molar concentration"],
              ]


def Iq(q,
       contrast_factor,
       bjerrum_length,
       virial_param,
       monomer_length,
       salt_concentration,
       ionization_degree,
       polymer_concentration):

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


def Iqxy(qx, qy, *args):
        iq = Iq(sqrt(qx**2 + qy**2), *args)

        return iq

Iqxy.vectorized = True  # Iqxy accepts an array of qx, qy values


demo = dict(scale=1, background=0.1,
            contrast_factor=10.0,
            bjerrum_length=7.1,
            virial_param=12.0,
            monomer_length=10.0,
            salt_concentration=0.0,
            ionization_degree=0.05,
            polymer_concentration=0.7)

oldname = "BEPolyelectrolyte"

oldpars = dict(background='background',
               contrast_factor='k',
               bjerrum_length='lb',
               virial_param='h',
               monomer_length='b',
               salt_concentration='cs',
               ionization_degree='alpha',
               polymer_concentration='c')

tests = [[{'contrast_factor':       10.0,
           'bjerrum_length':       100.0,
           'virial_param':           3.0,
           'monomer_length':         1.0,
           'salt_concentration':    10.0,
           'ionization_degree':      2.0,
           'polymer_concentration': 10.0,
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
           }, 200., 1.80664667511e-06],
         ]
