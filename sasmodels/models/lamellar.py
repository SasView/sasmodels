r"""
Polydispersity in the bilayer thickness can be applied from the GUI.

Definition
----------

The scattering intensity *I(Q)* is

.. math::

    I(Q) = 2\pi{P(Q) \over \delta Q^2}


The form factor is

.. math::

    P(Q) = {2\Delta\rho^2 \over Q^2}(1-cos(Q\delta))


where |delta| = bilayer thickness.

The 2D scattering intensity is calculated in the same way as 1D, where the $Q$ vector is defined as

.. math::

    Q = \sqrt{Q_x^2 + Q_y^2}



Our model uses the form factor calculations implemented in a c-library provided by the NIST Center for Neutron Research
(Kline, 2006).

.. figure:: img/lamellar_1d.jpg

    1D plot using the default values (w/1000 data point).


Reference
---------

F Nallet, R Laversanne, and D Roux, J. Phys. II France, 3, (1993) 487-502

also in J. Phys. Chem. B, 105, (2001) 11081-11088


"""

from numpy import pi, inf

name = "lamellar"
title = "Lyotropic lamellar phase with uniform SLD and random distribution"
description = """\
	[Dilute Lamellar Form Factor](from a lyotropic lamellar phase)
		I(q)= 2*pi*P(q)/(delta *q^(2)), where
		P(q)=2*(contrast/q)^(2)*(1-cos(q*delta))^(2))
		thickness = layer thickness
		sld = layer scattering length density
		sld_solvent = solvent scattering length density
		background = incoherent background
		scale = scale factor
"""

parameters = [
#   [ "name", "units", default, [lower, upper], "type",
#     "description" ],
    [ "sld", "1e-6/Ang^2", 1, [-inf,inf], "",
      "Layer scattering length density" ],
    [ "solvent_sld", "1e-6/Ang^2", 6, [-inf,inf], "",
      "Solvent scattering length density" ],
    [ "thickness", "Ang",  50, [0, inf], "volume",
      "Bilayer thickness" ],
    ]


# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return 1.0;
    """

Iq = """
    const double sub = sld - solvent_sld;
    const double qsq = q*q;
    return 4.0e-4*M_PI*sub*sub/qsq * (1.0-cos(q*thickness))
        / (thickness*qsq);
    """

Iqxy = """
    // never called since no orientation or magnetic parameters.
    return -1.0;
    """

# ER defaults to 0.0
# VR defaults to 1.0