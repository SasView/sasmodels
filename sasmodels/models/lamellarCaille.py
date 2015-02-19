# Note: model title and parameter table are inserted automatically
r"""
This model provides the scattering intensity, *I(q)* = *P(q)* \* *S(q)*, for a lamellar phase where a random
distribution in solution are assumed. Here a Caille S(Q) is used for the lamellar stacks.

The scattering intensity *I(q)* is

.. image:: img/lamellarCaille_139.PNG

The form factor is

.. image:: img/lamellarCaille_134.PNG

and the structure factor is

.. image:: img/lamellarCaille_.PNG

where

.. image:: img/lamellarCaille_.PNG

Here *d* = (repeat) spacing, |delta| = bilayer thickness, the contrast |drho| = SLD(headgroup) - SLD(solvent),
K = smectic bending elasticity, B = compression modulus, and N = number of lamellar plates (*n_plates*).

NB: **When the Caille parameter is greater than approximately 0.8 to 1.0, the assumptions of the model are incorrect.**
And due to a complication of the model function, users are responsible for making sure that all the assumptions are
handled accurately (see the original reference below for more details).

Non-integer numbers of stacks are calculated as a linear combination of results for the next lower and higher values.

The 2D scattering intensity is calculated in the same way as 1D, where the *q* vector is defined as

.. math::

    Q = \sqrt{Q_x^2 + Q_y^2}

The returned value is in units of |cm^-1|, on absolute scale.

==============  ========  =============
Parameter name  Units     Default value
==============  ========  =============
background      |cm^-1|   0.0
contrast        |Ang^-2|  5e-06
scale           None      1
delta           |Ang|     30
n_plates        None      20
spacing         |Ang|     400
caille          |Ang^-2|  0.1
==============  ========  =============

.. image:: img/lamellarPS_142.jpg

*Figure. 1D plot using the default values (w/6000 data point).*

Our model uses the form factor calculations implemented in a c-library provided by the NIST Center for Neutron Research
(Kline, 2006).

REFERENCE

F Nallet, R Laversanne, and D Roux, J. Phys. II France, 3, (1993) 487-502

also in J. Phys. Chem. B, 105, (2001) 11081-11088
"""
from numpy import pi, inf

name = "lamellarPS"
title = "Random lamellar sheet with Caille structure factor"
description = """\
	[Random lamellar phase with Caille  structure factor]
        randomly oriented stacks of infinite sheets
		with Caille S(Q), having polydisperse spacing.
    	sld = sheet scattering length density
		sld_solvent = solvent scattering length density
		background = incoherent background
		scale = scale factor
"""

parameters = [
#   [ "name", "units", default, [lower, upper], "type",
#     "description" ],
    [ "thickness", "Ang",  30.0, [0, inf], "volume",
      "sheet thickness" ],
    [ "Nlayers", "",  20, [0, inf], "",
      "Number of layers" ],
    [ "spacing", "Ang", 400., [0.0,inf], "volume",
      "d-spacing of Caille S(Q)" ],
    [ "Caille_parameter", "Ang^-2", 0.1, [0.0,0.8], "",
      "Caille parameter" ],
    [ "sld", "1e-6/Ang^2", 6.3, [-inf,inf], "",
      "layer scattering length density" ],
    [ "solvent_sld", "1e-6/Ang^2", 1.0, [-inf,inf], "",
      "Solvent scattering length density" ],
    ]

source = [ "lamellarCaille_kernel.c"]

# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return 1.0;
    """

Iqxy = """
    // never called since no orientation or magnetic parameters.
    return -1.0;
    """

# ER defaults to 0.0
# VR defaults to 1.0

demo = dict(
        scale=1, background=0,
		thickness=67.,Nlayers=3.75,spacing=200.,
        Caille_parameter=0.268,sld=1.0, solvent_sld=6.34,
		thickness_pd= 0.1, thickness_pd_n=100,
 		spacing_pd= 0.05, spacing_pd_n=40
         )

oldname = 'LamellarPSModel'
oldpars = dict(thickness='delta',Nlayers='N_plates',Caille_parameter='caille', sld='sld_bi',solvent_sld='sld_sol')
