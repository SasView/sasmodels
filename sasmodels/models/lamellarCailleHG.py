# Note: model title and parameter table are inserted automatically
r"""
This model provides the scattering intensity, *I(q)* = *P(q)* \* *S(q)*, for a lamellar phase where a random
distribution in solution are assumed. Here a Caille S(Q) is used for the lamellar stacks.

The scattering intensity *I(q)* is

.. image:: img/lamellarCailleHG_139.PNG

The form factor is

.. image:: img/lamellarCailleHG_143.PNG

and the structure factor is

.. image:: img/lamellarCailleHG_140.PNG

where

.. image:: img/lamellarCailleHG_141.PNG

where |delta|\ T = tail length (or *tail_length*), |delta|\ H = head thickness (or *h_thickness*),
|drho|\ H = SLD(headgroup) - SLD(solvent), and |drho|\ T = SLD(tail) - SLD(headgroup).
Here *d* = (repeat) spacing, *K* = smectic bending elasticity, *B* = compression modulus, and N = number of lamellar
plates (*n_plates*).

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
background      |cm^-1|   0.001
sld_head        |Ang^-2|  2e-06
scale           None      1
sld_solvent     |Ang^-2|  6e-06
deltaH          |Ang|     2
deltaT          |Ang|     10
sld_tail        |Ang^-2|  0
n_plates        None      30
spacing         |Ang|     40
caille          |Ang^-2|  0.001
==============  ========  =============

.. image:: img/lamellarCailleHG_142.jpg

*Figure. 1D plot using the default values (w/6000 data point).*

Our model uses the form factor calculations implemented in a c-library provided by the NIST Center for Neutron Research
(Kline, 2006).

REFERENCE

F Nallet, R Laversanne, and D Roux, J. Phys. II France, 3, (1993) 487-502

also in J. Phys. Chem. B, 105, (2001) 11081-11088
"""
from numpy import pi, inf

name = "lamellarCailleHG"
title = "Random lamellar sheet with Caille structure factor"
description = """\
	[Random lamellar phase with Caille  structure factor]
        randomly oriented stacks of infinite sheets
		with Caille S(Q), having polydisperse spacing.
		layer thickness =(H+T+T+H) = 2(Head+Tail)
		sld = Tail scattering length density
		sld_head = Head scattering length density
		sld_solvent = solvent scattering length density
		background = incoherent background
		scale = scale factor
"""
category = "shape:lamellae"

parameters = [
#   [ "name", "units", default, [lower, upper], "type",
#     "description" ],
    [ "tail_length", "Ang",  10, [0, inf], "volume",
      "Tail thickness" ],
    [ "head_length", "Ang",  2, [0, inf], "volume",
      "head thickness" ],
    [ "Nlayers", "",  30, [0, inf], "",
      "Number of layers" ],
    [ "spacing", "Ang", 40., [0.0,inf], "volume",
      "d-spacing of Caille S(Q)" ],
    [ "Caille_parameter", "", 0.001, [0.0,0.8], "",
      "Caille parameter" ],
    [ "sld", "1e-6/Ang^2", 0.4, [-inf,inf], "",
      "Tail scattering length density" ],
    [ "head_sld", "1e-6/Ang^2", 2.0, [-inf,inf], "",
      "Head scattering length density" ],
    [ "solvent_sld", "1e-6/Ang^2", 6, [-inf,inf], "",
      "Solvent scattering length density" ],
    ]

source = [ "lamellarCailleHG_kernel.c"]

# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return 1.0;
    """

Iqxy = """
    return Iq(sqrt(qx*qx+qy*qy), IQ_PARAMETERS);
    """

# ER defaults to 0.0
# VR defaults to 1.0

demo = dict(
        scale=1, background=0,
		Nlayers=20,spacing=200.,
        Caille_parameter=0.05,
		tail_length=15,head_length=10,
        sld=-1, head_sld=4.0, solvent_sld=6.0,
		tail_length_pd= 0.1, tail_length_pd_n=20,
		head_length_pd= 0.05, head_length_pd_n=30,
		spacing_pd= 0.2, spacing_pd_n=40
         )

oldname = 'LamellarPSHGModel'
oldpars = dict(tail_length='deltaT',head_length='deltaH',Nlayers='n_plates',Caille_parameter='caille', sld='sld_tail', head_sld='sld_head',solvent_sld='sld_solvent')
