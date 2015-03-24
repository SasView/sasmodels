# Note: model title and parameter table are inserted automatically
r"""
This model calculates the scattering from a stack of repeating lamellar
structures. The stacks of lamellae (infinite in lateral dimension) are
treated as a paracrystal to account for the repeating spacing. The repeat
distance is further characterized by a Gaussian polydispersity. **This model
can be used for large multilamellar vesicles.**

*2.1.33.1. Definition*

The scattering intensity *I(q)* is calculated as

.. image:: img/image145.jpg

The form factor of the bilayer is approximated as the cross section of an
infinite, planar bilayer of thickness *t*

.. image:: img/image146.jpg

Here, the scale factor is used instead of the mass per area of the
bilayer (*G*). The scale factor is the volume fraction of the material in
the bilayer, *not* the total excluded volume of the paracrystal.
*Z*\ :sub:`N`\ *(q)* describes the interference effects for aggregates
consisting of more than one bilayer. The equations used are (3-5)
from the Bergstrom reference below.

Non-integer numbers of stacks are calculated as a linear combination of
the lower and higher values

.. image:: img/image147.jpg

The 2D scattering intensity is the same as 1D, regardless of the orientation
of the *q* vector which is defined as

.. math::

    Q = \sqrt{Q_x^2 + Q_y^2}

The parameters of the model are *Nlayers* = no. of layers, and
*pd_spacing* = polydispersity of spacing.

==============  ========  =============
Parameter name  Units     Default value
==============  ========  =============
background      |cm^-1|   0
scale           None      1
Nlayers         None      20
pd_spacing      None      0.2
sld_layer       |Ang^-2|  1e-6
sld_solvent     |Ang^-2|  6.34e-6
spacing         |Ang|     250
thickness       |Ang|     33
==============  ========  =============

.. image:: img/image148.jpg

*Figure. 1D plot using the default values above (w/20000 data point).*

Our model uses the form factor calculations implemented in a C library
provided by the NIST Center for Neutron Research (Kline, 2006).

REFERENCE

M Bergstrom, J S Pedersen, P Schurtenberger, S U Egelhaaf,
*J. Phys. Chem. B*, 103 (1999) 9888-9897

"""

from numpy import inf

name = "lamellarPC"
title = "Random lamellar sheet with paracrystal structure factor"
description = """\
    [Random lamellar phase with paracrystal structure factor]
        randomly oriented stacks of infinite sheets
        with paracrytal S(Q), having polydisperse spacing.
        sld = sheet scattering length density
        sld_solvent = solvent scattering length density
        background = incoherent background
        scale = scale factor
"""
category = "shape:lamellae"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["thickness", "Ang", 33.0, [0, inf], "volume",
               "sheet thickness"],
              ["Nlayers", "", 20, [0, inf], "",
               "Number of layers"],
              ["spacing", "Ang", 250., [0.0, inf], "",
               "d-spacing of paracrystal stack"],
              ["spacing_polydisp", "Ang", 0.0, [0.0, inf], "",
               "d-spacing of paracrystal stack"],
              ["sld", "1e-6/Ang^2", 1.0, [-inf, inf], "",
               "layer scattering length density"],
              ["solvent_sld", "1e-6/Ang^2", 6.34, [-inf, inf], "",
               "Solvent scattering length density"],
             ]


source = ["lamellarPC_kernel.c"]

form_volume = """
    return 1.0;
    """

Iqxy = """
    return Iq(sqrt(qx*qx+qy*qy), IQ_PARAMETERS);
    """

# ER defaults to 0.0
# VR defaults to 1.0

demo = dict(scale=1, background=0,
            thickness=33, Nlayers=20, spacing=250, spacing_polydisp=0.2,
            sld=1.0, solvent_sld=6.34,
            thickness_pd=0.2, thickness_pd_n=40)

oldname = 'LamellarPCrystalModel'
oldpars = dict(spacing_polydisp='pd_spacing', sld='sld_layer',
               solvent_sld='sld_solvent')
