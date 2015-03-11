# Note: model title and parameter table are inserted automatically
r"""
This calculates the structure factor (the Fourier transform of the pair correlation function *g(r)*) for a system of
charged, spheroidal objects in a dielectric medium. When combined with an appropriate form factor (such as sphere,
core+shell, ellipsoid, etc), this allows for inclusion of the interparticle interference effects due to screened coulomb
repulsion between charged particles.

**This routine only works for charged particles**. If the charge is set to zero the routine will self-destruct!
For non-charged particles use a hard sphere potential.

The salt concentration is used to compute the ionic strength of the solution which in turn is used to compute the Debye
screening length. At present there is no provision for entering the ionic strength directly nor for use of any
multivalent salts. The counterions are also assumed to be monovalent.

For 2D data: The 2D scattering intensity is calculated in the same way as 1D, where the *q* vector is defined as

.. math::

    Q = \sqrt{Q_x^2 + Q_y^2}

==============  ========  =============
Parameter name  Units     Default value
==============  ========  =============
effect_radius   |Ang|     20.8
charge          *e*       19
volfraction     None      0.2
temperature     K         318
salt conc       M         0
dielectconst    None      71.1
==============  ========  =============

.. image:: img/HayterMSAsq_227.jpg

*Figure. 1D plot using the default values (in linear scale).*

REFERENCE

J B Hayter and J Penfold, *Molecular Physics*, 42 (1981) 109-118

J P Hansen and J B Hayter, *Molecular Physics*, 46 (1982) 651-656
"""

#  dp[0] = 2.0*effect_radius();
#  dp[1] = fabs(charge());
#  dp[2] = volfraction();
#  dp[3] = temperature();
#  dp[4] = saltconc();
#  dp[5] = dielectconst();

from numpy import inf

source = ["HayterMSAsq_kernel.c"]

name = "HayterMSAsq"
title = "Hayter-Penfold MSA charged sphere interparticle S(Q) structure factor"
description = """\
    [Hayter-Penfold MSA charged sphere interparticle S(Q) structure factor]
        Interparticle structure factor S(Q)for a charged hard spheres.
        Routine takes absolute value of charge, use HardSphere if charge goes to zero.
        In sasview the effective radius will be calculated from the
        parameters used in P(Q).
"""
#             [ "name", "units", default, [lower, upper], "type", "description" ],
parameters = [["effect_radius", "Ang", 20.75, [0, inf], "volume",
               "effective radius of hard sphere"],
              ["charge", "e", 19.0, [0, inf], "",
               "charge on sphere (in electrons)"],
              ["volfraction", "", 0.0192, [0, 0.74], "",
               "volume fraction of spheres"],
              ["temperature", "K", 318.16, [0, inf], "",
               "temperature, in Kelvin, for Debye length calculation"],
              ["saltconc", "M", 0.0, [-inf, inf], "",
               "conc of salt, 1:1 electolyte, for Debye length"],
              ["dielectconst", "", 71.08, [-inf, inf], "",
               "dielectric constant of solvent (default water), for Debye length"],
             ]

# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return 1.0;
    """
Iqxy = """
    // never called since no orientation or magnetic parameters.
    return Iq(sqrt(qx*qx+qy*qy), IQ_PARAMETERS);
    """
# ER defaults to 0.0
# VR defaults to 1.0

oldname = 'HayterMSAStructure'
oldpars = dict()
# default parameter set,  use  compare.sh -midQ -linear
# note the calculation varies in different limiting cases so a wide range of parameters will be required for a thorough test!
# odd that the default st has saltconc zero
demo = dict(effect_radius = 20.75,charge=19.0,volfraction = 0.0192,temperature=318.16,saltconc=0.05,dielectconst=71.08,effect_radius_pd = 0.1,effect_radius_pd_n = 40)

