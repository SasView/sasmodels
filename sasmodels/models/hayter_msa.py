# Note: model title and parameter table are inserted automatically
r"""
This calculates the structure factor (the Fourier transform of the pair
correlation function $g(r)$) for a system of charged, spheroidal objects
in a dielectric medium. When combined with an appropriate form factor
(such as sphere, core+shell, ellipsoid, etc), this allows for inclusion
of the interparticle interference effects due to screened coulomb repulsion
between charged particles.

**This routine only works for charged particles**. If the charge is set to
zero the routine may self-destruct! For non-charged particles use a hard
sphere potential.

The salt concentration is used to compute the ionic strength of the solution
which in turn is used to compute the Debye screening length. At present
there is no provision for entering the ionic strength directly nor for use
of any multivalent salts, though it should be possible to simulate the effect
of this by increasing the salt concentration. The counterions are also assumed to 
be monovalent.

In sasview the effective radius may be calculated from the parameters
used in the form factor $P(q)$ that this $S(q)$ is combined with.

The computation uses a Taylor series expansion at very small rescaled $qR$, to 
avoid some serious rounding error issues, this may result in a minor artefact 
in the transition region under some circumstances.

For 2D data, the scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


References
----------

J B Hayter and J Penfold, *Molecular Physics*, 42 (1981) 109-118

J P Hansen and J B Hayter, *Molecular Physics*, 46 (1982) 651-656
"""
from numpy import inf

category = "structure-factor"
structure_factor = True
single = False  # double precision only!

#  dp[0] = 2.0*radius_effective();
#  dp[1] = fabs(charge());
#  dp[2] = volfraction();
#  dp[3] = temperature();
#  dp[4] = saltconc();
#  dp[5] = dielectconst();




name = "hayter_msa"
title = "Hayter-Penfold rescaled MSA, charged sphere, interparticle S(Q) structure factor"
description = """\
    [Hayter-Penfold RMSA charged sphere interparticle S(Q) structure factor]
        Interparticle structure factor S(Q)for a charged hard spheres.
        Routine takes absolute value of charge, use HardSphere if charge
        goes to zero.
        In sasview the effective radius and volume fraction may be calculated 
        from the parameters used in P(Q).
"""


# pylint: disable=bad-whitespace, line-too-long
#             [ "name", "units", default, [lower, upper], "type", "description" ],
parameters = [
    ["radius_effective", "Ang", 20.75,   [0, inf],    "volume", "effective radius of charged sphere"],
    ["charge",        "e",   19.0,    [0, inf],    "", "charge on sphere (in electrons)"],
    ["volfraction",   "None",     0.0192, [0, 0.74],   "", "volume fraction of spheres"],
    ["temperature",   "K",  318.16,   [0, inf],    "", "temperature, in Kelvin, for Debye length calculation"],
    ["saltconc",      "M",    0.0,    [-inf, inf], "", "conc of salt, moles/litre, 1:1 electolyte, for Debye length"],
    ["dielectconst",  "None",    71.08,   [-inf, inf], "", "dielectric constant (relative permittivity) of solvent, default water, for Debye length"]
    ]
# pylint: enable=bad-whitespace, line-too-long

source = ["hayter_msa_kernel.c"]
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
#oldpars = dict(effect_radius="radius_effective",effect_radius_pd="radius_effective_pd",effect_radius_pd_n="radius_effective_pd_n")
oldpars = dict(radius_effective="effect_radius",radius_effective_pd="effect_radius_pd",radius_effective_pd_n="effect_radius_pd_n")
#oldpars = dict( )
# default parameter set,  use  compare.sh -midQ -linear
# note the calculation varies in different limiting cases so a wide range of
# parameters will be required for a thorough test!
# odd that the default st has saltconc zero
demo = dict(radius_effective=20.75,
            charge=19.0,
            volfraction=0.0192,
            temperature=318.16,
            saltconc=0.05,
            dielectconst=71.08,
            radius_effective_pd=0.1,
            radius_effective_pd_n=40)
#
# attempt to use same values as old sasview unit test at Q=.001 was 0.0712928, 
# then add lots new ones assuming values from new model are OK, need some low Q values to test the small Q Taylor expansion
tests = [
    [{'scale': 1.0,
      'background': 0.0,
      'radius_effective': 20.75,
      'charge': 19.0,
      'volfraction': 0.0192,
      'temperature': 298.0,
      'saltconc': 0,
      'dielectconst': 78.0,
      'radius_effective_pd': 0},
     [0.00001,0.0010,0.01,0.075], [0.0711646,0.0712928,0.0847006,1.07150]],
    [{'scale': 1.0,
      'background': 0.0,
      'radius_effective': 20.75,
      'charge': 19.0,
      'volfraction': 0.0192,
      'temperature': 298.0,
      'saltconc': 0.05,
      'dielectconst': 78.0,
      'radius_effective_pd': 0.1,
      'radius_effective_pd_n': 40},
     [0.00001,0.0010,0.01,0.075], [0.450272,0.450420,0.465116,1.039625]]
    ]
# ADDED by:  RKH  ON: 16Mar2016 converted from sasview, new Taylor expansion at smallest rescaled Q
