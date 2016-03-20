r"""
Definition
----------

The large and small spheres have their own SLD, as well as the solvent. The
surface coverage term is a fractional coverage (maximum of approximately 0.9
for hexagonally-packed spheres on a surface). Since not all of the small
spheres are necessarily attached to the surface, the excess free (small)
spheres scattering is also included in the calculation. The function calculate
follows equations (8)-(12) of the reference below, and the equations are not
reproduced here.

No inter-particle scattering is included in this model.


.. figure:: img/raspberry_geometry.jpg

    Schematic of the raspberry model
    
where *Ro* is the radius of the large sphere, *Rp* the radius of the smaller 
spheres on the surface and |delta| = the fractional penetration depth.

For 2D data: The 2D scattering intensity is calculated in the same way as 1D,
where the *q* vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


References
----------

K Larson-Smith, A Jackson, and D C Pozzo, *Small angle scattering model for Pickering emulsions and raspberry*
*particles*, *Journal of Colloid and Interface Science*, 343(1) (2010) 36-41

**Author:** Andrew jackson **on:** 2008

**Modified by:** Paul Butler **on:** March 18, 2016

**Reviewed by:** Paul Butler **on:** March 18, 2016
"""

from numpy import pi, inf

name = "raspberry"
title = "Calculates the form factor, *P(q)*, for a 'Raspberry-like' structure \
where there are smaller spheres at the surface of a larger sphere, such as the \
structure of a Pickering emulsion."
description = """
                RaspBerryModel:
                volf_Lsph = volume fraction large spheres
                radius_Lsph = radius large sphere (A)
                sld_Lsph = sld large sphere (A-2)
                volf_Ssph = volume fraction small spheres
                radius_Ssph = radius small sphere (A)
                surfrac_Ssph = fraction of small spheres at surface
                sld_Ssph = sld small sphere
                delta_Ssph = small sphere penetration (A) 
                sld_solv   = sld solvent
                background = background (cm-1)
            Ref: J. coll. inter. sci. (2010) vol. 343 (1) pp. 36-41."""
category = "shape:sphere"

#             [ "name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld_lg", "1e-6/Ang^2", -0.4, [-inf, inf], "",
               "large particle scattering length density"],
              ["sld_sm", "1e-6/Ang^2", 3.5, [-inf, inf], "",
               "small particle scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 6.36, [-inf, inf], "",
               "solvent scattering length density"],
              ["volfraction_lg", "", 0.05, [-inf, inf], "",
               "volume fraction of large spheres"],
              ["volfraction_sm", "", 0.005, [-inf, inf], "",
               "volume fraction of small spheres"],
              ["surf_fraction", "", 0.4, [-inf, inf], "",
               "fraction of small spheres at surface"],
              ["radius_lg", "Ang", 5000, [0, inf], "volume",
               "radius of large spheres"],
              ["radius_sm", "Ang", 100, [0, inf], "",
               "radius of small spheres"],
              ["penetration", "Ang", 0.0, [0, inf], "",
               "penetration depth of small spheres into large sphere"],
             ]

source = ["lib/sph_j1c.c", "raspberry.c"]

# parameters for demo
demo = dict(scale=1, background=0.001,
            sld_lg=-0.4, sld_sm=3.5, sld_solvent=6.36,
            volfraction_lg=0.05, volfraction_sm=0.005, surf_fraction=0.4,
            radius_lg=5000, radius_sm=100, penetration=0.0,
            radius_lg_pd=.2, radius_lg_pd_n=10)

# For testing against the old sasview models, include the converted parameter
# names and the target sasview model name.
oldname = 'RaspBerryModel'
oldpars = dict(sld_lg='sld_Lsph', sld_sm='sld_Ssph', sld_solvent='sld_solv',
               volfraction_lg='volf_Lsph', volfraction_sm='volf_Ssph',
               surf_fraction='surfrac_Ssph',
               radius_lg='radius_Lsph', radius_sm='radius_Ssph',
               penetration='delta_Ssph')


# NOTE: test results taken from values returned by SasView 3.1.2, with
# 0.001 added for a non-zero default background.
tests = [[{}, 0.0412755102041, 0.286669115234],
         [{}, 0.5, 0.00103818393658],
        ]
