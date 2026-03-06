r"""
Example of a pure python model with Fq() defined. This should support the β approximation
for polydisperse structure factor calculations.

Test using::

    python -m sasmodels.compare explore/spherepy.py@hardsphere,sphere@hardsphere -pars -midq structure_factor_mode=1 radius_pd=0.2
"""

import numpy as np
from numpy import cos, inf, log, pi, sin, sqrt

name = "sphere (python)"
title = "PAK testing ideas for Spheres with uniform scattering length density"
description = """\
P(q)=(scale/V)*[3V(sld-sld_solvent)*(sin(qr)-qr cos(qr))
                /(qr)^3]^2 + background
    r: radius of sphere
    V: The volume of the scatter
    sld: the SLD of the sphere
    sld_solvent: the SLD of the solvent
"""
category = "shape:sphere"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 1, [-inf, inf], "",
               "Layer scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 6, [-inf, inf], "",
               "Solvent scattering length density"],
              ["radius", "Ang", 50, [0, inf], "volume",
               "Sphere radius"],
             ]


def form_volume(radius):
    """Calculate volume for sphere"""
    return 1.333333333333333 * pi * radius ** 3

def radius_effective(mode, radius):
    """Calculate R_eff for sphere"""
    return radius if mode else 0.

vectorized = True  # For testing: toggle between vectorized and non-vectorized versions
# have_Fq = False  # For testing: uncomment to force Iq() rather than Fq()

def Fq(q, sld, sld_solvent, radius):
    """Calculate F(q), F^2(q) for sphere"""
    #print "q",q
    #print "sld,r",sld,sld_solvent,radius
    qr = q * radius
    sn, cn = sin(qr), cos(qr)
    if vectorized:
        with np.errstate(all='ignore'):
            bes = 3 * (sn - qr * cn) / qr ** 3 # may be 0/0 but we fix that next line
            bes[qr == 0] = 1
    else:
        bes = 3 * (sn-qr*cn)/qr**3 if qr != 0 else 1
    fq = bes * (1e-2 * (sld - sld_solvent) * form_volume(radius))
    return fq, fq**2
Fq.vectorized = vectorized  # Fq accepts an array of q value

def Iq(q, sld, sld_solvent, radius):
    """Calculate I(q) for sphere"""
    return Fq(q, sld, sld_solvent, radius)[1]
Iq.vectorized = vectorized  # Iq accepts an array of q value

def sesans(z, sld, sld_solvent, radius):
    """
    Calculate SESANS-correlation function for a solid sphere.

    Wim Bouwman after formulae Timofei Kruglov J.Appl.Cryst. 2003 article
    """
    d = z / radius
    g = np.zeros_like(z)
    g[d == 0] = 1.
    low = ((d > 0) & (d < 2))
    dlow = d[low]
    dlow2 = dlow ** 2
    g[low] = (sqrt(1 - dlow2/4.) * (1 + dlow2/8.)
              + dlow2/2.*(1 - dlow2/16.) * log(dlow / (2. + sqrt(4. - dlow2))))
    return g
sesans.vectorized = True  # sesans accepts an array of z values
