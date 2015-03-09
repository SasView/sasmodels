r"""
For information about polarised and magnetic scattering, click here_.

.. _here: polar_mag_help.html

Definition
----------

The 1D scattering intensity is calculated in the following way (Guinier, 1955)

.. math::

    I(Q) = \frac{\text{scale}}{V} \cdot \left[ \
        3V(\Delta\rho) \cdot \frac{\sin(QR) - QR\cos(QR))}{(QR)^3} \
        \right]^2 + \text{background}

where *scale* is a volume fraction, $V$ is the volume of the scatterer,
$R$ is the radius of the sphere, *background* is the background level and
*sld* and *solvent_sld* are the scattering length densities (SLDs) of the
scatterer and the solvent respectively.

Note that if your data is in absolute scale, the *scale* should represent
the volume fraction (which is unitless) if you have a good fit. If not,
it should represent the volume fraction times a factor (by which your data
might need to be rescaled).

The 2D scattering intensity is the same as above, regardless of the
orientation of $\vec q$.

Our model uses the form factor calculations as defined in the IGOR
package provided by the NIST Center for Neutron Research (Kline, 2006).

Validation
----------

Validation of our code was done by comparing the output of the 1D model
to the output of the software provided by the NIST (Kline, 2006).
Figure :num:`figure #sphere-comparison` shows a comparison of the output
of our model and the output of the NIST software.

.. _sphere-comparison:

.. figure:: img/sphere_comparison.jpg

    Comparison of the DANSE scattering intensity for a sphere with the
    output of the NIST SANS analysis software. The parameters were set to:
    *scale* = 1.0, *radius* = 60 |Ang|, *contrast* = 1e-6 |Ang^-2|, and
    *background* = 0.01 |cm^-1|.


Reference
---------

A Guinier and G. Fournet, *Small-Angle Scattering of X-Rays*,
John Wiley and Sons, New York, (1955)

*2013/09/09 and 2014/01/06 - Description reviewed by S King and P Parker.*
"""

import numpy as np
from numpy import pi, inf, sin, cos, sqrt, log

name = "sphere"
title = "Spheres with uniform scattering length density"
description = """\
P(q)=(scale/V)*[3V(sld-solvent_sld)*(sin(qR)-qRcos(qR))
                /(qR)^3]^2 + background
    R: radius of sphere
    V: The volume of the scatter
    sld: the SLD of the sphere
    solvent_sld: the SLD of the solvent
"""
category = "shape:sphere"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 1, [-inf, inf], "",
               "Layer scattering length density"],
              ["solvent_sld", "1e-6/Ang^2", 6, [-inf, inf], "",
               "Solvent scattering length density"],
              ["radius", "Ang", 50, [0, inf], "volume",
               "Sphere radius"],
             ]


def form_volume(radius):
    return 1.333333333333333 * pi * radius ** 3

def Iq(q, sld, solvent_sld, radius):
    #print "q",q
    #print "sld,r",sld,solvent_sld,radius
    qr = q * radius
    sn, cn = sin(qr), cos(qr)
    # FOR VECTORIZED VERSION, UNCOMMENT THE NEXT TWO LINES
    bes = 3 * (sn - qr * cn) / qr ** 3 # may be 0/0 but we fix that next line
    bes[qr == 0] = 1
    # FOR NON VECTORIZED VERSION, UNCOMMENT THE NEXT LINE
    #bes = 3 * (sn-qr*cn)/qr**3 if qr>0 else 1
    fq = bes * (sld - solvent_sld) * form_volume(radius)
    return 1.0e-4 * fq ** 2
# FOR VECTORIZED VERSION, UNCOMMENT THE NEXT LINE
Iq.vectorized = True

def Iqxy(qx, qy, sld, solvent_sld, radius):
    return Iq(sqrt(qx ** 2 + qy ** 2), sld, solvent_sld, radius)
Iqxy.vectorized = True

def sesans(z, sld, solvent_sld, radius):
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
    g[low] = sqrt(1 - dlow2 / 4.) * (1 + dlow2 / 8.) + dlow2 / 2.*(1 - dlow2 / 16.) * log(dlow / (2. + sqrt(4. - dlow2)))
    return g
sesans.vectorized = True

def ER(radius):
    return radius

# VR defaults to 1.0

demo = dict(scale=1, background=0,
            sld=6, solvent_sld=1,
            radius=120,
            radius_pd=.2, radius_pd_n=45)
oldname = "SphereModel"
oldpars = dict(sld='sldSph', solvent_sld='sldSolv', radius='radius')
