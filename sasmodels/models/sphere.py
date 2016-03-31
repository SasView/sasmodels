r"""
For information about polarised and magnetic scattering, see 
the :doc:`magnetic help <../sasgui/perspectives/fitting/mag_help>` documentation.

Definition
----------

The 1D scattering intensity is calculated in the following way (Guinier, 1955)

.. math::

    I(q) = \frac{\text{scale}}{V} \cdot \left[
        3V(\Delta\rho) \cdot \frac{\sin(qr) - qr\cos(qr))}{(qr)^3}
        \right]^2 + \text{background}

where *scale* is a volume fraction, $V$ is the volume of the scatterer,
$r$ is the radius of the sphere, *background* is the background level and
*sld* and *sld_solvent* are the scattering length densities (SLDs) of the
scatterer and the solvent respectively.

Note that if your data is in absolute scale, the *scale* should represent
the volume fraction (which is unitless) if you have a good fit. If not,
it should represent the volume fraction times a factor (by which your data
might need to be rescaled).

The 2D scattering intensity is the same as above, regardless of the
orientation of $\vec q$.

Validation
----------

Validation of our code was done by comparing the output of the 1D model
to the output of the software provided by the NIST (Kline, 2006).


References
----------

A Guinier and G. Fournet, *Small-Angle Scattering of X-Rays*,
John Wiley and Sons, New York, (1955)

*2013/09/09 and 2014/01/06 - Description reviewed by S King and P Parker.*
"""

from numpy import inf

name = "sphere"
title = "Spheres with uniform scattering length density"
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

source = ["lib/sph_j1c.c", "lib/sphere_form.c"]

# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return sphere_volume(radius);
    """

Iq = """
    return sphere_form(q, radius, sld, sld_solvent);
    """

Iqxy = """
    // never called since no orientation or magnetic parameters.
    //return -1.0;
    return Iq(sqrt(qx*qx + qy*qy), sld, sld_solvent, radius);
    """

def ER(radius):
    """
        Return equivalent radius (ER)
    """
    return radius

# VR defaults to 1.0

demo = dict(scale=1, background=0,
            sld=6, sld_solvent=1,
            radius=120,
            radius_pd=.2, radius_pd_n=45)
oldname = "SphereModel"
oldpars = dict(sld='sldSph', sld_solvent='sldSolv', radius='radius')
