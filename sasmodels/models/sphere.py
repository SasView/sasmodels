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

from numpy import pi, inf

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

parameters = [
#   [ "name", "units", default, [lower, upper], "type",
#     "description" ],
    [ "sld", "1e-6/Ang^2", 1, [-inf,inf], "",
      "Layer scattering length density" ],
    [ "solvent_sld", "1e-6/Ang^2", 6, [-inf,inf], "",
      "Solvent scattering length density" ],
    [ "radius", "Ang",  50, [0, inf], "volume",
      "Sphere radius" ],
    ]


# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return 1.333333333333333*M_PI*radius*radius*radius;
    """

Iq = """
    const double qr = q*radius;
    double sn, cn;
    SINCOS(qr, sn, cn);
    const double bes = qr==0.0 ? 1.0 : 3.0*(sn-qr*cn)/(qr*qr*qr);
    const double fq = bes * (sld - solvent_sld) * form_volume(radius);
    return 1.0e-4*fq*fq;
    """


Iqxy = """
    // never called since no orientation or magnetic parameters.
    //return -1.0;
    return Iq(sqrt(qx*qx + qy*qy), sld, solvent_sld, radius);
    """

def ER(radius):
    return radius

# VR defaults to 1.0