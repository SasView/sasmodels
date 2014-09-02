r"""

"""

from numpy import pi, inf

name = "sphere"
title = "Spheres with uniform SLD"
description = """\
	[Sphere Form Factor]
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
    return REAL(1.333333333333333)*M_PI*radius*radius*radius;
    """

Iq = """
    real sn, cn;
    const real qr = q*radius;
    SINCOS(qr, sn, cn);
    const real bes = (qr==REAL(0.0) ? REAL(1.0) : (REAL(3.0)*(sn-qr*cn))/(qr*qr*qr));
    const real f = bes * (sld - solvent_sld) * form_volume(radius);
    return REAL(1.0e-4)*f*f;
    """


Iqxy = """
    // never called since no orientation or magnetic parameters.
    return REAL(-1.0);
    """

def ER(radius):
    return radius

# VR defaults to 1.0