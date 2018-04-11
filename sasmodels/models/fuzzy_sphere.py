r"""
For information about polarised and magnetic scattering, see
the :ref:`magnetism` documentation.

Definition
----------

The scattering intensity $I(q)$ is calculated as:

.. math::

    I(q) = \frac{\text{scale}}{V}(\Delta \rho)^2 A^2(q) S(q)
           + \text{background}


where the amplitude $A(q)$ is given as the typical sphere scattering convoluted
with a Gaussian to get a gradual drop-off in the scattering length density:

.. math::

    A(q) = \frac{3\left[\sin(qR) - qR \cos(qR)\right]}{(qR)^3}
           \exp\left(\frac{-(\sigma_\text{fuzzy}q)^2}{2}\right)

Here $A(q)^2$ is the form factor, $P(q)$. The scale is equivalent to the
volume fraction of spheres, each of volume, $V$. Contrast $(\Delta \rho)$
is the difference of scattering length densities of the sphere and the
surrounding solvent.

Poly-dispersion in radius and in fuzziness is provided for, though the
fuzziness must be kept much smaller than the sphere radius for meaningful
results.

From the reference:

  The "fuzziness" of the interface is defined by the parameter
  $\sigma_\text{fuzzy}$. The particle radius $R$ represents the radius of the
  particle where the scattering length density profile decreased to 1/2 of the
  core density. $\sigma_\text{fuzzy}$ is the width of the smeared particle
  surface; i.e., the standard deviation from the average height of the fuzzy
  interface. The inner regions of the microgel that display a higher density
  are described by the radial box profile extending to a radius of
  approximately $R_\text{box} \sim R - 2 \sigma$. The profile approaches
  zero as $R_\text{sans} \sim R + 2\sigma$.

For 2D data: The 2D scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math:: q = \sqrt{{q_x}^2 + {q_y}^2}

References
----------

M Stieger, J. S Pedersen, P Lindner, W Richtering, *Langmuir*,
20 (2004) 7283-7292
"""

import numpy as np
from numpy import inf

name = "fuzzy_sphere"
title = "Scattering from spherical particles with a fuzzy surface."
description = """\
scale: scale factor times volume fraction,
or just volume fraction for absolute scale data
radius: radius of the solid sphere
fuzziness = the standard deviation of the fuzzy interfacial
thickness (ie., so-called interfacial roughness)
sld: the SLD of the sphere
solvend_sld: the SLD of the solvent
background: incoherent background
Note: By definition, this function works only when fuzziness << radius.
"""
category = "shape:sphere"

# pylint: disable=bad-whitespace,line-too-long
# ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld",         "1e-6/Ang^2",  1, [-inf, inf], "sld",    "Particle scattering length density"],
              ["sld_solvent", "1e-6/Ang^2",  3, [-inf, inf], "sld",    "Solvent scattering length density"],
              ["radius",      "Ang",        60, [0, inf],    "volume", "Sphere radius"],
              ["fuzziness",   "Ang",        10, [0, inf],    "",       "std deviation of Gaussian convolution for interface (must be << radius)"],
             ]
# pylint: enable=bad-whitespace,line-too-long

source = ["lib/sas_3j1x_x.c"]

# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return M_4PI_3*cube(radius);
    """

Iq = """
    const double qr = q*radius;
    const double bes = sas_3j1x_x(qr);
    const double qf = q*fuzziness;
    const double fq = bes * (sld - sld_solvent) * form_volume(radius) * exp(-0.5*qf*qf);
    return 1.0e-4*fq*fq;
    """

def ER(radius):
    """
    Return radius
    """
    return radius

# VR defaults to 1.0

def random():
    radius = 10**np.random.uniform(1, 4.7)
    fuzziness = 10**np.random.uniform(-2, -0.5)*radius  # 1% to 31% fuzziness
    pars = dict(
        radius=radius,
        fuzziness=fuzziness,
    )
    return pars

demo = dict(scale=1, background=0.001,
            sld=1, sld_solvent=3,
            radius=60,
            fuzziness=10,
            radius_pd=.2, radius_pd_n=45,
            fuzziness_pd=.2, fuzziness_pd_n=0)

tests = [
    # Accuracy tests based on content in test/utest_models_new1_3.py
    #[{'background': 0.001}, 1.0, 0.001],

    [{}, 0.00301005, 359.2315],

    ]
