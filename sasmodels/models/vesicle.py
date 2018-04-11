r"""
Definition
----------

The 1D scattering intensity is calculated in the following way (Guinier, 1955)

.. math::

    P(q) = \frac{\phi}{V_\text{shell}} \left[
           \frac{3V_{\text{core}}({\rho_{\text{solvent}}
           - \rho_{\text{shell}})j_1(qR_{\text{core}})}}{qR_{\text{core}}}
           + \frac{3V_{\text{tot}}(\rho_{\text{shell}}
           - \rho_{\text{solvent}}) j_1(qR_{\text{tot}})}{qR_{\text{tot}}}
           \right]^2 + \text{background}


where $\phi$ is the volume fraction of shell material, $V_{shell}$ is the volume
of the shell, $V_{\text{cor}}$ is the volume of the core, $V_{\text{tot}}$ is
the total volume, $R_{\text{core}}$ is the radius of the core, $R_{\text{tot}}$
is the outer radius of the shell, $\rho_{\text{solvent}}$ is the scattering
length density of the solvent (which is the same as for the core in this case),
$\rho_{\text{scale}}$ is the scattering length density of the shell, background
is a flat background level (due for example to incoherent scattering in the
case of neutrons), and $j_1$ is the spherical bessel function
$j_1 = (\sin(x) - x \cos(x))/ x^2$.

The functional form is identical to a "typical" core-shell structure, except
that the scattering is normalized by the volume that is contributing to the
scattering, namely the volume of the shell alone, the scattering length density
of the core is fixed the same as that of the solvent, the scale factor when the
data are on an absolute scale is equivalent to the volume fraction of material
in the shell rather than the entire core+shell sphere, and the parameterization
is done in terms of the core radius = $R_{\text{core}}$ and the shell
thickness = $R_{\text{tot}} - R_{\text{core}}$.

.. figure:: img/vesicle_geometry.jpg

    Vesicle geometry.

The 2D scattering intensity is the same as *P(q)* above, regardless of the
orientation of the *q* vector which is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


NB: The outer most radius (= *radius* + *thickness*) is used as the effective
radius for *S(Q)* when *P(Q)* \* *S(Q)* is applied.


References
----------

A Guinier and G. Fournet, *Small-Angle Scattering of X-Rays*, John Wiley and
Sons, New York, (1955)

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Paul Butler **Date:** March 20, 2016
* **Last Reviewed by:** Paul Butler **Date:** March 20, 2016
"""

import numpy as np
from numpy import pi, inf

name = "vesicle"
title = "This model provides the form factor, *P(q)*, for an unilamellar \
    vesicle. This is model is effectively identical to the hollow sphere \
    reparameterized to be more intuitive for a vesicle and normalizing the \
    form factor by the volume of the shell."
description = """
    Model parameters:
        radius : the core radius of the vesicle
        thickness: the shell thickness
        sld: the shell SLD
        sld_solvent: the solvent (and core) SLD
        background: incoherent background
        volfraction: shell volume fraction
        scale : scale factor = 1 if on absolute scale"""
category = "shape:sphere"

#             [ "name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld", "1e-6/Ang^2", 0.5, [-inf, inf], "sld",
               "vesicle shell scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 6.36, [-inf, inf], "sld",
               "solvent scattering length density"],
              ["volfraction", "", 0.05, [0, 1.0], "",
               "volume fraction of shell"],
              ["radius", "Ang", 100, [0, inf], "volume",
               "vesicle core radius"],
              ["thickness", "Ang", 30, [0, inf], "volume",
               "vesicle shell thickness"],
             ]

source = ["lib/sas_3j1x_x.c", "vesicle.c"]

def ER(radius, thickness):
    '''
    returns the effective radius used in the S*P calculation

    :param radius: core radius
    :param thickness: shell thickness
    '''
    return radius + thickness

def VR(radius, thickness):
    '''
    returns the volumes of the shell and of the whole sphere including the
    core plus shell - is used to normalize when including polydispersity.

    :param radius: core radius
    :param thickness: shell thickness
    :return whole: volume of core and shell
    :return whole-core: volume of the shell
    '''

    whole = 4./3. * pi * (radius + thickness)**3
    core = 4./3. * pi * radius**3
    return whole, whole - core

def random():
    total_radius = 10**np.random.uniform(1.3, 5)
    radius = total_radius * np.random.uniform(0, 1)
    thickness = total_radius - radius
    volfraction = 10**np.random.uniform(-3, -1)
    pars = dict(
        #background=0,
        scale=1,  # volfraction is part of the model, so scale=1
        radius=radius,
        thickness=thickness,
        volfraction=volfraction,
    )
    return pars

# parameters for demo
demo = dict(sld=0.5, sld_solvent=6.36,
            volfraction=0.05,
            radius=100, thickness=30,
            radius_pd=.2, radius_pd_n=10,
            thickness_pd=.2, thickness_pd_n=10)

# NOTE: test results taken from values returned by SasView 3.1.2, with
# 0.001 added for a non-zero default background.
tests = [[{}, 0.0005, 859.916526646],
         [{}, 0.100600200401, 1.77063682331],
         [{}, 0.5, 0.00355351388906],
         [{}, 'ER', 130.],
         [{}, 'VR', 0.54483386436],
        ]
