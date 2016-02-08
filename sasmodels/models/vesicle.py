r"""
Definition
----------

The 1D scattering intensity is calculated in the following way (Guinier, 1955)

.. math::

    P(q) = \frac{\text{scale}}{V_\text{shell}} \left[
           \frac{3V_{\text{core}}({\rho_{\text{solvent}}
           - \rho_{\text{shell}})j_1(qR_{\text{core}})}}{qR_{\text{core}}}
           + \frac{3V_{\text{tot}}(\rho_{\text{shell}}
           - \rho_{\text{solvent}}) j_1(qR_{\text{tot}})}{qR_{\text{tot}}}
           \right]^2 + \text{background}


where scale is a scale factor equivalent to the volume fraction of shell
material if the data is on an absolute scale, $V_{shell}$ is the volume of the
shell, $V_{\text{cor}}$ is the volume of the core, $V_{\text{tot}}$ is the
total volume, $R_{\text{core}}$ is the radius of the core, $R_{\text{tot}}$ is
the outer radius of the shell, $\rho_{\text{solvent}}$ is the scattering length
density of the solvent (which is the same as for the core in this case),
$\rho_{\text{scale}}$ is the scattering length density of the shell, background
is a flat background level (due for example to incoherent scattering in the
case of neutrons), and $j_1$ is the spherical bessel function
$j_1 = (sin(x) - x cos(x))/ x^2$.

The functional form is identical to a "typical" core-shell structure, except
that the scattering is normalized by the volume that is contributing to the
scattering, namely the volume of the shell alone, the scattering length density
of the core is fixed the same as that of the solvent, the scale factor when the
data are on an absolute scale is equivalent to the volume fraction of material
in the shell rather than the entire core+shell sphere, and the parameterization
is done in terms of the core radius = $R_{\text{core}}$ and the shell
thickness = $R_{\text{tot}} - R_{\text{core}}$.

.. figure: img/vesicle_geometry.jpg

The 2D scattering intensity is the same as *P(q)* above, regardless of the
orientation of the *q* vector which is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


NB: The outer most radius (= *radius* + *thickness*) is used as the effective
radius for *S(Q)* when *P(Q)* \* *S(Q)* is applied.

.. image:: img/vesicle_1d.jpg

*Figure. 1D plot using the default values given in the table
(w/200 data point). Polydispersity and instrumental resolution normally
will smear out most of the rapidly oscillating features.*

REFERENCE

A Guinier and G. Fournet, *Small-Angle Scattering of X-Rays*, John Wiley and
Sons, New York, (1955)
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
        solvent_sld: the solvent (and core) SLD
        background: incoherent background
        scale : scale factor = shell volume fraction if on absolute scale"""
category = "shape:sphere"

#             [ "name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld", "1e-6/Ang^2", 0.5, [-inf, inf], "",
               "vesicle shell scattering length density"],
              ["solvent_sld", "1e-6/Ang^2", 6.36, [-inf, inf], "",
               "solvent scattering length density"],
              ["radius", "Ang", 100, [0, inf], "volume",
               "vesicle core radius"],
              ["thickness", "Ang", 30, [0, inf], "volume",
               "vesicle shell thickness"],
             ]

source = ["lib/sph_j1c.c", "vesicle.c"]

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

    whole = 4. * pi * (radius + thickness) ** 3. / 3.
    core = 4. * pi * radius ** 3. / 3.
    return whole, whole - core


# parameters for demo
demo = dict(scale=1, background=0,
            sld=0.5, solvent_sld=6.36,
            radius=100, thickness=30,
            radius_pd=.2, radius_pd_n=10,
            thickness_pd=.2, thickness_pd_n=10)

# For testing against the old sasview models, include the converted parameter
# names and the target sasview model name.
oldname = 'VesicleModel'
oldpars = dict(sld='shell_sld', solvent_sld='solv_sld')


# NOTE: test results taken from values returned by SasView 3.1.2
tests = [[{}, 0.0010005303255, 17139.8268799],
         [{}, 0.200027832249, 0.130387268704 ],
         [{}, 'ER', 130.],
         [{}, 'VR', 0.54483386436],
        ]
