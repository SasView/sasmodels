r"""

This model provides the form factor, $P(q)$, for a spherical particle with
a core-shell structure, with damping of the scattering for "fuzzy" interfaces. 
The form factor is normalized by the particle volume.

BEWARE - this version of the model is over-parametrized by allowing calculation
of sld_core and sld_shell via local fractions of solvent and a dry sld for each.
Along with the overall "scale" volume fraction you may only be able to adjust
two of these parameters at any one time. Multiple sets of best fit parameters 
may exist.

ALSO - if you multiply this model by an S(Q), any increase in effective 
particle size due to "fuzziness" is ignored in effective radius, beta(Q) etc 
adjustments in the calculation, which proceeds as for normal core shell sphers.

For information about polarised and magnetic scattering, see
the :ref:`magnetism` documentation.

Definition
----------

The 1D scattering intensity is calculated in the following way (Guinier, 1955)

.. math::

    P(q) = \frac{\text{scale}}{V} F^2(q) + \text{background}

where

.. math::

    F(q) = \frac{3}{V_s}\left[
    V_c(\rho_c-\rho_s)\frac{\sin(qr_c)-qr_c\cos(qr_c)}{(qr_c)^3} +
    V_s(\rho_s-\rho_\text{solv})\frac{\sin(qr_s)-qr_s\cos(qr_s)}{(qr_s)^3}
    \right]\exp\left(\frac{-(\sigma_\text{fuzzy}q)^2}{2}\right)


where $V_s$ is the volume of the whole particle, $V_c$ is the volume of the
core, $r_s$ = $radius$ + $thickness$ is the radius of the particle, $r_c$
is the radius of the core, $\rho_c$ is the scattering length density of the
core, $\rho_s$ is the scattering length density of the shell,
$\rho_\text{solv}$, is the scattering length density of the solvent.

From the Steiger et.al. reference:

  The "fuzziness" of the interface is defined by the parameter
  $\sigma_\text{fuzzy}$. The particle radius $R$ represents the radius of the
  particle where the scattering length density profile decreased to 1/2 of the
  core density. $\sigma_\text{fuzzy}$ is the width of the smeared particle
  surface; i.e., the standard deviation from the average height of the fuzzy
  interface. The inner regions of the microgel that display a higher density
  are described by the radial box profile extending to a radius of
  approximately $R_\text{box} \sim R - 2 \sigma$. The profile approaches
  zero as $R_\text{sans} \sim R + 2\sigma$.

In general the "fuzziness" approximation only works well when $\sigma_\text{fuzzy}$ 
is very small compared to the particle size. For exact I(Q) calculation 
using different radial sld profiles see the :ref:`onion` and :ref:`spherical-sld` 
models.

The 2D scattering intensity is the same as $P(q)$ above, regardless of the
orientation of the $q$ vector.

NB: The outer most radius (ie, = radius + thickness) is used as the
effective radius for $S(Q)$ when $P(Q) \cdot S(Q)$ is applied.


References
----------

#.  A Guinier and G Fournet, *Small-Angle Scattering of X-Rays*, John Wiley and Sons, New York, (1955)
#.  M Stieger, J. S Pedersen, P Lindner, W Richtering, *Langmuir*, 20 (2004) 7283-7292

Authorship and Verification
----------------------------

* **Author:**           RKH    create fuzzy version of core_shell_sphere **Date:** May 28, 2020
* **Last Modified by:**
* **Last Reviewed by:**
"""

import numpy as np
from numpy import pi, inf

name = "core_shell_fuzzy_sphere"
title = "Form factor for core shell spherical particle with damping for fuzzy interfaces"
description = """
    Form factor for core shell spherical particle with damping for fuzzy interfaces
"""
category = "shape:sphere"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["radius",      "Ang",        60.0, [0, inf],    "volume", "Sphere core radius"],
              ["thickness",   "Ang",        10.0, [0, inf],    "volume", "Sphere shell thickness"],
              ["f_solv_core",  "",           0.2,  [0.0,0.999999], "",    "local fraction solvent in core"],
              ["f_solv_shell", "",           0.6,  [0.0,0.999999], "",    "local fraction solvent in shell"],
              ["sld_dry_core",    "1e-6/Ang^2", 1.0,  [-inf, inf], "sld",    "dry core scattering length density"],
              ["sld_dry_shell",   "1e-6/Ang^2", 2.0,  [-inf, inf], "sld",    "dry shell scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 3.0,  [-inf, inf], "sld",    "Solvent scattering length density"],
              ["fuzziness",   "Ang",        10.0, [0, inf],    "volume",   "std deviation of convolution for interface (must be << radius)"]
              ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sas_3j1x_x.c", "lib/core_shell.c", "core_shell_fuzzy_sphere.c"]
have_Fq = True
radius_effective_modes = ["outer radius", "core radius"]

def random():
    """Return a random parameter set for the model."""
    outer_radius = 10**np.random.uniform(1.3, 4.3)
    # Use a distribution with a preference for thin shell or thin core
    # Avoid core,shell radii < 1
    radius = np.random.beta(0.5, 0.5)*(outer_radius-2) + 1
    thickness = outer_radius - radius
    pars = dict(
        radius=radius,
        thickness=thickness,
    )
    return pars

tests = [
    [{'radius': 20.0, 'thickness': 10.0}, 0.1, None, None, 30.0, 4.*pi/3*30**3, 1.0],

    # The SasView test result was 0.00169, with a background of 0.001
    [{'radius': 60.0, 'thickness': 10.0, 'f_solv_core':0.0, 'f_solv_shell':0.0, 'sld_dry_core': 1.0, 'sld_dry_shell': 2.0,
      'sld_solvent': 3.0, 'fuzziness': 0.0, 'background': 0.0}, 0.4, 0.000698838],
]
