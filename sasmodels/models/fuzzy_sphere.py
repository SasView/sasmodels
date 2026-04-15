r"""
For information about polarised and magnetic scattering, see
the :ref:`magnetism` documentation.

Definition
----------
Scattering from spheres with a Gaussian interface.

The scattering intensity $I(q)$ is calculated as:

.. math::

    I(q) = \frac{\text{scale}}{V}(\Delta \rho)^2 A(q)^2 S(q)
           + \text{background}

where the amplitude $A(q)$ describes the radial scattering length density
profile of a homogeneous sphere convoluted with a Gaussian in order to give
a function with a gradual drop-off in scattering length density (SLD) towards
the interface (i.e. a sphere with a diffuse, or "fuzzy", interface):

.. math::

    A(q) = \frac{3\left[\sin(qR) - qR \cos(qR)\right]}{(qR)^3}
           \exp\left(\frac{-(\sigma_\text{fuzzy}q)^2}{2}\right)

Here $A(q)^2$ is the form factor, $P(q)$. The $scale$ is equivalent to the
volume fraction of spheres, each of volume, $V$. And the contrast $(\Delta \rho)$
is the difference in SLD between a sphere and the surrounding medium.

In this model, **$R$ represents the radius at which the SLD has decreased to
half of its value at the core, not the overall radius of a sphere**. This is
a frequent source of confusion when applying this model.
  
$\sigma_\text{fuzzy}$ is then the width of the fuzzy interface; strictly, the
standard deviation from the average thickness of the interface.

From Reference [1]:

  The inner regions ... that display a higher ... [SLD] are described by
  a radial box profile extending to a radius of approximately $R_\text{box} \sim R - 2 \sigma_\text{fuzzy}$.
  The profile approaches zero at $R_\text{SANS} \sim R + 2 \sigma_\text{fuzzy}$. Therefore,
  the overall size of the fuzzy sphere is approximated by ... $R_\text{SANS}$.

For this model to give meaningful results it is important $\sigma_\text{fuzzy} \ll R$.
It is for the User to ensure that this condition is maintained, especially if
applyingnpolydispersity to one or both length scales.
    
This model has been widely applied to the scattering from polymer microgel
particles as illustrated below, where $R_\text{h}$ is the hydrodynamic radius.

.. figure:: img/fuzzy_sphere_geometry.jpg

    fuzzy_sphere model applied to a microgel particle (adapted [2], Fig 5).

Although the fuzzy sphere model often provides a good description of scattering
data from such systems, advances in measurement techniques have highlighted
that the $real-space$ density profile can be far more complex than this model
assumes [3].

This model is not suitable for describing spherical particles decorated with
so-called polymer 'brushes' (where the SLD profile follows a parabolic decay) or
spherical particles with terminally-attached polymer chains (where the SLD
profile is expected to exhibit a maximum before the Gaussian decay).

To model more complex SLD profiles, see the :ref:`onion` and :ref:`spherical_sld` models.

For 2D data: The 2D scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math:: q = \sqrt{{q_x}^2 + {q_y}^2}

References
----------

#. M Stieger, J. S Pedersen, P Lindner, W Richtering,
   *Langmuir*, 20 (2004) 7283-7292

#. E Ponomareva, B Tadgell, M Hildebrandt, M Krüsmann, S Prévost, P Mulvaney, M Karg,
   *Soft Matter*, 18 (2022) 807-825
   
#. F Scheffold,
   *Soft Matter*, 20 (2024) 8181-8184

Authorship and Verification
----------------------------

* **Author:**
* **Last Modified by: Steve King Date: April 15, 2026**
* **Last Reviewed by:**
"""

import numpy as np
from numpy import inf

name = "fuzzy_sphere"
title = "Scattering from spheres with a Gaussian interface."
description = """\
scale: scale factor times volume fraction,
or just volume fraction for absolute scale data
radius: pseudo-radius of the sphere (read the model help)
fuzziness = standard deviation of the fuzzy interfacial
thickness
sld: the SLD of the sphere
sld_solvent: the SLD of the solvent
background: background
Note: By definition, this model only works when fuzziness << radius.
"""
category = "shape:sphere"

# pylint: disable=bad-whitespace,line-too-long
# ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld",         "1e-6/Ang^2",  1, [-inf, inf], "sld",    "Particle scattering length density"],
              ["sld_solvent", "1e-6/Ang^2",  3, [-inf, inf], "sld",    "Solvent scattering length density"],
              ["radius",      "Ang",        60, [0, inf],    "volume", "Particle pseudo-radius (read the model help)"],
              ["fuzziness",   "Ang",        10, [0, inf],    "volume",       "Std deviation from average thickness of the fuzzy interface (must be << radius)"],
             ]
# pylint: enable=bad-whitespace,line-too-long

source = ["lib/sas_3j1x_x.c", "fuzzy_sphere.c"]
have_Fq = True
radius_effective_modes = ["radius", "radius + fuzziness"]

def random():
    """Return a random parameter set for the model."""
    radius = 10**np.random.uniform(1, 4.7)
    fuzziness = 10**np.random.uniform(-2, -0.5)*radius  # 1% to 31% fuzziness
    pars = dict(
        radius=radius,
        fuzziness=fuzziness,
    )
    return pars

tests = [
    # Accuracy tests based on content in test/utest_models_new1_3.py
    #[{'background': 0.001}, 1.0, 0.001],

    [{}, 0.00301005, 359.2315],

    ]
