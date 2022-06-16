r"""

Definition
----------


This model is an extension of the fuzzy_sphere model to a core-double-shell structure.  It considers an inner and an
outer fuzzieness of the sphere. The inner fuzzieness is of importance in the case of a hollow double layered sphere
(removed core).


The scattering amplitude $A(q)$ is calculated as:

.. math::

    A(q) = \Delta \rho_\text{sh,out} V_\text{sh_out} \Phi_\text{sh_out} (q,R_\text{sh,out},\sigma_\text{out})
           + (\Delta \rho_{sh,in} - \Delta \rho_\text{sh,out}) V_\text{sh,in} \Phi_\text{sh,in}
           (q,R_\text{sh,in},\sigma_\text{int}) + (\Delta \rho_\text{core} - \Delta \rho_\text{sh,in})
           V_\text{core} \Phi_\text{core} (q,R_\text{in},\sigma_\text{in})

    R_\text{sh,out} = W_\text{core} + 2\sigma_{in} + W_\text{sh,in} + 2\sigma_\text{int}
                      + W_\text{sh,out} + \sigma_\text{out}

    R_\text{sh,in} = W_\text{core} + 2\sigma_{in} + W_\text{sh,in} + 2\sigma_\text{int}

where $A(q)^2$ is proportional to the form factor $P(q)$. \Delta \rho_\text{sh,out}, \Delta \rho_\text{sh,in} and
\Delta \rho_\text{core} are the differences between the scattering length density of the solvent and the inner shell
(sh,in), the outer shell (sh,out) and the core (core) respectively. $V$ describes their partial volumes and $W$ the
shell widths. $\Phi(q,R,\sigma)$ is the normalized Fourier transform of the radial density profile:

.. math::

    \Phi(q,R,\sigma) = \frac{1}{V_n} \left[ \left( \frac{R}{\sigma^2} + \frac{1}{\sigma} \right)
    \frac{\cos [q(R + \sigma)]}{q^4} + \left( \frac{R}{\sigma^2} - \frac{1}{\sigma} \right)
    \frac{\cos [q(R - \sigma)]}{q^4} - \frac{3\sin [q(R - \sigma)]}{q^5\sigma^2} -
    \frac{2R \cos(qR)}{q^4\sigma^2} + \frac{6\sin (qR)}{q^5\sigma^2} \right]

    V_n = \frac{R^3}{3} + \frac{R\sigma^2}{6}

From the reference:

    The model describes a structure with an interpenetrating layer of core and inner shell of the length
    2$\sigma_\text{in}$, a intermediate interpenetration layer between the shells 2$\sigma_\text{int}$, and a fuzzy
    outer surface with an extension $\sigma_\text{out}$.

References
----------

#. Brugnoni, M., Scotti, A., Rudov, A. A., Gelissen, A. P., Caumanns, T.,
   Radulescu, A., ... & Richtering, W. (2018). Macromolecules, 51(7), 2662-2671.

Authorship and Verification
----------------------------

* **Author:** Tom Höfken **Date:** June 16, 2022
* **Last Modified by:** Tom Höfken **Date:** June 16, 2022
* **Last Reviewed by:**  **Date:**
"""

import numpy as np
from numpy import inf

name = "core_double_fuzzy_shell"
title = "Scattering from spherical particles with a fuzzy core-double-shell structure."
description = """\
    background: background,
    sld_core: the SLD of the core
	sld_inner_shell: the SLD of the inner shell
	sld_outer_shell: the SLD of the outer shell
    sld_solv: the SLD of the solvent
	rad_core0: radius of sphere(core)
	thick_inner_shell: the thickness of the inner shell
	thick_outer_shell: the thickness of the outer shell
	inner_fuzzieness : the fuzzieness between core and inner shell
	intermediate_fuzzieness : the fuzzieness between inner and outer shell
	outer_fuzzieness : the outer fuzzieness


"""
category = "shape:sphere"

# pylint: disable=bad-whitespace,line-too-long
# ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld_core",         "1e-6/Ang^2",  0, [-inf, inf], "sld",    "Core scattering length density"],
              ["sld_inner_shell",  "1e-6/Ang^2",  1, [-inf, inf], "sld",    "Inner shell scattering length density"],
              ["sld_outer_shell",  "1e-6/Ang^2",  1, [-inf, inf], "sld",    "Outer shell scattering length density"],
              ["sld_solvent",      "1e-6/Ang^2",  3, [-inf, inf], "sld",    "Solvent scattering length density"],
              ["core_radius",      "Ang",        60, [0, inf],    "volume", "Core radius"],
              ["thick_inner_shell","Ang",        40, [0, inf],    "volume", "Thickness of the inner shell"],
              ["thick_outer_shell","Ang",        40, [0, inf],    "volume", "Thickness of the outer shell"],
              ["inner_fuzziness",  "Ang",        5, [0, inf],    "volume", "Fuzzieness between core and inner shell"],
              ["inter_fuzziness",  "Ang",        5, [0, inf], "volume", "Fuzzieness between inner and outer shell"],
              ["outer_fuzziness",  "Ang",        5, [0, inf], "volume", "Fuzzieness at the outer shell"],
              ]
# pylint: enable=bad-whitespace,line-too-long

def phi(q,R,sigma):
    phi = (R/sigma**2 + 1/sigma)*(np.cos(q*(R + sigma)))/(q**4) + (R/sigma**2 - 1/sigma)*(np.cos(q*(R - sigma)))/(q**4) \
    - (3*np.cos(q*(R - sigma)))/(q**5*sigma**2) - (2*R*np.cos(q*R))/(q**4*sigma**2) + (6*np.sin(q*R)/(q**5*sigma**2))
    return phi

def I(q,sld_core = 1,
        sld_inner_shell = 1,
        sld_outer_shell = 1,
        sld_solvent =  3,
        core_radius = 60,
        thick_inner_shell = 60,
        thick_outer_shell = 60,
        inner_fuzziness = 10,
        inter_fuzziness = 10,
        outer_fuzziness = 10,):
    contrast_core  = sld_core - sld_solvent
    contrast_outer_shell = sld_outer_shell - sld_solvent
    contrast_inner_shell  = sld_inner_shell - sld_solvent
    radius_inner_shell = core_radius + 2*inner_fuzziness + thick_inner_shell + inter_fuzziness
    radius_outer_shell = core_radius + 2 * inner_fuzziness + thick_inner_shell + inter_fuzziness + thick_outer_shell + outer_fuzziness
    amplitude = contrast_outer_shell * phi(q,radius_outer_shell,outer_fuzziness) \
                + (contrast_inner_shell - contrast_outer_shell) * phi(q,radius_inner_shell,inter_fuzziness) \
                + (contrast_core - contrast_inner_shell) * phi(q,core_radius,inner_fuzziness)
    inten = amplitude**2
    return inten

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
