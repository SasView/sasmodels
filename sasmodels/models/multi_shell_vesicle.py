r"""
Definition
----------

This model calculates the scattering from a spherical, solvent-filled core
surrounded by up to 10 concentric shells of individually specified thickness
and scattering length density (a "multi-shell vesicle").

It is a trivial extension of the :ref:`vesicle` model to a larger number of
shells. As for the vesicle model, the SLD of the core is constrained to be
equal to that of the solvent, and the volume fraction refers to the material
in the shells (not the full enclosed sphere).

The 1D scattering intensity is calculated as

.. math::

    I(q) = \frac{\phi_s}{V_s}\,\Big[
        \sum_{k=1}^{N}\, V(r_k)\,(\rho_k - \rho_{k-1})\,\frac{3 j_1(q r_k)}{q r_k}
    \Big]^2 + \text{background}

where $j_1(x) = (\sin x - x\cos x)/x^2$ is the spherical Bessel function,
$\rho_0 = \rho_{solvent}$ is the solvent (and core) SLD, $\rho_k$ is the SLD of
shell $k$, $r_k$ is the radius at the outer edge of shell $k$, $V(r_k) =
4\pi r_k^3/3$ is the volume enclosed by that radius, $V_s$ is the total volume
of the shell material, and $\phi_s$ is the volume fraction of shell material.
The last term of the sum (with $\rho_{k} = \rho_{solvent}$) accounts for the
contrast between the outermost shell and the solvent.

The 2D scattering intensity is the same as $P(q)$ above, regardless of the
orientation of the $\vec q$ vector which is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}

.. note:: **Be careful!** The SLDs and volfraction can be highly correlated.
    Hold as many of these parameters fixed as possible. Since *volfraction*
    is part of the model, *scale* should be fixed at 1.0.

.. note:: The number of shells is given by the parameter *n*. For $n = 1$ and
    a uniform shell SLD this model is identical to the :ref:`vesicle` model.

.. note:: The outer most radius (= *radius* + sum of *thickness_k*) is used as
    the effective radius for $S(Q)$ when $P(Q)*S(Q)$ is applied (mode 1); the
    core radius is available as mode 2.

For information about polarised and magnetic scattering, see
the :ref:`magnetism` documentation.

References
----------

.. [#Kline2006] S R Kline, *J Appl. Cryst.*, 39 (2006) 895
.. [#Feigin1987] L A Feigin and D I Svergun, *Structure Analysis by
   Small-Angle X-Ray and Neutron Scattering*, Plenum Press, New York, 1987.

Authorship and Verification
----------------------------
* **Author:** NIST IGOR/DANSE/Paul Kienzle **Date:** September 12, 2016
* **Last Modified by:** Lionel Porcar **Date:** May 3, 2023
* **Last Reviewed by:** Ursula Perez-Salas **Date:** November 4, 2026
"""

import numpy as np
from numpy import inf

name = "multi_shell_vesicle"
title = "Scattering from a solvent-core sphere with 1 to 10 concentric shells \
(multi-shell vesicle). Each shell can have a unique thickness and SLD. The \
core SLD is equal to the solvent SLD and volfraction refers to the shell \
material."

description = """
Form factor for a multi-shell (up to 10) vesicle with a solvent-filled core,
normalized by the volume of the shell material.
Each shell can have a unique thickness and sld.

    volfraction: volume fraction of shell material
    radius: the radius of the (solvent-filled) core
    n: number of shells
    thickness_k: the thickness of the k'th shell from the core
    sld_k: the SLD of the k'th shell from the core
    sld_solvent: the SLD of the solvent (and of the core)
    background: incoherent background
    scale: should be fixed at 1 (volfraction is part of the model)
"""

category = "shape:sphere"


#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["radius", "Ang", 100., [0, inf], "volume",
               "Radius of the solvent-filled core"],
              ["sld_solvent", "1e-6/Ang^2", 6.36, [-inf, inf], "sld",
               "Solvent (and core) scattering length density"],
              ["n", "", 1, [1, 10], "volume",
               "number of shells"],
              ["sld[n]", "1e-6/Ang^2", 0.5, [-inf, inf], "sld",
               "scattering length density of shell k"],
              ["thickness[n]", "Ang", 30., [0, inf], "volume",
               "Thickness of shell k"],
              ["volfraction", "", 0.05, [0, 1.0], "",
               "volume fraction of shell material"],
             ]

source = ["lib/sas_3j1x_x.c", "multi_shell_vesicle.c"]

have_Fq = True
radius_effective_modes = ["outer radius", "core radius"]


def random():
    """Return a random parameter set for the model."""
    num_shells = np.minimum(np.random.poisson(3)+1, 10)
    total_radius = 10**np.random.uniform(1.7, 4)
    thickness = np.random.exponential(size=num_shells+1)
    thickness *= total_radius/np.sum(thickness)
    pars = dict(
        #background=0,
        scale=1,  # volfraction is part of the model, so scale=1
        volfraction=10**np.random.uniform(-3, -1),
        n=num_shells,
        radius=thickness[0],
    )
    for k, v in enumerate(thickness[1:]):
        pars['thickness%d' % (k+1)] = v
    return pars


def profile(radius, sld_solvent, n, sld, thickness, volfraction):
    """
    Returns the SLD profile *r* (Ang), and *rho* (1e-6/Ang^2).
    """
    n = int(n+0.5)
    z = []
    rho = []

    # add in the solvent-filled core (SLD = solvent)
    z.append(0)
    rho.append(sld_solvent)
    z.append(radius)
    rho.append(sld_solvent)

    # add in the shells
    for k in range(int(n)):
        # Left side of each shell
        z.append(z[-1])
        rho.append(sld[k])
        z.append(z[-1] + thickness[k])
        rho.append(sld[k])
    # add in the solvent
    z.append(z[-1])
    rho.append(sld_solvent)
    z.append(z[-1]*1.25)
    rho.append(sld_solvent)

    return np.asarray(z), np.asarray(rho)


# Equivalence with the vesicle model for a single shell (n=1): these are the
# vesicle.py test values, computed with the vesicle defaults
# (radius=100, thickness=30, sld=0.5, sld_solvent=6.36, volfraction=0.05,
# background=0.001). The multi-shell model with n=1 and matching parameters
# must reproduce them exactly. VERIFY with `sascomp` before submitting.
_VESICLE_DEFAULTS = dict(radius=100., n=1, thickness1=30., sld1=0.5,
                         sld_solvent=6.36, volfraction=0.05)

tests = [
    [_VESICLE_DEFAULTS, 0.0005, 859.916526646],
    [_VESICLE_DEFAULTS, 0.100600200401, 1.77063682331],
    [_VESICLE_DEFAULTS, 0.5, 0.00355351388906],
    # R_eff (outer radius mode), form:shell volume ratio for n=1 vesicle case
    [_VESICLE_DEFAULTS, 0.1, None, None, 130., None, 1./0.54483386436],
]
