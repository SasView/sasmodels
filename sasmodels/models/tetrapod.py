r"""
Definition
----------

Calculates the scattering from a tetrapod-shaped structure. A tetrapod consists
of four cylindrical arms radiating from a central point, oriented along the
(1,1,1), (-1,-1,1), (-1,1,-1), and (1,-1,-1) directions.

The scattering intensity is calculated as an average over all
orientations:

.. math::

    I(q) = \frac{(\Delta \rho)^2}{4\pi}
        \int_0^{\pi} \int_0^{2\pi}
        \left|\sum_{n=1}^{4} F_n(q, \theta, \varphi)\right|^2
        \sin\theta \, d\theta \, d\varphi

where $F_n$ is the core-shell form factor amplitude of the $n$-th arm:

.. math::

    F_n(q, \theta, \varphi) = \text{sinc}\!\left(\frac{q u_n L}{2}\right)
        \left[
            (\rho_\text{core} - \rho_\text{shell})\, V_\text{core}\,
            \frac{2 J_1(q \mu_n R)}{q \mu_n R}
          + (\rho_\text{shell} - \rho_\text{solvent})\, V_\text{outer}\,
            \frac{2 J_1(q \mu_n (R+t))}{q \mu_n (R+t)}
        \right]

with $u_n = \hat{q} \cdot \hat{a}_n$, $\mu_n = \sqrt{1 - u_n^2}$, $L$ the arm
length, $R$ the arm core radius, $R + t$ the outer radius ($t$ = shell
thickness), and $V_\text{core} = \pi R^2 L$, $V_\text{outer} = \pi (R+t)^2 L$.
Expanding the squared modulus into a double sum gives:

.. math::

    I(q) = \frac{1}{4\pi}
        \int \sum_{n=1}^{4}\sum_{m=1}^{4}
        F_n F_m \cos\!\left(\frac{q(u_n-u_m)L}{2}\right)
        \sin\theta \, d\theta \, d\varphi

The cosine factor is the interference term between the centres of arms $n$
and $m$, which are displaced by $\tfrac{L}{2}\hat{a}_n$ from the junction.

Geometry
--------

The four arms are oriented along tetrahedral directions.  With
$A = 109.5 /2$ (the half-angle between arms), the arm unit vectors and the
corresponding projections $u_n$ are

.. math::

    u_n = s_n \cos A \cos\theta + \sin A \sin\theta \cos(\varphi - \varphi_n)

where $(s_n, \varphi_n) = (+1,\ 0),\ (-1,\ \pi/2),\ (+1,\ \pi),\ (-1,\ 3\pi/2)$
for $n = 1, 2, 3, 4$ respectively.

Each arm has length $L$, core radius $R$, shell thickness $t$, and hence outer
radius $R + t$.

.. note::

    Each of the four arms is modelled as a complete cylinder of length $L$
    extending from the central junction, so the arms overlap near the origin.
    This overlap is neglected in two ways, following the treatment used in the
    reference. First, the particle volume is overestimated (the effect being
    largest when the arms are short and wide), which affects the calculated
    intensity through the volume normalisation. Second, the scattering amplitude
    is the plain sum of the four cylinder amplitudes, so the overlapping region
    near the origin is counted more than once; because this region is small
    compared with the arms, the resulting artefacts are expected to appear
    mainly in the high-$q$ range. Consequently the model is only valid for
    long-arm tetrapods, i.e. for arm lengths much larger than the arm width,
    $L \gg R + t$.

References
----------

#. Seoki Kyoo Seo *Korean J. Chem. Eng.* 34(2017) 1192-1198 DOI:10.1007/s11814-016-0341-x

Authorship and Verification
----------------------------

* **Author:** Yuhei Yamada (Github user name: Indigo Carmine, https://orcid.org/0009-0003-9780-4135)
* **Last Modified by:**
"""

from numpy import inf

name = "tetrapod"
title = "Core-shell tetrapod with four cylindrical arms"
description = """
    Calculates the scattering from a core-shell tetrapod structure with four
    cylindrical arms radiating from a central point. Each arm has a cylindrical
    core and a coaxial shell.
"""
category = "shape:cylinder"

#             [ "name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["length", "Ang", 1000, [0, inf], "volume", "Arm length"],
    ["radius", "Ang", 50, [0, inf], "volume", "Arm core radius"],
    ["thickness", "Ang", 10, [0, inf], "volume", "Arm shell thickness"],
    ["sld_core", "1e-6/Ang^2", 1, [-inf, inf], "sld", "Arm core scattering length density"],
    ["sld_shell", "1e-6/Ang^2", 0.5, [-inf, inf], "sld", "Arm shell scattering length density"],
    ["sld_solvent", "1e-6/Ang^2", 0, [-inf, inf], "sld", "Solvent scattering length density"],
]

source = ["lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c", "tetrapod.c"]
have_Fq = False
opencl = True


radius_effective_modes = [
    "equivalent volume sphere",
    "length of tetrapod arms (L)",
]


test = [
    # thickness=0 reduces to uniform cylinder: contrast = sld_core - sld_solvent = 1
    [{"thickness": 0}, 1.099643429303014e-03, 2.746377685546875e03],
    [{"thickness": 0}, 1.033727919136565e-01, 4.374285638332367e-01],
    [{"thickness": 0}, 3.145051218117226e-01, 3.527995198965073e-03],
    # default core-shell parameters
    [{}, 1.099643429303014e-03, 1.846798618097820e03],
    [{}, 1.033727919136565e-01, 1.811197691988146e-01],
    [{}, 3.145051218117226e-01, 1.080396557066592e-03],
]
