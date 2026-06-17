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

where $\Delta\rho$ is the SLD contrast and $F_n$ is the form factor amplitude
of the $n$-th cylindrical arm:

.. math::

    F_n(q, \theta, \varphi) = \text{sinc}\!\left(\frac{q u_n L}{2}\right)
        \cdot \frac{2 J_1(q \mu_n R)}{q \mu_n R}

with $u_n = \hat{q} \cdot \hat{a}_n$ the projection of $\hat{q}$ onto the
arm axis, $\mu_n = \sqrt{1 - u_n^2}$, $L$ the arm length, and $R$ the arm
radius.  Expanding the squared modulus into a double sum and exploiting the
reality of $F_n$ gives:

.. math::

    I(q) = \frac{(\Delta \rho)^2}{4\pi}
        \int \sum_{n=1}^{4}\sum_{m=1}^{4}
        F_n F_m \cos\!\left(\frac{q(u_n-u_m)L}{2}\right)
        \sin\theta \, d\theta \, d\varphi

The cosine factor is the interference term between the centres of arms $n$
and $m$, which are displaced by $\tfrac{L}{2}\hat{a}_n$ from the junction.

Geometry
--------

The four arms are oriented along tetrahedral directions.  With
$A = 109.5°/2$ (the half-angle between arms), the arm unit vectors and the
corresponding projections $u_n$ are

.. math::

    u_n = s_n \cos A \cos\theta + \sin A \sin\theta \cos(\varphi - \varphi_n)

where $(s_n, \varphi_n) = (+1,\ 0),\ (-1,\ \pi/2),\ (+1,\ \pi),\ (-1,\ 3\pi/2)$
for $n = 1, 2, 3, 4$ respectively.

Each arm has length $L$ and radius $R$.

References
----------

#. Seoki Kyoo Seo *Korean J. Chem. Eng.* 34(2017) 1192-1198

Authorship and Verification
----------------------------

* **Author:** Yuhei Yamada (Github user name: Indigo Carmine, https://orcid.org/0009-0003-9780-4135)
* **Last Modified by:**
"""

from numpy import inf

name = "tetrapod"
title = "Tetrapod with four cylindrical arms"
description = """
    Calculates the scattering from a tetrapod structure with four cylindrical
    arms radiating from a central point.
"""
category = "shape:cylinder"

#             [ "name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["length", "Ang", 1000, [0, inf], "volume", "Cylindrical arm length"],
    ["radius", "Ang", 50, [0, inf], "volume", "Cylindrical arm radius"],
    ["sld", "1e-6/Ang^2", 1, [-inf, inf], "sld", "Tetrapod scattering length density"],
    ["sld_solvent", "1e-6/Ang^2", 0, [-inf, inf], "sld", "Solvent scattering length density"],
]

source = ["lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c", "tetrapod.c"]
have_Fq = False
opencl = True


test = [
    [{}, 1.099643429303014e-03, 2.746377685546875e03],
    [{}, 1.033727919136565e-01, 4.374285638332367e-01],
    [{}, 3.145051218117226e-01, 3.527995198965073e-03],
]
