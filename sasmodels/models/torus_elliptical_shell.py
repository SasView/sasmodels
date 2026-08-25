r"""
Definition
----------

This model computes the scattering from a torus with an elliptical tube
cross-section and a concentric shell.

The geometric parameters are:

* $R$: major (ring) radius of the torus centerline
* $a$: core minor radius in the radial direction
* $\nu_{core}$: aspect ratio of the elliptical core cross-section
* $\nu_{shell}$: aspect ratio of the elliptical shell cross-section
* $t$: shell thickness

.. figure:: img/torus_elliptical_shell_geometry.png

    Schematic geometry of the torus with elliptical tube cross-section.

The core and shell have independent aspect ratios, so the core semi-axes are
$a$ and $\nu_{core} \cdot a$, and the outer semi-axes are $a+t$ and
$\nu_{core} \cdot a + \nu_{shell} \cdot t$.

For a given orientation angle $\theta$ between the torus symmetry axis and
$\vec q$, the kernel evaluates the amplitude using a numerical integral over
the tube section:

.. math::

    F(q,\theta; x, \Delta\rho, \nu)
    = 4\pi\,\Delta\rho\int_{R-x}^{R+x}
      r\,J_0\!\left(qr\sin\theta\right)
      \frac{\sin\!\left(q\,\gamma\cos\theta\right)}{q\cos\theta}\,dr


with

.. math::

    \gamma = \nu\sqrt{x^2-(r-R)^2}


The core-shell amplitude is formed from outer and inner contributions:

.. math::

    F_{cs}(q,\theta) =
    F\!\left(q,\theta;a + t,\rho_{shell}-\rho_{solvent},\nu_{outer}\right)
    -F\!\left(q,\theta;a,\rho_{shell}-\rho_{core},\nu_{core}\right)


where $\nu_{outer} = \frac{\nu_{shell}t + \nu_{core}a}{a+t}$ is the aspect ratio of the outer ellipse.

and the orientationally averaged intensity is

.. math::

    I(q) = \int_0^{\pi/2}\!\left|F_{cs}(q,\theta)\right|^2\sin\theta\,d\theta


References
----------
#.  T. Kawaguchi, *J. Appl. Crystallogr*, 34(2001) 580-584
#.  S. Förster, *J. Phys. Chem.*, 103(1999) 6657-6668

Authorship and Verification
---------------------------

* **Author:** Itsuki Tajima and Yuhei Yamada (Github user name: Indigo Carmine, https://orcid.org/0009-0003-9780-4135)
* **Last Modified by:**
* **Last Reviewed by:**
"""

from numpy import inf

# ---- SasView model metadata -------------------------------------------------
name = "torus_elliptical_shell"
title = "core-shell torus with elliptical cross-section"
description = "Core-shell torus with elliptical tube cross-section"
category = "shape:cylinder"
parameters = [
    # name          units        default  [min,   max]  type     description
    ["radius", "Ang", 100.0, [0, inf], "volume", "Torus major radius R"],
    [
        "radius_core",
        "Ang",
        5.0,
        [0, inf],
        "volume",
        "Elliptical core minor radius a",
    ],
    ["thickness", "Ang", 2.0, [0, inf], "volume", "Shell thickness t"],
    ["core_nu", "", 1.0, [0.1, 10.0], "volume", "Aspect ratio of core cross-section"],
    ["shell_nu", "", 1.0, [0.1, 10.0], "volume", "Aspect ratio of shell cross-section"],
    ["sld_core", "1e-6/Ang^2", 0.0, [-inf, inf], "sld", "Core SLD"],
    ["sld_shell", "1e-6/Ang^2", 1.0, [-inf, inf], "sld", "Shell SLD"],
    ["sld_solvent", "1e-6/Ang^2", 0.0, [-inf, inf], "sld", "Solvent SLD"],
]

valid = "radius >= radius_core + thickness"

# -- tell sasmodels that a C kernel is provided -------------------------------
source = [
    "lib/polevl.c",
    "lib/sas_J0.c",
    "lib/gauss76.c",
    "torus_elliptical_shell.c",
]  # compiled together with default libs

radius_effective_models = ["outer radius", "equivalent volume sphere", "max distance from center", "radius"]
have_Fq = True

tests = [
    [{}, 1.047452236000633e-03, 2.312834927839701e00],
    [{}, 2.099653600438444e-03, 2.287247365381240e00],
    [{}, 4.208826990208894e-03, 2.186990768206136e00],
]
