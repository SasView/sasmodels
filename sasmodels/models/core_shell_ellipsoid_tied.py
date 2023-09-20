r"""
Definition
----------
Core shell ellipsoid with constrained shell thickness 
and contrast, starting from (dry shell)/ core volume ratio and local fraction 
of solvent in shell. It solves a cubic equation on the fly for the shell thickness.

This may be useful for block copolymer or surfactant micelles where to first 
approximation the shell/core volume ratio is constant, so with polydispersity 
of the core the larger particles would have thicker shells. The calculation of 
I(Q) for this model, unlike many other sasviw models, will adjust the shell 
thickness at all times during integrations over polydispersity.**

Sasview currently has no way to tell the user what the values are for the 
shell thickness and the shell sld say at least at the mean particle size.
(This is in part due to the software architecture, which is designed to 
run fast on gpu etc. )

As a work around set $f\_solvent\_in\_shell$ to be negative, then every time
the I(Q) is calulated some diagnostics, including the shell thickness are sent 
to the sasview log file at the bottom of the screen. This will generate a great 
deal of output, so turn off all fitting first! If polydispersity is being used 
the last shown values in the log will likely be for the largest particle size 
sampled by the numerical integration. The I(Q) calculation uses the absolute 
value of $f\_solvent\_in\_shell$, so should not be affected. You may need to 
change the allowed range from the default 0.01 to 0.99. If 
$f\_solvent\_in\_shell  = 0$ then set  it to say -1.0e-12, any value less than 
-1.0e-24 will count.

.. figure:: img/core_shell_ellipsoid_geometry.png

The geometric parameters of this model are shown in the diagram above, which
shows (a) a cut through at the circular equator and (b) a cross section through
the poles, of a prolate ellipsoid.

When *X_core < 1* the core is oblate; when *X_core > 1* it is prolate.
*X_core = 1* is a spherical core.

For a fixed shell thickness *XpolarShell = 1*, to scale the shell thickness
pro-rata with the radius set or constrain *XpolarShell = X_core*.

If SAS data are in absolute units, and the SLDs are correct, then scale should
be the total volume fraction of the "outer particle". When $S(q)$ is introduced
this moves to the $S(q)$ volume fraction, and scale should then be 1.0, or
contain some other units conversion factor (for example, if you have SAXS data).

The calculation of intensity follows that for the solid ellipsoid, but
with separate terms for the core-shell and shell-solvent boundaries.

.. math::

    P(q,\alpha) = \frac{\text{scale}}{V} F^2(q,\alpha) + \text{background}

where

.. math::
    :nowrap:

    \begin{align*}
    F(q,\alpha) = &f(q,radius\_equat\_core,radius\_equat\_core.x\_core,\alpha) \\
    &+ f(q,radius\_equat\_core + thick\_shell,
         radius\_equat\_core.x\_core + \\ 
      &+ thick\_shell.x\_polar\_shell,\alpha)
    \end{align*}

where

.. math::

    f(q,R_e,R_p,\alpha) = \frac{3 \Delta \rho V (\sin[qr(R_p,R_e,\alpha)]
                - \cos[qr(R_p,R_e,\alpha)])}
                {[qr(R_p,R_e,\alpha)]^3}

and

.. math::

    r(R_e,R_p,\alpha) = \left[ R_e^2 \sin^2 \alpha
        + R_p^2 \cos^2 \alpha \right]^{1/2}


$\alpha$ is the angle between the axis of the ellipsoid and $\vec q$,
$V = (4/3)\pi R_pR_e^2$ is the volume of the ellipsoid , $R_p$ is the
polar radius along the rotational axis of the ellipsoid, $R_e$ is the
equatorial radius perpendicular to the rotational axis of the ellipsoid
and $\Delta \rho$ (contrast) is the scattering length density difference,
either $(sld\_core - sld\_shell)$ or $(sld\_shell - sld\_solvent)$, where 

.. math::

   sld\_shell = sld\_dry\_shell . ( 1 - f\_solvent\_in\_shell) + 
                  sld\_solvent . f\_solvent\_in\_shell.

Note that though the parameterisation here allows "solvent in shell" the 
contrast steps are stil sharp, thus the shell has a uniform scattering 
length density. Also that $x\_polar\_shell$ multiplies $thick\_shell$.

For randomly oriented particles:

.. math::

   F^2(q)=\int_{0}^{\pi/2}{F^2(q,\alpha)\sin(\alpha)d\alpha}

For oriented ellipsoids the *theta*, *phi* and *psi* orientation parameters
will appear when fitting 2D data, see the :ref:`elliptical-cylinder` model
for further information.

References
----------
see for example:
Kotlarchyk, M.; Chen, S.-H. J. Chem. Phys., 1983, 79, 2461.
Berr, S.  J. Phys. Chem., 1987, 91, 4760.

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Richard Heenan (reparametrised model) **Date:** 2015
* **Last Reviewed by:** Richard Heenan **Date:** October 6, 2016
"""

import numpy as np
from numpy import inf, sin, cos, pi

name = "core_shell_ellipsoid_tied"
title = "Form factor for an ellipsoidal particle with a constrained core shell structure."
description = """
        [core_shell_ellipsoid_tied] Calculates the form factor for an spheroid
        ellipsoid particle with a core_shell structure.
        The form factor is averaged over all possible
        orientations of the ellipsoid such that P(q)
        = scale*<f^2>/Vol + bkg, where f is the
        single particle scattering amplitude.
    """
category = "shape:ellipsoid"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["radius_equat_core",       "Ang",       20,   [0, inf],   "volume",    "Equatorial radius of core"],
    ["x_core",                  "None",       3,   [0, inf],    "volume",         "axial ratio of core, X = r_polar/r_equatorial"],
    ["vol_dry_shell_over_core", "None",    0.75,   [0, inf],    "volume",         "volume ratio of dry shell to core"],
    ["x_polar_shell",           "None",       1,   [0, inf],    "volume",   "ratio of thickness of shell at pole to that at equator"],
    ["sld_core",                "1e-6/Ang^2", 2,   [-inf, inf], "sld",      "Core scattering length density"],
    ["sld_dry_shell",           "1e-6/Ang^2", 1,   [-inf, inf], "sld",      "Dry shell scattering length density"],
    ["sld_solvent",             "1e-6/Ang^2", 6.3, [-inf, inf], "sld",      "Solvent scattering length density"],
    ["f_solvent_in_shell",      "None",       0.3, [0.0, 0.99], "volume",      "Local volume fraction of solvent in wet shell"],
    ["theta",                   "degrees",    0,   [-360, 360], "orientation", "ellipsoid axis to beam angle"],
    ["phi",                     "degrees",    0,   [-360, 360], "orientation", "rotation about beam"]
    ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sas_3j1x_x.c", "lib/gauss76.c", "lib/cubic_solve_reparam3.c", "lib/cs_ellipsoid_funcs.c",
 "core_shell_ellipsoid_tied.c"]
# TODO get core_shell_ellipsoid to use the new cs/_ellipsoid_funcs
have_Fq = True
radius_effective_modes = [
    "average outer curvature", "equivalent volume sphere",
    "min outer radius", "max outer radius",
    ]


q = 0.1
# tests had in old coords theta=0, phi=0; new coords theta=90, phi=0
qx = q*cos(pi/6.0)
qy = q*sin(pi/6.0)
# TESTS - work in progress!  default params have thick=7.19915, excel says 7.199326;  try 7.88356 with x-polar_shell=0.5
# TODO - add tests with S(Q) included
tests = [
    # first replicate all five tests from core_shell_ellipsoid
    [{'radius_equat_core': 200.0,
      'x_core': 0.1,
      'vol_dry_shell_over_core':1.34375,
      #'thick_shell': 50.0,
      'x_polar_shell': 0.2,
      'sld_core': 2.0,
      'sld_dry_shell': 1.0,
      'sld_solvent': 6.3,
      'f_solvent_in_shell':0.0,
      'background': 0.001,
      'scale': 1.0,
     }, 1.0, 0.00189402],

    # Additional tests with larger range of parameters as per core_shell_ellipsoid
    [{'radius_equat_core': 20.0,
      'x_core': 3.0,
      'vol_dry_shell_over_core':8.3751,
      #'thick_shell': 30.0,
      'x_polar_shell': 1.0,
      'sld_core': 2.0,
      'sld_dry_shell': 1.0,
      'sld_solvent': 6.3,
      'f_solvent_in_shell':0.0,
      'background': 0.01,
      'scale': 1.0,}, 0.1, 11.6915],

    [{'radius_equat_core': 20.0,
      'x_core': 200.0,
      'vol_dry_shell_over_core':13.2444500,
#      'thick_shell': 54.0, - this case has +ve determinant
      'x_polar_shell': 3.0,
      'sld_core': 20.0,
      'sld_dry_shell': 10.0,
      'sld_solvent': 6.0,
      'f_solvent_in_shell':0.0,
      'background': 0.0,
      'scale': 1.0,
     }, 0.01, 8688.53],

    # 2D tests  as per core_shell_ellipsoid
    [{'radius_equat_core': 20.0,
      'x_core': 3.0,
      'vol_dry_shell_over_core':8.375001,
      #'thick_shell': 30.0,
      'x_polar_shell': 1.0,
      'sld_core': 2.0,
      'sld_dry_shell': 1.0,
      'sld_solvent': 6.3,
      'f_solvent_in_shell':0.0,
      'background': 0.001,
      'theta': 90.0,
      'phi': 0.0,
     }, (0.4, 0.5), 0.00690673],

    [{'radius_equat_core': 20.0,
      'x_core': 200.0,
      'vol_dry_shell_over_core':13.2444500,
      #'thick_shell': 54.0,
      'x_polar_shell': 3.0,
      'sld_core': 20.0,
      'sld_dry_shell': 10.0,
      'sld_solvent': 6.0,
      'f_solvent_in_shell':0.0,
      'background': 0.01,
      'scale': 0.01,
      'theta': 90.0,
      'phi': 0.0,
     }, (0.0866025403, 0.05), 0.01000025],

     # 1D using new parameters - not checked by other means
    [{'radius_equat_core': 20.0,
      'x_core': 3.0,
      'vol_dry_shell_over_core':0.75,
#      'thick_shell': try 7.88356 with x-polar_shell=0.5 - this case has +ve determinant
      'x_polar_shell': 0.5,
      'sld_core': 2.0,
      'sld_dry_shell': 1.0,
      'sld_solvent': 6.3,
      'f_solvent_in_shell':0.3,
      'background': 0.01,
      'scale':0.1,
     }, 0.025, 26.45088379],
     # 1D using new parameters - not checked by other means
    [{'radius_equat_core': 20.0,
      'x_core': 3.0,
      'vol_dry_shell_over_core':8.0,
#      'thick_shell': try 50.0 with x-polar_shell=0.0 - this case is quadratic
      'x_polar_shell': 0.0,
      'sld_core': 2.0,
      'sld_dry_shell': 1.0,
      'sld_solvent': 6.3,
      'f_solvent_in_shell':0.2888889,
      'background': 0.01,
      'scale':0.1,
     }, 0.05, 11.83784036],
]
