r"""
An alternative version of $P(q)$ for the core_shell_ellipsoid
having as parameters the core axial ratio X and a shell thickness, 
which are more often what we would like to determine.

This model is also better behaved when polydispersity is applied than the four
independent radii in core_shell_ellipsoid model.

Definition
----------

.. figure:: img/core_shell_ellipsoid_geometry.png

The geometric parameters of this model are

*equat_core =* equatorial core radius *= R_minor_core*

*X_core = polar_core / equat_core = Rmajor_core / Rminor_core*

*T_shell = equat_outer - equat_core = Rminor_outer - Rminor_core*

*XpolarShell = Tpolar_shell / T_shell = (Rmajor_outer - Rmajor_core)/
(Rminor_outer - Rminor_core)*

In terms of the original radii

*polar_core = equat_core * X_core*

*equat_shell = equat_core + T_shell*

*polar_shell = equat_core * X_core + T_shell * XpolarShell*

(where we note that "shell" perhaps confusingly, relates to the outer radius)
When *X_core < 1* the core is oblate; when *X_core > 1* it is prolate.
*X_core = 1* is a spherical core.

For a fixed shell thickness *XpolarShell = 1*, to scale the shell thickness
pro-rata with the radius *XpolarShell = X_core*.

When including an $S(q)$, the radius in $S(q)$ is calculated to be that of
a sphere with the same 2nd virial coefficient of the outer surface of the
ellipsoid. This may have some undesirable effects if the aspect ratio of the
ellipsoid is large (ie, if $X << 1$ or $X >> 1$ ), when the $S(q)$
- which assumes spheres - will not in any case be valid.

If SAS data are in absolute units, and the SLDs are correct, then scale should
be the total volume fraction of the "outer particle". When $S(q)$ is introduced
this moves to the $S(q)$ volume fraction, and scale should then be 1.0,
or contain some other units conversion factor (for example, if you have SAXS data).

References
----------

R K Heenan, 2015, reparametrised the core_shell_ellipsoid model

"""

from numpy import inf, sin, cos, pi

name = "core_shell_ellipsoid_xt"
title = "Form factor for an spheroid ellipsoid particle with a core shell structure."
description = """
        [core_shell_ellipsoid_xt] Calculates the form factor for an spheroid
        ellipsoid particle with a core_shell structure.
        The form factor is averaged over all possible
        orientations of the ellipsoid such that P(q)
        = scale*<f^2>/Vol + bkg, where f is the
        single particle scattering amplitude.
        [Parameters]:
        equat_core = equatorial radius of core,
        x_core = ratio of core polar/equatorial radii,
        t_shell = equatorial radius of outer surface,
        x_polar_shell = ratio of polar shell thickness to equatorial shell thickness,
        sld_core = SLD_core
        sld_shell = SLD_shell
        sld_solvent = SLD_solvent
        background = Incoherent bkg
        scale =scale
        Note:It is the users' responsibility to ensure
        that shell radii are larger than core radii.
        oblate: polar radius < equatorial radius
        prolate :  polar radius > equatorial radius - this new model will make this easier
        and polydispersity integrals more logical (as previously the shell could disappear).
    """
category = "shape:ellipsoid"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["equat_core",    "Ang",       20,   [0, inf],    "volume",      "Equatorial radius of core"],
    ["x_core",        "None",       3,   [0, inf],    "volume",      "Polar radius of core"],
    ["t_shell",       "Ang",       30,   [0, inf],    "volume",      "Equatorial radius of shell"],
    ["x_polar_shell", "",           1,   [0, inf],    "volume",      "Polar radius of shell"],
    ["sld_core",      "1e-6/Ang^2", 2,   [-inf, inf], "",            "Core scattering length density"],
    ["sld_shell",     "1e-6/Ang^2", 1,   [-inf, inf], "",            "Shell scattering length density"],
    ["sld_solvent",   "1e-6/Ang^2", 6.3, [-inf, inf], "",            "Solvent scattering length density"],
    ["theta",         "degrees",    0,   [-inf, inf], "orientation", "Oblate orientation wrt incoming beam"],
    ["phi",           "degrees",    0,   [-inf, inf], "orientation", "Oblate orientation in the plane of the detector"],
    ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sph_j1c.c", "lib/gfn.c", "lib/gauss76.c",
          "core_shell_ellipsoid_xt.c"]

demo = dict(scale=0.05, background=0.001,
            equat_core=20.0,
            x_core=3.0,
            t_shell=30.0,
            x_polar_shell=1.0,
            sld_core=2.0,
            sld_shell=1.0,
            sld_solvent=6.3,
            theta=0,
            phi=0)

oldname = 'CoreShellEllipsoidXTModel'

oldpars = dict(x_polar_shell='XpolarShell',
               x_core='X_core',
               t_shell='T_shell',
               sld_core='sld_core',
               sld_shell='sld_shell',
               sld_solvent='sld_solvent',
               theta='axis_theta',
               phi='axis_phi')

q = 0.1
phi = pi/6
qx = q*cos(phi)
qy = q*sin(phi)

tests = [
    # Accuracy tests based on content in test/utest_coreshellellipsoidXTmodel.py
    [{'equat_core': 200.0,
      'x_core': 0.1,
      't_shell': 50.0,
      'x_polar_shell': 0.2,
      'sld_core': 2.0,
      'sld_shell': 1.0,
      'sld_solvent': 6.3,
      'background': 0.001,
      'scale': 1.0,
     }, 1.0, 0.00189402],

    # Additional tests with larger range of parameters
    [{'background': 0.01}, 0.1, 11.6915],

    [{'equat_core': 20.0,
      'x_core': 200.0,
      't_shell': 54.0,
      'x_polar_shell': 3.0,
      'sld_core': 20.0,
      'sld_shell': 10.0,
      'sld_solvent': 6.0,
      'background': 0.0,
      'scale': 1.0,
     }, 0.01, 8688.53],

    [{'background': 0.001}, (0.4, 0.5), 0.00690673],

    [{'equat_core': 20.0,
      'x_core': 200.0,
      't_shell': 54.0,
      'x_polar_shell': 3.0,
      'sld_core': 20.0,
      'sld_shell': 10.0,
      'sld_solvent': 6.0,
      'background': 0.01,
      'scale': 0.01,
     }, (qx, qy), 0.0100002],
    ]
