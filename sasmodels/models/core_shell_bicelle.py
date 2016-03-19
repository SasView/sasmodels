r"""

Definition
----------
This model provides the form factor for a circular cylinder with a core-shell
scattering length density profile. The form factor is normalized by the
particle volume.

.. _core-shell-bicelle-geometry:

.. figure:: img/core_shell_bicelle_geometry.png

    (Graphic from DOI: 10.1039/C0NP00002G)


The output of the 1D scattering intensity function for randomly oriented
cylinders is then given by the equation above.

The *theta* and *phi* parameters are not used for the 1D output.
Our implementation of the scattering kernel and the 1D scattering intensity
use the c-library from NIST.

.. figure:: img/cylinder_angle_definition.jpg

    Definition of the angles for the oriented core shell bicelle tmodel.

.. figure:: img/cylinder_angle_projection.jpg

    Examples of the angles for oriented pp against the detector plane.

References
----------

L A Feigin and D I Svergun,
*Structure Analysis by Small-Angle X-Ray and Neutron Scattering,*
Plenum Press, New York, (1987)

"""

from numpy import inf, sin, cos

name = "core_shell_bicelle"
title = "Circular cylinder with a core-shell scattering length density profile.."
description = """
    P(q,alpha)= scale/Vs*f(q)^(2) + bkg,  where: f(q)= 2(core_sld
    - solvant_sld)* Vc*sin[qLcos(alpha/2)]
    /[qLcos(alpha/2)]*J1(qRsin(alpha))
    /[qRsin(alpha)]+2(shell_sld-solvent_sld)
    *Vs*sin[q(L+T)cos(alpha/2)][[q(L+T)
    *cos(alpha/2)]*J1(q(R+T)sin(alpha))
    /q(R+T)sin(alpha)]

    alpha:is the angle between the axis of
    the cylinder and the q-vector
    Vs: the volume of the outer shell
    Vc: the volume of the core
    L: the length of the core
    shell_sld: the scattering length density
    of the shell
    solvent_sld: the scattering length density
    of the solvent
    bkg: the background
    T: the thickness
    R+T: is the outer radius
    L+2T: The total length of the outershell
    J1: the first order Bessel function
    theta: axis_theta of the cylinder
    phi: the axis_phi of the cylinder...
        """
category = "shape:cylinder"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["radius",         "Ang",       20, [0, inf],    "volume",      "Cylinder core radius"],
    ["rim_thickness",  "Ang",       10, [0, inf],    "volume",      "Rim shell thickness"],
    ["face_thickness", "Ang",       10, [0, inf],    "volume",      "Cylinder face thickness"],
    ["length",         "Ang",      400, [0, inf],    "volume",      "Cylinder length"],
    ["core_sld",       "1e-6/Ang^2", 1, [-inf, inf], "",            "Cylinder core scattering length density"],
    ["face_sld",       "1e-6/Ang^2", 4, [-inf, inf], "",            "Cylinder face scattering length density"],
    ["rim_sld",        "1e-6/Ang^2", 4, [-inf, inf], "",            "Cylinder rim scattering length density"],
    ["solvent_sld",    "1e-6/Ang^2", 1, [-inf, inf], "",            "Solvent scattering length density"],
    ["theta",          "degrees",   90, [-inf, inf], "orientation", "In plane angle"],
    ["phi",            "degrees",    0, [-inf, inf], "orientation", "Out of plane angle"],
    ]

# pylint: enable=bad-whitespace, line-too-long

source = ["lib/Si.c","lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c", "core_shell_bicelle.c"]

demo = dict(scale=1, background=0,
            radius=20.0,
            rim_thickness=10.0,
            face_thickness=10.0,
            length=400.0,
            core_sld=1.0,
            face_sld=4.0,
            rim_sld=4.0,
            solvent_sld=1.0,
            theta=90,
            phi=0)

oldname = 'CoreShellBicelleModel'

oldpars = dict(rim_thickness='rim_thick',
               face_thickness='face_thick',
               theta='axis_theta',
               phi='axis_phi')

qx, qy = 0.4 * cos(90), 0.5 * sin(0)
tests = [
    # Accuracy tests based on content in test/utest_other_models.py
    [{'radius': 20.0,
      'rim_thickness': 10.0,
      'face_thickness': 10.0,
      'length': 400.0,
      'core_sld': 1.0,
      'face_sld': 4.0,
      'rim_sld': 4.0,
      'solvent_sld': 1.0,
      'background': 0.0,
     }, 0.001, 353.550],

    [{'radius': 20.0,
      'rim_thickness': 10.0,
      'face_thickness': 10.0,
      'length': 400.0,
      'core_sld': 1.0,
      'face_sld': 4.0,
      'rim_sld': 4.0,
      'solvent_sld': 1.0,
      'theta': 90.0,
      'phi': 0.0,
      'background': 0.00,
     }, (qx, qy), 24.9167],

    # Additional tests with larger range of parameters
    [{'radius': 3.0,
      'rim_thickness': 100.0,
      'face_thickness': 100.0,
      'length': 1200.0,
      'core_sld': 5.0,
      'face_sld': 41.0,
      'rim_sld': 42.0,
      'solvent_sld': 21.0,
     }, 0.05, 1670.1828],
    ]
