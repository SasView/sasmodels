r"""
Calculates the scattering from a **simple cubic lattice** with
paracrystalline distortion. Thermal vibrations are considered to be
negligible, and the size of the paracrystal is infinitely large.
Paracrystalline distortion is assumed to be isotropic and characterized
by a Gaussian distribution.

Definition
----------

The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = \frac{scale}{V_p}V_{lattice}P(q)Z(q)

where scale is the volume fraction of spheres, $V_p$ is the volume of
the primary particle, $V_{lattice}$ is a volume correction for the crystal
structure, $P(q)$ is the form factor of the sphere (normalized), and
$Z(q)$ is the paracrystalline structure factor for a simple cubic structure.

Equation (16) of the 1987 reference is used to calculate $Z(q)$, using
equations (13)-(15) from the 1987 paper for Z1, Z2, and Z3.

The lattice correction (the occupied volume of the lattice) for a simple
cubic structure of particles of radius *R* and nearest neighbor separation *D* is

.. math::

    V_{lattice}=\frac{4\pi}{3}\frac{R^3}{D^3}

The distortion factor (one standard deviation) of the paracrystal is included
in the calculation of $Z(q)$

.. math::

    \Delta a = gD

where *g* is a fractional distortion based on the nearest neighbor distance.

The simple cubic lattice is

.. figure:: img/sc_crystal_geometry.jpg

For a crystal, diffraction peaks appear at reduced q-values given by

.. math::

    \frac{qD}{2\pi} = \sqrt{h^2+k^2+l^2}

where for a simple cubic lattice any h, k, l are allowed and none are
forbidden. Thus the peak positions correspond to (just the first 5)

.. math::
    :nowrap:

    \begin{align*}
    q/q_0 \quad & \quad 1
                & \sqrt{2} \quad
                & \quad  \sqrt{3} \quad
                & \sqrt{4} \quad
                & \quad \sqrt{5}\quad \\
    Indices \quad & (100)
                  & \quad (110) \quad
                  & \quad (111)
                  & (200) \quad
                  & \quad (210)
    \end{align*}

.. note::

    The calculation of *Z(q)* is a double numerical integral that must be
    carried out with a high density of points to properly capture the sharp
    peaks of the paracrystalline scattering.
    So be warned that the calculation is SLOW. Go get some coffee.
    Fitting of any experimental data must be resolution smeared for any
    meaningful fit. This makes a triple integral. Very, very slow.
    Go get lunch!

The 2D (Anisotropic model) is based on the reference below where *I(q)* is
approximated for 1d scattering. Thus the scattering pattern for 2D may not
be accurate. Note that we are not responsible for any incorrectness of the 2D
model computation.

.. figure:: img/sc_crystal_angle_definition.jpg

Reference
---------
Hideki Matsuoka et. al. *Physical Review B,* 36 (1987) 1754-1765
(Original Paper)

Hideki Matsuoka et. al. *Physical Review B,* 41 (1990) 3854 -3856
(Corrections to FCC and BCC lattice structure calculation)

"""

from numpy import inf

name = "sc_paracrystal"
title = "Simple cubic lattice with paracrystalline distortion"
description = """
        P(q)=(scale/Vp)*V_lattice*P(q)*Z(q)+bkg where scale is the volume
        fraction of sphere,
        Vp = volume of the primary particle,
        V_lattice = volume correction for
        for the crystal structure,
        P(q)= form factor of the sphere (normalized),
        Z(q)= paracrystalline structure factor
        for a simple cubic structure.
        [Simple Cubic ParaCrystal Model]
        Parameters;
        scale: volume fraction of spheres
        bkg:background, R: radius of sphere
        dnn: Nearest neighbor distance
        d_factor: Paracrystal distortion factor
        radius: radius of the spheres
        sldSph: SLD of the sphere
        sldSolv: SLD of the solvent
        """
category = "shape:paracrystal"
single = False
# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["dnn",         "Ang",       220.0,  [0.0, inf],  "",            "Nearest neighbor distance"],
              ["d_factor",    "",            0.06, [-inf, inf], "",            "Paracrystal distortion factor"],
              ["radius",      "Ang",        40.0,  [0.0, inf],  "volume",      "Radius of sphere"],
              ["sld",  "1e-6/Ang^2",  3.0,  [0.0, inf],  "",            "Sphere scattering length density"],
              ["sld_solvent", "1e-6/Ang^2",  6.3,  [0.0, inf],  "",            "Solvent scattering length density"],
              ["theta",       "degrees",     0.0,  [-inf, inf], "orientation", "Orientation of the a1 axis w/respect incoming beam"],
              ["phi",         "degrees",     0.0,  [-inf, inf], "orientation", "Orientation of the a2 in the plane of the detector"],
              ["psi",         "degrees",     0.0,  [-inf, inf], "orientation", "Orientation of the a3 in the plane of the detector"],
             ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sph_j1c.c", "lib/sphere_form.c", "lib/gauss150.c", "sc_paracrystal_kernel.c"]

demo = dict(scale=1, background=0,
            dnn=220.0,
            d_factor=0.06,
            radius=40.0,
            sld=3.0,
            sld_solvent=6.3,
            theta=0.0,
            phi=0.0,
            psi=0.0)

oldname = 'SCCrystalModel'

oldpars = dict(sld='sldSph',
               sld_solvent='sldSolv')

tests = [
    # Accuracy tests based on content in test/utest_extra_models.py
    [{}, 0.001, 10.3048],
    [{}, 0.215268, 0.00814889],
    [{}, (0.414467), 0.001313289]
    ]


