#fcc paracrystal model
#note model title and parameter table are automatically inserted
#note - calculation requires double precision
r"""
Calculates the scattering from a **face-centered cubic lattice** with
paracrystalline distortion. Thermal vibrations are considered to be
negligible, and the size of the paracrystal is infinitely large.
Paracrystalline distortion is assumed to be isotropic and characterized by
a Gaussian distribution.

Definition
----------

The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = \frac{\text{scale}}{V_p} V_\text{lattice} P(q) Z(q)

where *scale* is the volume fraction of spheres, $V_p$ is the volume of
the primary particle, $V_\text{lattice}$ is a volume correction for the crystal
structure, $P(q)$ is the form factor of the sphere (normalized), and $Z(q)$
is the paracrystalline structure factor for a face-centered cubic structure.

Equation (1) of the 1990 reference is used to calculate $Z(q)$, using
equations (23)-(25) from the 1987 paper for $Z1$, $Z2$, and $Z3$.

The lattice correction (the occupied volume of the lattice) for a
face-centered cubic structure of particles of radius $R$ and nearest
neighbor separation $D$ is

.. math::

   V_\text{lattice} = \frac{16\pi}{3}\frac{R^3}{\left(D\sqrt{2}\right)^3}

The distortion factor (one standard deviation) of the paracrystal is
included in the calculation of $Z(q)$

.. math::

    \Delta a = gD

where $g$ is a fractional distortion based on the nearest neighbor distance.

.. figure:: img/fcc_geometry.jpg

    Face-centered cubic lattice.

For a crystal, diffraction peaks appear at reduced q-values given by

.. math::

    \frac{qD}{2\pi} = \sqrt{h^2 + k^2 + l^2}

where for a face-centered cubic lattice $h, k , l$ all odd or all
even are allowed and reflections where $h, k, l$ are mixed odd/even
are forbidden. Thus the peak positions correspond to (just the first 5)

.. math::

    \begin{array}{cccccc}
    q/q_0 & 1 & \sqrt{4/3} & \sqrt{8/3} & \sqrt{11/3} & \sqrt{4} \\
    \text{Indices} & (111)  & (200) & (220) & (311) & (222)
    \end{array}

**NB**: The calculation of $Z(q)$ is a double numerical integral that
must be carried out with a high density of points to properly capture
the sharp peaks of the paracrystalline scattering. So be warned that the
calculation is SLOW. Go get some coffee. Fitting of any experimental data
must be resolution smeared for any meaningful fit. This makes a triple
integral. Very, very slow. Go get lunch!

The 2D (Anisotropic model) is based on the reference below where $I(q)$ is
approximated for 1d scattering. Thus the scattering pattern for 2D may not
be accurate. Note that we are not responsible for any incorrectness of the
2D model computation.

.. figure:: img/parallelepiped_angle_definition.png

    Orientation of the crystal with respect to the scattering plane, when 
    $\theta = \phi = 0$ the $c$ axis is along the beam direction (the $z$ axis).

References
----------

Hideki Matsuoka et. al. *Physical Review B*, 36 (1987) 1754-1765
(Original Paper)

Hideki Matsuoka et. al. *Physical Review B*, 41 (1990) 3854 -3856
(Corrections to FCC and BCC lattice structure calculation)
"""

from numpy import inf, pi

name = "fcc_paracrystal"
title = "Face-centred cubic lattic with paracrystalline distortion"
description = """
    Calculates the scattering from a **face-centered cubic lattice** with paracrystalline distortion. Thermal vibrations
    are considered to be negligible, and the size of the paracrystal is infinitely large. Paracrystalline distortion is
    assumed to be isotropic and characterized by a Gaussian distribution.
    """
category = "shape:paracrystal"

single = False

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["dnn", "Ang", 220, [-inf, inf], "", "Nearest neighbour distance"],
              ["d_factor", "", 0.06, [-inf, inf], "", "Paracrystal distortion factor"],
              ["radius", "Ang", 40, [0, inf], "volume", "Particle radius"],
              ["sld", "1e-6/Ang^2", 4, [-inf, inf], "sld", "Particle scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld", "Solvent scattering length density"],
              ["theta",       "degrees",    60,    [-360, 360], "orientation", "c axis to beam angle"],
              ["phi",         "degrees",    60,    [-360, 360], "orientation", "rotation about beam"],
              ["psi",         "degrees",    60,    [-360, 360], "orientation", "rotation about c axis"]
             ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sas_3j1x_x.c", "lib/gauss150.c", "lib/sphere_form.c", "fcc_paracrystal.c"]

# parameters for demo
demo = dict(scale=1, background=0,
            dnn=220, d_factor=0.06, sld=4, sld_solvent=1,
            radius=40,
            theta=60, phi=60, psi=60,
            radius_pd=.2, radius_pd_n=0.2,
            theta_pd=15, theta_pd_n=0,
            phi_pd=15, phi_pd_n=0,
            psi_pd=15, psi_pd_n=0,
           )
# april 6 2017, rkh add unit tests, NOT compared with any other calc method, assume correct!
# add 2d test later
q =4.*pi/220.
tests = [
    [{ },
     [0.001, q, 0.215268], [0.275164706668, 5.7776842567, 0.00958167119232]],
]
