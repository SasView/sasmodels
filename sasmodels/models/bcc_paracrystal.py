r"""
Definition
----------

Calculates the scattering from a **body-centered cubic lattice** with
paracrystalline distortion. Thermal vibrations are considered to be negligible,
and the size of the paracrystal is infinitely large. Paracrystalline distortion
is assumed to be isotropic and characterized by a Gaussian distribution.

The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = \frac{\text{scale}}{V_p} V_\text{lattice} P(q) Z(q)


where *scale* is the volume fraction of spheres, $V_p$ is the volume of the
primary particle, $V_\text{lattice}$ is a volume correction for the crystal
structure, $P(q)$ is the form factor of the sphere (normalized), and $Z(q)$
is the paracrystalline structure factor for a body-centered cubic structure.

Equation (1) of the 1990 reference\ [#CIT1990]_ is used to calculate $Z(q)$,
using equations (29)-(31) from the 1987 paper\ [#CIT1987]_ for $Z1$, $Z2$, and
$Z3$.

The lattice correction (the occupied volume of the lattice) for a
body-centered cubic structure of particles of radius $R$ and nearest neighbor
separation $D$ is

.. math::

    V_\text{lattice} = \frac{16\pi}{3} \frac{R^3}{\left(D\sqrt{2}\right)^3}


The distortion factor (one standard deviation) of the paracrystal is included
in the calculation of $Z(q)$

.. math::

    \Delta a = g D

where $g$ is a fractional distortion based on the nearest neighbor distance.


.. figure:: img/bcc_geometry.jpg

    Body-centered cubic lattice.

For a crystal, diffraction peaks appear at reduced q-values given by

.. math::

    \frac{qD}{2\pi} = \sqrt{h^2 + k^2 + l^2}

where for a body-centered cubic lattice, only reflections where
$(h + k + l) = \text{even}$ are allowed and reflections where
$(h + k + l) = \text{odd}$ are forbidden. Thus the peak positions
correspond to (just the first 5)

.. math::

    \begin{array}{lccccc}
    q/q_o          &   1   & \sqrt{2} & \sqrt{3} & \sqrt{4} & \sqrt{5} \\
    \text{Indices} & (110) &    (200) & (211)    & (220)    & (310)    \\
    \end{array}

**NB**: The calculation of $Z(q)$ is a double numerical integral that must
be carried out with a high density of points to properly capture the sharp
peaks of the paracrystalline scattering. So be warned that the calculation
is SLOW. Go get some coffee. Fitting of any experimental data must be
resolution smeared for any meaningful fit. This makes a triple integral.
Very, very slow. Go get lunch!

This example dataset is produced using 200 data points,
*qmin* = 0.001 |Ang^-1|, *qmax* = 0.1 |Ang^-1| and the above default values.

The 2D (Anisotropic model) is based on the reference below where $I(q)$ is
approximated for 1d scattering. Thus the scattering pattern for 2D may not
be accurate.

.. figure:: img/bcc_angle_definition.png

    Orientation of the crystal with respect to the scattering plane.

References
----------

.. [#CIT1987] Hideki Matsuoka et. al. *Physical Review B*, 36 (1987) 1754-1765
   (Original Paper)
.. [#CIT1990] Hideki Matsuoka et. al. *Physical Review B*, 41 (1990) 3854 -3856
   (Corrections to FCC and BCC lattice structure calculation)

Authorship and Verification
----------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Paul Butler **Date:** September 29, 2016
* **Last Reviewed by:** Richard Heenan **Date:** March 21, 2016
"""

from numpy import inf, pi

name = "bcc_paracrystal"
title = "Body-centred cubic lattic with paracrystalline distortion"
description = """
    Calculates the scattering from a **body-centered cubic lattice** with
    paracrystalline distortion. Thermal vibrations are considered to be
    negligible, and the size of the paracrystal is infinitely large.
    Paracrystalline distortion is assumed to be isotropic and characterized
    by a Gaussian distribution.
    """
category = "shape:paracrystal"

#note - calculation requires double precision
single = False

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description" ],
parameters = [["dnn",         "Ang",       220,    [-inf, inf], "",            "Nearest neighbour distance"],
              ["d_factor",    "",            0.06, [-inf, inf], "",            "Paracrystal distortion factor"],
              ["radius",      "Ang",        40,    [0, inf],    "volume",      "Particle radius"],
              ["sld",         "1e-6/Ang^2",  4,    [-inf, inf], "sld",         "Particle scattering length density"],
              ["sld_solvent", "1e-6/Ang^2",  1,    [-inf, inf], "sld",         "Solvent scattering length density"],
              ["theta",       "degrees",    60,    [-inf, inf], "orientation", "In plane angle"],
              ["phi",         "degrees",    60,    [-inf, inf], "orientation", "Out of plane angle"],
              ["psi",         "degrees",    60,    [-inf, inf], "orientation", "Out of plane angle"]
             ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sas_3j1x_x.c", "lib/gauss150.c", "lib/sphere_form.c", "bcc_paracrystal.c"]

# parameters for demo
demo = dict(
    scale=1, background=0,
    dnn=220, d_factor=0.06, sld=4, sld_solvent=1,
    radius=40,
    theta=60, phi=60, psi=60,
    radius_pd=.2, radius_pd_n=2,
    theta_pd=15, theta_pd_n=0,
    phi_pd=15, phi_pd_n=0,
    psi_pd=15, psi_pd_n=0,
    )
# april 6 2017, rkh add unit tests, NOT compared with any other calc method, assume correct!
q =4.*pi/220.
tests = [
    [{ },
     [0.001, q, 0.215268], [1.46601394721, 2.85851284174, 0.00866710287078]],
]
