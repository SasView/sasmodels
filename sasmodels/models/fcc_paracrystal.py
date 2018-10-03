#fcc paracrystal model
#note model title and parameter table are automatically inserted
#note - calculation requires double precision
r"""
.. warning:: This model and this model description are under review following 
             concerns raised by SasView users. If you need to use this model, 
             please email help@sasview.org for the latest situation. *The 
             SasView Developers. September 2018.*

Definition
----------

Calculates the scattering from a **face-centered cubic lattice** with
paracrystalline distortion. Thermal vibrations are considered to be
negligible, and the size of the paracrystal is infinitely large.
Paracrystalline distortion is assumed to be isotropic and characterized by
a Gaussian distribution.

The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = \frac{\text{scale}}{V_p} V_\text{lattice} P(q) Z(q)

where *scale* is the volume fraction of spheres, $V_p$ is the volume of
the primary particle, $V_\text{lattice}$ is a volume correction for the crystal
structure, $P(q)$ is the form factor of the sphere (normalized), and $Z(q)$
is the paracrystalline structure factor for a face-centered cubic structure.

Equation (1) of the 1990 reference\ [#CIT1990]_ is used to calculate $Z(q)$,
using equations (23)-(25) from the 1987 paper\ [#CIT1987]_ for $Z1$, $Z2$, and
$Z3$.

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

.. note::

  The calculation of $Z(q)$ is a double numerical integral that
  must be carried out with a high density of points to properly capture
  the sharp peaks of the paracrystalline scattering.
  So be warned that the calculation is slow. Fitting of any experimental data
  must be resolution smeared for any meaningful fit. This makes a triple integral
  which may be very slow.

The 2D (Anisotropic model) is based on the reference below where $I(q)$ is
approximated for 1d scattering. Thus the scattering pattern for 2D may not
be accurate particularly at low $q$. For general details of the calculation
and angular dispersions for oriented particles see :ref:`orientation` .
Note that we are not responsible for any incorrectness of the
2D model computation.

.. figure:: img/parallelepiped_angle_definition.png

    Orientation of the crystal with respect to the scattering plane, when
    $\theta = \phi = 0$ the $c$ axis is along the beam direction (the $z$ axis).

References
----------

.. [#CIT1987] Hideki Matsuoka et. al. *Physical Review B*, 36 (1987) 1754-1765
   (Original Paper)
.. [#CIT1990] Hideki Matsuoka et. al. *Physical Review B*, 41 (1990) 3854 -3856
   (Corrections to FCC and BCC lattice structure calculation)

Authorship and Verification
---------------------------

* **Author:** NIST IGOR/DANSE **Date:** pre 2010
* **Last Modified by:** Paul Butler **Date:** September 29, 2016
* **Last Reviewed by:** Richard Heenan **Date:** March 21, 2016
"""

import numpy as np
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

def random():
    # copied from bcc_paracrystal
    radius = 10**np.random.uniform(1.3, 4)
    d_factor = 10**np.random.uniform(-2, -0.7)  # sigma_d in 0.01-0.7
    dnn_fraction = np.random.beta(a=10, b=1)
    dnn = radius*4/np.sqrt(2)/dnn_fraction
    pars = dict(
        #sld=1, sld_solvent=0, scale=1, background=1e-32,
        dnn=dnn,
        d_factor=d_factor,
        radius=radius,
    )
    return pars

# april 10 2017, rkh add unit tests, NOT compared with any other calc method, assume correct!
# TODO: fix the 2d tests
q = 4.*pi/220.
tests = [
    [{}, [0.001, q, 0.215268], [0.275164706668, 5.7776842567, 0.00958167119232]],
    #[{}, (-0.047, -0.007), 238.103096286],
    #[{}, (0.053, 0.063), 0.863609587796],
]
