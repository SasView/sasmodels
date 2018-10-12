r"""
.. warning:: This model and this model description are under review following 
             concerns raised by SasView users. If you need to use this model, 
             please email help@sasview.org for the latest situation. *The 
             SasView Developers. September 2018.*
             
Definition
----------

Calculates the scattering from a **simple cubic lattice** with
paracrystalline distortion. Thermal vibrations are considered to be
negligible, and the size of the paracrystal is infinitely large.
Paracrystalline distortion is assumed to be isotropic and characterized
by a Gaussian distribution.

The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = \text{scale}\frac{V_\text{lattice}P(q)Z(q)}{V_p} + \text{background}

where scale is the volume fraction of spheres, $V_p$ is the volume of
the primary particle, $V_\text{lattice}$ is a volume correction for the crystal
structure, $P(q)$ is the form factor of the sphere (normalized), and
$Z(q)$ is the paracrystalline structure factor for a simple cubic structure.

Equation (16) of the 1987 reference\ [#CIT1987]_ is used to calculate $Z(q)$,
using equations (13)-(15) from the 1987 paper\ [#CIT1987]_ for $Z1$, $Z2$, and
$Z3$.

The lattice correction (the occupied volume of the lattice) for a simple cubic
structure of particles of radius *R* and nearest neighbor separation *D* is

.. math::

    V_\text{lattice}=\frac{4\pi}{3}\frac{R^3}{D^3}

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
    So be warned that the calculation is slow. Fitting of any experimental data
    must be resolution smeared for any meaningful fit. This makes a triple
    integral which may be very slow.

The 2D (Anisotropic model) is based on the reference below where *I(q)* is
approximated for 1d scattering. Thus the scattering pattern for 2D may not
be accurate particularly at low $q$. For general details of the calculation
and angular dispersions for oriented particles see :ref:`orientation` .
Note that we are not responsible for any incorrectness of the
2D model computation.

.. figure:: img/parallelepiped_angle_definition.png

    Orientation of the crystal with respect to the scattering plane, when
    $\theta = \phi = 0$ the $c$ axis is along the beam direction (the $z$ axis).

Reference
---------

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
parameters = [["dnn",         "Ang",       220.0, [0.0, inf],  "",            "Nearest neighbor distance"],
              ["d_factor",    "",           0.06, [-inf, inf], "",            "Paracrystal distortion factor"],
              ["radius",      "Ang",        40.0, [0.0, inf],  "volume",      "Radius of sphere"],
              ["sld",  "1e-6/Ang^2",         3.0, [0.0, inf],  "sld",         "Sphere scattering length density"],
              ["sld_solvent", "1e-6/Ang^2",  6.3, [0.0, inf],  "sld",         "Solvent scattering length density"],
              ["theta",       "degrees",    0,    [-360, 360], "orientation", "c axis to beam angle"],
              ["phi",         "degrees",    0,    [-360, 360], "orientation", "rotation about beam"],
              ["psi",         "degrees",    0,    [-360, 360], "orientation", "rotation about c axis"]
             ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sas_3j1x_x.c", "lib/sphere_form.c", "lib/gauss150.c", "sc_paracrystal.c"]

def random():
    # copied from bcc_paracrystal
    radius = 10**np.random.uniform(1.3, 4)
    d_factor = 10**np.random.uniform(-2, -0.7)  # sigma_d in 0.01-0.7
    dnn_fraction = np.random.beta(a=10, b=1)
    dnn = radius*4/np.sqrt(4)/dnn_fraction
    pars = dict(
        #sld=1, sld_solvent=0, scale=1, background=1e-32,
        dnn=dnn,
        d_factor=d_factor,
        radius=radius,
    )
    return pars

tests = [
    # Accuracy tests based on content in test/utest_extra_models.py, 2d tests added April 10, 2017
    [{}, 0.001, 10.3048],
    [{}, 0.215268, 0.00814889],
    [{}, 0.414467, 0.001313289],
    [{'theta': 10.0, 'phi': 20, 'psi': 30.0}, (0.045, -0.035), 18.0397138402],
    [{'theta': 10.0, 'phi': 20, 'psi': 30.0}, (0.023, 0.045), 0.0177333171285],
    ]
