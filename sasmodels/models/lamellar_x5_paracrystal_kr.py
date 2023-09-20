# Note: model title and parameter table are inserted automatically
#2345 789012345678901234567890123456789012345678901234567890123456789012
r"""
THESE DOCS ARE INCOMPLETE AND STILL IN PROGRESS
Implementation of Kotlarchyk & Ritzau paracrystal, which has a Lorentz 
curvature term at small Q, polydispersity and diffuseness on the layer.
See J.Appl.Cryst. 24(1991)753-758 with some "corrections" as detailed 
in the .c code - this is the version as used in the "FISH" program.

This x5 version has explicit amounts of N = 1,2,3,4 or 5 layers, 
so "scale" gets even more complicated to interpret!
Thus you may design your own "polydispersity" in N, but beware the best 
fit may be very poorly determined.
NOTE there may also be redundancy in the scale parameters, do not adjust 
too many at once!
This is an all in one P(Q)S(Q) model, ideally the paracrystal S(Q) 
should be separated.

EQUATIONS HERE ARE CURRENTLY FROM PREVIOUS PEDERSEN VERSION - they need 
work!
This model calculates the scattering from a stack of repeating lamellar
structures. The stacks of lamellae (infinite in lateral dimension) are
treated as a paracrystal to account for the repeating spacing. The repeat
distance is further characterized by a Gaussian polydispersity. 
**This model can be used for large multilamellar vesicles.**

Definition
----------

In the equations below,

- *scale* is used instead of the mass per area of the bilayer $\Gamma_m$
  (this corresponds to the volume fraction of the material in the bilayer,
  *not* the total excluded volume of the paracrystal),

- *sld* $-$ *sld_solvent* is the contrast $\Delta \rho$,

- *thickness* is the layer thickness $t$,

- *Nlayers* is the number of layers $N$,

- *d_spacing* is the average distance between adjacent layers
  $\langle D \rangle$, and

- *sigma_d* is the relative standard deviation of the Gaussian
  layer distance distribution $\sigma_D / \langle D \rangle$.


The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = 2\pi\Delta\rho^2\Gamma_m\frac{P_\text{bil}(q)}{q^2} Z_N(q)

The form factor of the bilayer is approximated as the cross section of an
infinite, planar bilayer of thickness $t$ (compare the equations for the
lamellar model).

.. math::

    P_\text{bil}(q) = \left(\frac{\sin(qt/2)}{qt/2}\right)^2

$Z_N(q)$ describes the interference effects for aggregates
consisting of more than one bilayer. The equations used are (3-5)
from the Bergstrom reference:

.. math::


    Z_N(q) = \frac{1 - w^2}{1 + w^2 - 2w \cos(q \langle D \rangle)}
        + x_N S_N + (1 - x_N) S_{N+1}

where

.. math::

    S_N(q) = \frac{a_N}{N}[1 + w^2 - 2 w \cos(q \langle D \rangle)]^2

and

.. math::

    a_N &= 4w^2 - 2(w^3 + w) \cos(q \langle D \rangle) \\
        &\quad - 4w^{N+2}\cos(Nq \langle D \rangle)
        + 2 w^{N+3}\cos[(N-1)q \langle D \rangle]
        + 2w^{N+1}\cos[(N+1)q \langle D \rangle]

for the layer spacing distribution $w = \exp(-\sigma_D^2 q^2/2)$.

Non-integer numbers of stacks are calculated as a linear combination of
the lower and higher values

.. math::

    N_L = x_N N + (1 - x_N)(N+1)

The 2D scattering intensity is the same as 1D, regardless of the orientation
of the $q$ vector which is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


Reference
---------

M Bergstrom, J S Pedersen, P Schurtenberger, S U Egelhaaf,
*J. Phys. Chem. B*, 103 (1999) 9888-9897
"""

import numpy as np
from numpy import inf

name = "lamellar_x5_paracrystal_kr"
title = "Random lamellar sheet plus stacks with paracrystal structure factors"
description = """\
    [Random lamellar sheet plus stacks with paracrystal structure factors]
       work in progress
"""
category = "shape:lamellae"

single = False

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["thickness", "Ang", 33.0, [0, inf], "volume",
               "mean sheet thickness"],
              ["sigma_t_by_t", "", 0.1, [0.0, inf], "",
               "Sigma (relative polydispersity) of the sheet thickness"],
              ["interface_t", "Ang", 0.0, [0.0, inf], "",
               "interfacial thickness, sqrt(2.pi).sigma, damps form factor"],
              ["rsig_lorentz", "Ang", 0.0, [0.0, inf], "",
               "Lorentz factor for curvature or orientation"],
              ["scale_N1", "", 0.5, [0.0, inf], "",
               "relative amount of N = 1, monolayer"],
              ["scale_N2", "", 0.0, [0.0, inf], "",
               "relative amount of N = 2, bilayer"],
              ["scale_N3", "", 0.5, [0.0, inf], "",
               "relative amount of N = 3, three layers"],
              ["scale_N4", "", 0.0, [0.0, inf], "",
               "relative amount of N = 4, four layers"],
              ["scale_N5", "", 0.0, [0.0, inf], "",
               "relative amount of N = 5, five layers"],
              ["d_spacing", "Ang", 75., [0.0, inf], "",
               "lamellar spacing of paracrystal stack"],
              ["sigma_d_by_d", "", 0.05, [0.0, inf], "",
               "Sigma (relative polydispersity) of the lamellar spacing"],
              ["sld", "1e-6/Ang^2", 1.0, [-inf, inf], "sld",
               "layer scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 6.34, [-inf, inf], "sld",
               "Solvent scattering length density"]
             ]


source = ["lamellar_x5_paracrystal_kr.c"]

form_volume = """
    return 1.0;
"""

def random():
    total_thickness = 10**np.random.uniform(2, 4.7)
    scale_N1 = np.random.randint(0, 1)
    scale_N2 = np.random.randint(0, 1)
    scale_N3 = np.random.randint(0, 1)
    scale_N4 = np.random.randint(0, 1)
    scale_N5 = np.random.randint(0, 1)
    thickness = np.random.uniform(10, 1000)
    d_spacing = thickness*np.random.uniform(1.2,5.0)
    # Let polydispersity peak around 15%; 95% < 0.4; max=100%
    sigma_d = np.random.beta(1.5, 7)
    pars = dict(
        thickness=thickness,
        scale_N1=scale_N1,
        scale_N2=scale_N2,
        scale_N3=scale_N3,
        scale_N4=scale_N4,
        scale_N5=scale_N5,
        d_spacing=d_spacing,
        sigma_d=sigma_d,
    )
    return pars

# ER defaults to 0.0
# VR defaults to 1.0

demo = dict(scale=1, background=0,
            thickness=33, scale_N1=0.5, scale_N2=0.0, scale_N3=0.5, scale_N4=0.0, d_spacing=250, sigma_d=0.2,
            sld=1.0, sld_solvent=6.34,
            thickness_pd=0.2, thickness_pd_n=40, plus_monolayer=0.5)

# needs new unit tests here 
'''
tests = [
    [{'scale': 1.0, 'background': 0.0, 'thickness': 33., 'Nlayers': 20.0,
      'd_spacing': 250., 'sigma_d': 0.2, 'sld': 1.0,
      'sld_solvent': 6.34, 'thickness_pd': 0.0, 'thickness_pd_n': 40},
     [0.001, 0.215268], [21829.3, 0.00487686]],
]
'''
# ADDED by: RKH  ON: June 11, 2020

