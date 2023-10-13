# Note: model title and parameter table are inserted automatically
#2345 789012345678901234567890123456789012345678901234567890123456789012
r"""

Implementation of Kotlarchyk & Ritzau paracrystal stack model, 
J.Appl.Cryst. 24(1991)753-758, which has Shultz polydisperse sheets 
with diffuse interfaces, multiplied by a paracrystalline struture 
factor. Additionally here a Lorentz term is included at small Q.to 
allow for some curvature or waviness of the layers.

This x5 version has explicit amounts of N = 1,2,3,4 or 5 layers, 
plus an additional stack with an adjustable number of M\_adj layers.

Thus a distribution in N may be designed and/or the number of layers 
fit. Beware that the best fit may be poorly determined, so do not 
adjust too many parameters at once! In practice this model has proven 
to be very effective in systems with a mixture of 
unilammellar and multilammellar vesicles.

The vesicles must be sufficiently large (or flexible) that the 
interference in Q space from opposite faces, across their diameter, is 
well below the minimum Q values of the data. If in doubt about this 
and/or to help understand both the absolute intensity and 
rsig\_lorentz distance, compare results from the :ref:`vesicle` model 
which uses the exact form factor for rigid spheres.

This is a first attempt at the equations, reports of any issues 
with them or the calculations in the code will be gratefully received.

Further explanations of the scale parameters will follow. It will 
be best to leave the overall scale parameter included by sasview set 
to 1.0 or some other convenient value.

Definition
----------

In the equations below,

- *scale* for each of $Ni = 1, 2, 3, 4, 5$ and then $Madj$ layer stacks determines 
   the amount of each in the final intensity,

- *sld* $-$ *sld_solvent* is the contrast $\Delta \rho$,

- *thickness* is the mean layer thickness, the same for all 
  layers $\langle t \rangle$,

- *sigma_t* is the standard deviation of the Schultz
  layer distance distribution from $\sigma_t / \langle t \rangle$.

- *interface_t* is the width of the interfaces, note $\sigma_{interface} = interface\_t/\sqrt{2\pi}$,

- *The number of layers* in each stack is $Ni = 1, 2, 3, 4, 5$, then $Madj >= 2$,

- *d_spacing* is the average distance between adjacent layers
  $d = \langle D \rangle$, its initial value will need to be close to $2\pi/Qpeak$ 
  for a fit to work,

- *sigma_d* is the standard deviation of the Gaussian
  layer distance distribution from $\sigma_D / \langle D \rangle$, and

- *rsig\_lorentz* in the Lorentz term allows approximately for a finite
  radius of curvature and/or lateral waviness in the layers, as its
  value increases from zero the scattering moves rapidly from the 
  $Q^{-2}$ dependence of scatter from an oriented flat sheet to the 
  $Q^{-4}$ dependence from randomly oriented sheets.


The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = \frac{\Delta\rho^2P_{layer}(q)}{4(1+rsig{\_}lorentz^2q^2/2)}( scale\_N1 + \sum_{i=2}^{5} \left[ scale\_Ni.S(q,i)\right] )


This will then be multiplied by the leading sasview scale parameter and a 
flat background added. The form factor of the layer, $P\_{layer}(q)$ follows 
eqn (17) in the reference (but has been re-arranged in the c++ code to help 
avoid computational errors). 

There is some uncertainty regarding the leading factor of 4, which now seems 
to cancel with the 4 in the numerator of the previous equation, see notes 
in c++ code.

.. math::

    P\_{layer}(q) = \frac{4(2\alpha^{z+1}(2\alpha^{(1-z)} -
       ((2\alpha)^2+4)^{-(z-1)/2})}{2z(z-1)}
       \cos((z-1)tan^{-1}(1/\alpha))\exp{(-q^2\sigma_{interface}^2/2)}

where

.. math::

    z+1 = \frac{1}{(\sigma_t/t)^2}

and

.. math::

    \alpha = \frac{(z+1)}{Qt} = \frac{1}{Qt(\sigma_t/t)^2}

The structure factor for each stack of $N$ layers follows equations (10-12) of 
the reference. (Symbols $I$ and $F$ here follow those in the reference, which are 
different to their more usual meanings.)

.. math::

    S(q,N) = Z(Q) + \frac{I_c(q,N)}{N}

where

.. math::

    I_c(q,N) = -2F&\biggl\{(1+F^2)\cos(qd)-2F-F^N\cos((N+1)qd) \\
    & +2F^{(N+1)}\cos(Nqd)-F^{(N+2)}\cos((N-1)qd) \biggr\}/(1-2F\cos(qd)+F^2)^2


where

.. math::

    Z(q) = \frac{(1-F^2)}{(1-2F\cos(qd)+F^2)}

and

.. math::

    F=exp{(-\sigma_d^2q^2/2)}


There is currently no 2d intensity calculation for this model.


Source
------

Reference
---------
"Paracrystal Model of the High-Temperature Lamellar Phase of a Ternary Microemulsion
System", M.Kotlarchyk & S.M.Ritzau, J.Appl.Cryst. 24(1991)753-758

Authorship and Verification
----------------------------

* **Author:** Richard Heenan **Date:** April 24, 2019
* **Last Modified by:Richard Heenan, Sept 2023**
* **Last Reviewed by:**
"""

import numpy as np
from numpy import inf

name = "lamellar_x5_paracrystal_kr"
title = "Random lamellar sheet plus stacks with paracrystal structure factors"
description = """\
    [lamellar sheet plus stacks with paracrystal structure factors]

"""
category = "shape:lamellae"
have_Fq = False
single = False

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["thickness", "Ang", 33.0, [0, inf], "volume",
               "mean sheet thickness"],
              ["sigma_t_by_t", "", 0.1, [0.0, inf], "",
               "Sigma (relative polydispersity) of the sheet thickness"],
              ["interface_t", "Ang", 5.0, [0.0, inf], "",
               "interfacial thickness, sqrt(2.pi).sigma, damps form factor"],
              ["rsig_lorentz", "Ang", 100.0, [0.0, inf], "",
               "Lorentz factor for curvature or orientation"],
              ["scale_N1", "", 2.5, [0.0, inf], "",
               "relative amount of N = 1, monolayer"],
              ["scale_N2", "", 5.0, [0.0, inf], "",
               "relative amount of N = 2, bilayer"],
              ["scale_N3", "", 0.5, [0.0, inf], "",
               "relative amount of N = 3, three layers"],
              ["scale_N4", "", 0.0, [0.0, inf], "",
               "relative amount of N = 4, four layers"],
              ["scale_N5", "", 0.0, [0.0, inf], "",
               "relative amount of N = 5, five layers"],
              ["scale_Madj", "", 0.0, [0.0, inf], "",
               "relative amount of Madj layers"],
              ["Madj", "", 2.0, [2.0, 50], "",
               "adjustable number of layers"],
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
#TODO test this random function
def random():
    scale_N1 = np.random.randint(0, 1)
    scale_N2 = np.random.randint(0, 1)
    scale_N3 = np.random.randint(0, 1)
    scale_N4 = np.random.randint(0, 1)
    scale_N5 = np.random.randint(0, 1)
    scale_Madj = np.random.randint(0, 1)
    Madj = np.random.randint(6, 50)
    thickness = np.random.uniform(10, 1000)
    d_spacing = thickness*np.random.uniform(1.2,5.0)
    # Let polydispersity peak around 15%; 95% < 0.4; max=100%
    sigma_d_by_d = np.random.beta(1.5, 7)
    pars = dict(
        thickness=thickness,
        scale_N1=scale_N1,
        scale_N2=scale_N2,
        scale_N3=scale_N3,
        scale_N4=scale_N4,
        scale_N5=scale_N5,
        scale_Madj=scale_Madj,
        Madj=Madj,
        d_spacing=d_spacing,
        sigma_d=sigma_d_by_d,
    )
    return pars

# ER defaults to 0.0
# VR defaults to 1.0

#TODO check the demo plot
demo = dict(scale=1, background=0,
            thickness=33, scale_N1=0.5, scale_N2=0.0, scale_N3=0.5, scale_N4=0.0, scale_Madj=0.75, Madj=40,
            d_spacing=250, sigma_d_by_d=0.2,
            sld=1.0, sld_solvent=6.34)

#TODO needs new unit tests here  
'''
# this from a different model - edit it!
tests = [
    [{'scale': 1.0, 'background': 0.0, 'thickness': 33., 'Nlayers': 20.0,
      'd_spacing': 250., 'sigma_d': 0.2, 'sld': 1.0,
      'sld_solvent': 6.34, 'thickness_pd': 0.0, 'thickness_pd_n': 40},
     [0.001, 0.215268], [21829.3, 0.00487686]],
]
'''


