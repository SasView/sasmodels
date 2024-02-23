r"""

Definition
----------
This plug-in model calculates chains of oriented parallelepipeds, with the option of 
adding a magnetic SLD. The chain scattering is the incoherent 
sum of a user-defined combination of singletons, dimers, trimers, 
quadramers, and pentamers. Note that no matter the numerical values 
selected for the amount of each chain type, the fraction of each will be 
normalized such that the sum of the chain type fractions is unity. 

The chains are oriented along the x-direction. 

The magnetism of the chains is given by a uniform magnetic sld within each parallelepiped particle.

The scattering amplitude form factor is calculated in same way as the parallelepiped (s. https://www.sasview.org/docs/user/models/parallelepiped.html), and it is then multiplied by a complex structure factor that depends on chain length:

.. math::

    I(q)= scale*V*((nuc_sld - sld_solvent)+mag_sld)^2*P(q,alpha)*(  ext{real  phase}^2 + 	ext{img  phase}^2)  +background
        P(q,alpha) = integral from 0 to 1 of ...
           phi(mu*sqrt(1-sigma^2),a) * S(mu*c*sigma/2)^2 * dsigma
        with
            phi(mu,a) = integral from 0 to 1 of ..
            (S((mu/2)*cos(pi*u/2))*S((mu*a/2)*sin(pi*u/2)))^2 * du
            S(x) = sin(x)/x
            mu = q*B
        V: Volume of the rectangular parallelepiped
        alpha: angle between the long axis of the
            parallelepiped and the q-vector for 1D

and
.. math::
    	ext{real  phase} =  1.0 + sum_{k=0}^{4} sum_{n=0}^{k} cos(k*Length*(Q_X*ChainProjX))
and
.. math::
    	ext{img  phase} =  0.0 + sum_{k=0}^{4} sum_{n=0}^{k} sin(k*Length*(Q_X*ChainProjX ))


Here $V$ is the volume of one parallelepiped, $nuc_sld$ is the univorm nuclear sld of one parallelepiped, $sld_solvent$ the surrounding solvent sld, $sld_mag$ the uniform magnetic sld, $k$ defines the number of parallelepipeds in the chain which is oriented along Q_X. This model is a combination of a single parallelepiped (s. https://www.sasview.org/docs/user/models/parallelepiped.html) together with the chain model described in https://marketplace.sasview.org/models/143/, but with fixed oriented direction along x, without polydisperse orientation of the chain, and with a uniform magnetic sld.
"""
import numpy as np
from sasmodels.special import *
from numpy import inf

name = "ParallelepipedChain"
title = "Base-script: Rectangular parallelepiped. Add-on: Chain of parallelepipeds along x-direction"
description = """User model for chains of parallelepipeds oriented along X-axis"""

category = "shape:parallelepiped"

parameters = [["sld", "1e-6/Ang^2", 4, [-inf, inf], "sld",
               "Parallelepiped scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["length_a", "Ang", 35, [0, inf], "length",
               "Shorter side of the parallelepiped"],
              ["length_b", "Ang", 75, [0, inf], "length",
               "Second side of the parallelepiped"],
              ["length_c", "Ang", 400, [0, inf], "length",
               "Larger side of the parallelepiped"],
              ["Length", "Ang", 400, [0, inf], "length",
               "Length of the distance between parallelepipeds"],
              ['singlets', 'Fraction of singlets', 1.0, [0, 100], '', ''],
              ['doublets', 'Fraction of doubles', 1.0, [0, 100], '', ''],
              ['trimers', 'Fraction of trimers', 1.0, [0, 100], '', ''],
              ['quadramers', 'Fraction of quadramers', 1.0, [0, 100], '', ''],
              ['pentamers', 'Fraction of pentamers', 1.0, [0, 100], '', ''],
             ]

source = ["lib/gauss76.c", "ParallelepipedChain.c"]
have_Fq = True
radius_effective_modes = [
    "equivalent cylinder excluded volume", "equivalent volume sphere",
    "half length_a", "half length_b", "half length_c",
    "equivalent circular cross-section", "half ab diagonal", "half diagonal",
    ]

def random():
    """Return a random parameter set for the model."""
    length = 10**np.random.uniform(1, 4.7, size=3)
    pars = dict(
        length_a=length[0],
        length_b=length[1],
        length_c=length[2],
    )
    return pars



