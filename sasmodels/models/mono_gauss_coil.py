#mono_gauss_coil model
#conversion of DebyeModel.py
#converted by Steve King, Mar 2016
r"""
This Debye Gaussian coil model strictly describes the scattering from
*monodisperse* polymer chains in theta solvents or polymer melts, conditions
under which the distances between segments follow a Gaussian distribution.
Provided the number of segments is large (ie, high molecular weight polymers)
the single-chain form factor P(Q) is that described by Debye (1947).

To describe the scattering from *polydisperse* polymer chains see the
:ref:`poly-gauss-coil` model.

Definition
----------

.. math::

     I(q) = \text{scale} \cdot I_0 \cdot P(q) + \text{background}

where

.. math::

     I_0 &= \phi_\text{poly} \cdot V
            \cdot (\rho_\text{poly} - \rho_\text{solv})^2 \\
     P(q) &= 2 [\exp(-Z) + Z - 1] / Z^2 \\
     Z &= (q R_g)^2 \\
     V &= M / (N_A \delta)

Here, $\phi_\text{poly}$ is the volume fraction of polymer, $V$ is the
volume of a polymer coil, *M* is the molecular weight of the polymer,
$N_A$ is Avogadro's Number, $\delta$ is the bulk density of the polymer,
$\rho_\text{poly}$ is the sld of the polymer, $\rho\text{solv}$ is the
sld of the solvent, and $R_g$ is the radius of gyration of the polymer coil.

The 2D scattering intensity is calculated in the same way as the 1D,
but where the *q* vector is redefined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}

References
----------

P Debye, *J. Phys. Colloid. Chem.*, 51 (1947) 18.

R J Roe, *Methods of X-Ray and Neutron Scattering in Polymer Science*,
Oxford University Press, New York (2000).

http://www.ncnr.nist.gov/staff/hammouda/distance_learning/chapter_28.pdf
"""

import numpy as np
from numpy import inf, exp, errstate

name = "mono_gauss_coil"
title = "Scattering from monodisperse polymer coils"

description = """
    Evaluates the scattering from
    monodisperse polymer chains.
    """
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["i_zero", "1/cm", 70.0, [0.0, inf], "", "Intensity at q=0"],
    ["rg", "Ang", 75.0, [0.0, inf], "", "Radius of gyration"],
    ]
# pylint: enable=bad-whitespace, line-too-long

# NB: Scale and Background are implicit parameters on every model
def Iq(q, i_zero, rg):
    # pylint: disable = missing-docstring
    z = (q * rg)**2

    with errstate(invalid='ignore'):
        inten = (i_zero * 2.0) * (exp(-z) + z - 1.0)/z**2
        inten[q == 0] = i_zero
    return inten
Iq.vectorized = True # Iq accepts an array of q values

def random():
    rg = 10**np.random.uniform(0, 4)
    #rg = 1e3
    pars = dict(
        #scale=1, background=0,
        i_zero=1e7, # i_zero is a simple scale
        rg=rg,
    )
    return pars

demo = dict(scale=1.0, i_zero=70.0, rg=75.0, background=0.0)

# these unit test values taken from SasView 3.1.2
tests = [
    [{'scale': 1.0, 'i_zero': 70.0, 'rg': 75.0, 'background': 0.0},
     [0.0106939, 0.469418], [57.1241, 0.112859]],
    ]
