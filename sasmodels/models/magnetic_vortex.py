r"""
Definition
----------
This model describes the analytically approximated scattering of a magnetic vortex in a flat ferromagnetic cylinders made of isotropic material [#Metlov2016]_. 
The circular cylinder with radius $R$ and length $L$ is assumed thin enough so that the magnetization vector $\matbf{M}(r)$ only depends on the radial coordinates $r$. The parameter $b$ denotes a vortex-center displacement from the center.
The parameter $R_M$ allows one to describe the magnetization state [#Metlov2006]_. It can be considered as the distance from the center where a magnetic charge is located to which the magnetization points normally. In the vortex state when the magnetization is always tangential to the boundary ($R_M = R$).
For the quasiuniform state, $p > R$. The signature of the vortex core has been neglected, solutions for the highly inhomogeneous core region can be found in the SI of [#Metlov2016]_.
The magnetic scattering length density (SLD) is defined as $\rho_{\mathrm{mag}}=b_H M_S$, where $b_H= 2.91*10^{8}A^{-1}m^{-1}$ and $M_S$ is the saturation magnetisation (in $A/m$).

The magnetic field is oriented with an inclination of $alpha$ to the neutron beam ans rotated by $beta$. The magnetic scattering length density (SLD) is defined as $\rho_{\mathrm{mag}}=b_H M_S$, where $b_H= 2.91*10^{8}A^{-1}m^{-1}$ and $M_S$ is the saturation magnetisation (in $A/m$).
The fraction of "upward" neutrons before ('up_frac_i') and after the sample ('up_frac_f') must range between 0 to 1, with 0.5 denoting an unpolarised beam. Note that a fit may result in a negative magnetic SLD, and hence magnetisation, when the polarisation state is inverted, i.e. if you have analysed for a $I_{00}$ state wheras your data are $I_{11}$. The model allows to construct the 4 spin-resolved cross sections (non-spin-flip I_{00}, I_{11} and spin-flip, here I_{01}=I_{10}), half-polarised SANS (SANSpol, incoming polarised beam I_0 and I_1, no analysis after sample 'up_frac_f'$=0.5$), and unpolarised beam ('up_frac_i'$=$'up_frac_f'$=0.5$). Differences and other combinations between polarised scattering cross section, e.g. to obtain the nuclear-magnetic interference scattering, or subtraction of the residual scattering of the high field reference state can be constructed with a custom model (Fitting> Add/Multiply Model) and using approbriate scales. For dense systems, special care has to be taken as the nculear structure factor (arrangement of particles) does not need to be identical with the magnetic microstructure e.g. local textures and correlations between easy axes (see [#Honecker2020]_ for further details). The use of structure model is therefore strongly discouraged. Better $I_nuc$, $S_K$ and $S_M$ are fit independent from each other in a model-free way. 


Validation
----------

The model needs testing and validation. Your feedback is kindly requested.



References
----------

.. [#Metlov2016] K.L. Metlov and A. Michels, Sci. Reports 6, 25055 (2016)
.. [#Metlov2006] K.L. Metlov, Phys. Rev. Lett. 97, 127205 (2006)




Authorship and Verification
----------------------------

* **Author: Dirk Honecker **Date:** March 07, 2021
* **Last Modified by:**
* **Last Reviewed by:**

"""

import numpy as np
from numpy import pi, inf

name = "magnetic_vortex_in_disc"
title = "SANS by magnetic vortices in thin submicron-sized softferromagnetic cylinders"
description = """
    I(q) = A F_N^2(q)+ C F_N F_M + D F_M^2 
            A: weighting function =1 for unpolarised beam and non-neutron-spin-flip scattering, zero for spin-flip scattering.
            D: weighting function for purely magnetic scattering of the inhomogeneous magnetisation distribution in the disk.
            C: weighting function for nuclear-magnetic interference scattering.
            The weighting function differ for the various possible spin-resolved scattering cross sections.
            F_N: nuclear form factor
            F_M: magnetic form factor
"""
category = "shape:cylinder"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["radius", "Ang", 50., [0, inf], "volume", "Radius of the core"],
              ["length", "Ang", 40., [0, inf], "volume", "Thickness of disc"],    
              ["nuc_sld_core", "1e-6/Ang^2", 1.0, [-inf, inf], "", "Core scattering length density"],
              ["nuc_sld_solvent", "1e-6/Ang^2", 6.4, [-inf, inf], "", "Solvent scattering length density"],
              ["magnetic_sld_disc",  "1e-6/Ang^2", 1.0,  [-inf, inf], "",    "Magnetic scattering length density of disc"],
              ["b", "", 0.0, [0, inf], "",    "vortex center displacement"],
              ["R_m", "Angstroms", 50,  [0, inf], "",    "position of magnetic charge"],                          
              ["up_i",          "None",  0.5, [0, 1], "", "Polarisation incoming beam"],
              ["up_f",          "None",  0.5, [0, 1], "", "Polarisation outgoing beam"],
              ["alpha",          "None",  90, [0, 180], "", "inclination of field to neutron beam"],
              ["beta",          "None",  0, [0, 360], "", "rotation of field around neutron beam"],             
             ]
# pylint: enable=bad-whitespace, line-too-long




source = ["lib/polevl.c","lib/sas_3j1x_x.c", "lib/sas_gamma.c","lib/gauss76.c","lib/sas_J0.c","lib/sas_J1.c","lib/sas_JN.c","jv.c", "struve.c", "magnetic_vortex.c"]
structure_factor = False
have_Fq = False
single=False






# tests = [
    # [{'radius': 20.0, 'thickness': 10.0}, 0.1, None, None, 30.0, 4.*pi/3*30**3, 1.0],

    # # The SasView test result was 0.00169, with a background of 0.001
    # [{'radius': 60.0, 'thickness': 10.0, 'sld_core': 1.0, 'sld_shell': 2.0,
      # 'sld_solvent': 3.0, 'background': 0.0}, 0.4, 0.000698838],
# ]
