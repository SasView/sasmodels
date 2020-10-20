r"""
Definition
----------
This model describes the SANS of individual (dilute), superparamagnetic particles. The magnetic field is oriented with an inclination of $alpha$ to the neutron beam ans rotated by $beta$.
The model is based on spherical particles with a core multishell structure, see :ref:`core_multi_shell`.
The magnetic scattering length density (SLD) is defined as $\rho_{\mathrm{mag}}=b_H M_S$, 
where $b_H= 2.91*10^{8}A^{-1}m^{-1}$ and $M_S$ is the saturation magnetisation (in $A/m$).

The alignment of the particle moments along the magnetic field is disturbed by thermal fluctuations like Brownian motion of freely rotating particles in a liquid. 
The magnetic SLD $\rho_{\mathrm{mag}}_mag=b_H M_S$ is hence weighted for each magnetisation component with the corresponding average of a Boltzmann orientation distribution 
with $\eta=\mu_0 \mu H/ k_B T$ denotes the Langevin parameter with $\mu_0$ the vacuum permeability and $\mu$ the apparent magnetic moment, and $k_B T$ the thermal energy  [#Wiedenmann2005]_.
The model describes a magnetisation distribution that is equally distributed around the magnetic field direction (no effective prefered easy axis of the ensemble).
In particular, the transversal magnetisation components are randomly oriented at least from one particle to the other.

The resulting magnetisation averages are:

.. math::
    \langle M_{\perp}\rangle=0
    \langle M_{\parallel}\rangle=F_m L(\eta)
    \langle M_{\perp}^2\rangle=F_m^2 L(\eta)/\eta
    \langle M_{\parallel}^2\rangle=F_m^2 (1-2L(\eta)/\eta)

with $F_m$ the magnetic scattering amplitude, the Langevin function $\langle \cos \alpha\rangle= L(\eta)=\coth(\eta)-1/\eta$ and $M_{\perp}$ and $M_\parallel$ the two transversal and the in-field magnetisation components, respectively (for a detailed discussion see [#Muehlbauer2019]_). The transversal magnetisation components of core, shell and solvent are often considered coaligned, corresponding to a coupling parameter $\delta =1$, and 0 if the transversal magnetisation of core and shell are uncorrelated (random orientational relation), i.e. the shell magnetisation component (perpendicular to the field) points in an arbitrary direction with respect to the core spin deflection. The scattering intensity reduces in this case to the sum over the individual squared form factors. The crossed terms $F_{M,i}\cdot F_{M,j}$ disappear as the scattering amplitudes $F_{M,i}$ and $F_{M,j}$ of two different parts can take both signs with equal propability.

The fraction of "upward" neutrons before ('up_frac_i') and after the sample ('up_frac_f') must range between 0 to 1, with 0.5 denoting an unpolarised beam. Note that a fit may result in a negative magnetic SLD, and hence magnetisation, when the polarisation state is inverted, i.e. if you have analysed for a $I_{00}$ state wheras your data are $I_{11}$. The model allows to construct the 4 spin-resolved cross sections (non-spin-flip I_{00}, I_{11} and spin-flip, here I_{01}=I_{10}), half-polarised SANS (SANSpol, incoming polarised beam I_0 and I_1, no analysis after sample 'up_frac_f'$=0.5$), and unpolarised beam ('up_frac_i'$=$'up_frac_f'$=0.5$).
Differences and other combinations between polarised scattering cross section, e.g. to obtain the nuclear-magnetic interference scattering, can be constructed with a custom model (Fitting> Add/Multiply Model) and using approbriate scales.

Every component in the core-multishell-matrix model has a (material dependent) saturation magnetisation and
a Langevin term $\eta$.
The Langevin term (i.e. mean orientation of the magnetic moment to the magnetic field) should be in the simplest case identical for core and shell, and related to the 
average particle moment $\mu=V*M_S$. However, it may vary due to different magnetic anisotropy, spin canting or other disturbances 
acting as a modified effective magnetic field from simple coaligned state.

For dense systems, where a structure factor of the arrangement might be relevant, special care has to be taken to describe 
the local environment around a particle and the structure factor is highly field-dependent due to varying correlations between the particles (see [#Honecker2020]_ for further details). The use of structure model is therefore strongly discouraged. A (micromagnetic) model is often the more approbriate approach.



Validation
----------

The model needs testing and validation.Your feedback is kindly requested.

Demonstrate that the model reproduces known results.

Test correct weighting of spin resolved cross sections to reconstruct unpolarised etc, also with imperfect neutron optics.

Test that 2d gives same result as 1d.



References
----------

.. [#Wiedenmann2005] A. Wiedenmann, *Physica B* 356, 246 (2005)
.. [#Muehlbauer2019] S. Muehlbauer, A. Michels et al., *Rev. Mod. Phys.* 91, 015004 (2019)
.. [#Honecker2020] D. Honecker, L. Fernandez Barguin, and P. Bender, *Phys. Rev. B* 101, 134401 (2020)



Authorship and Verification
----------------------------
* **Author: Dirk Honecker **Date:** September 04, 2020
* **Last Modified by:**
* **Last Reviewed by:**

"""

import numpy as np
from numpy import pi, inf

name = "magnetic_langevin_core_shell_sphere_3D"
title = "Dilute magnetic core-shell particles in a matrix with the magnetisation relaxing with respect to magnetic field"
description = """
    I(q) = A_ij (F_N^2(q)+ C_ij(L(H)) F_N F_M ) +B_ij(L(H)) F_M^2(q) 
            A_ij: weighting function =1 for unpolarised beam and non-neutron-spin-flip scattering, zero for spin-flip scattering
            B_ij(L(H)): weighting function for purely magnetic scattering, different for the various possible spin-resolved scattering cross sections
            C_ij(L(H)): weighting function for nuclear-magnetic interference scattering
            L(H): Langevin term describing the orientational distribution of magnetic moments with respect to magnetic field
            F_N: nuclear form factor
            F_M: magnetic form factor
The core and each shell can have a unique thickness, Langevin parameter and magnetic and nuclear sld.
"""
category = "shape:sphere"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["nuc_sld_core", "1e-6/Ang^2", 1.0, [-inf, inf], "", "Core scattering length density"],
              ["magnetic_sld_core",    "1e-6/Ang^2", 1.0,  [-inf, inf], "",    "Magnetic core scattering length density"],
              ["eta_core", "None", 3.0,  [0, inf], "",    "Langevin parameter of core"],
              ["radius", "Ang", 50., [0, inf], "volume", "Radius of the core"],
              ["nuc_sld_solvent", "1e-6/Ang^2", 6.4, [-inf, inf], "", "Solvent scattering length density"],
              ["magnetic_sld_solvent", "1e-6/Ang^2", 3.0,  [-inf, inf], "",    "Magnetic Solvent scattering length density"],
              ["eta_solvent", "None ", 3.0,  [0, inf], "",    "Langevin parameter of solvent of solvent"],
              ["delta_solvent", "None", 1.0,  [0, 1], "",    "Disorder coupling parameter of matrix to core"],
              ["n", "", 1, [0, 10], "volume", "Number of shells"],
              ["nuc_sld_shell[n]", "1e-6/Ang^2", 1.7, [-inf, inf], "", "Scattering length density of shell k"],
              ["magnetic_sld_shell[n]", "1e-6/Ang^2", 1.7, [-inf, inf], "", "Magnetic scattering length density of shell k"],
              ["eta[n]", "None", 3.0,  [0, inf], "",    "Langevin parameter of shell k"],
              ["delta[n]", "None", 1.0,  [0, 1], "",    "Disorder coupling parameter of shell to core"],
              ["thickness[n]", "Ang", 40., [0, inf], "volume", "Thickness of shell k"],
              ["up_i",          "None",  0, [0, 1], "", "Polarisation incoming beam"],
              ["up_f",          "None",  0, [0, 1], "", "Polarisation outgoing beam"],
              ["alpha",          "None",  0, [0, 180], "", "inclination of field to neutron beam"],
              ["beta",          "None",  0, [0, 360], "", "rotation of field around neutron beam"],             
             ]
# pylint: enable=bad-whitespace, line-too-long




source = ["lib/sas_3j1x_x.c","lib/gauss76.c", "magnetic_langevin_core_shell_3D.c"]
structure_factor = False
have_Fq = False
single=False


effective_radius_type = ["outer radius", "core radius"]

def random():
    """Return a random parameter set for the model."""
    num_shells = np.minimum(np.random.poisson(3)+1, 10)
    total_radius = 10**np.random.uniform(1.7, 4)
    thickness = np.random.exponential(size=num_shells+1)
    thickness *= total_radius/np.sum(thickness)
    pars = dict(
        #background=0,
        n=num_shells,
        radius=thickness[0],
    )
    for k, v in enumerate(thickness[1:]):
        pars['thickness%d'%(k+1)] = v
    return pars

def profile(nuc_sld_core, magnetic_sld_core,eta_core,radius, nuc_sld_solvent,magnetic_sld_solvent,eta_solvent,delta_solvent, n, nuc_sld_shell,magnetic_sld_shell,eta,delta, thickness,up_i,up_f, alpha,beta):
    """
    Returns the SLD profile *r* (Ang), and *rho* (1e-6/Ang^2).
    """
    n = int(n+0.5)
    z = []
    rho = []
    rho_mag = []

    # add in the core
    z.append(0)
    rho.append(nuc_sld_core)
    rho_mag.append(magnetic_sld_core)
    z.append(radius)
    rho.append(nuc_sld_core)
    rho_mag.append(magnetic_sld_core)

    # add in the shells
    for k in range(int(n)):
        # Left side of each shells
        z.append(z[-1])
        rho.append(nuc_sld_shell[k])
        rho_mag.append(magnetic_sld_shell[k])
        z.append(z[-1] + thickness[k])
        rho.append(nuc_sld_shell[k])
    # add in the solvent
    z.append(z[-1])
    rho.append(nuc_sld_solvent)
    rho_mag.append(magnetic_sld_solvent)
    z.append(z[-1]*1.25)
    rho.append(nuc_sld_solvent)
    rho_mag.append(magnetic_sld_solvent)

    return np.asarray(z), np.asarray(rho)
    







# tests = [
    # [{'radius': 20.0, 'thickness': 10.0}, 0.1, None, None, 30.0, 4.*pi/3*30**3, 1.0],

    # # The SasView test result was 0.00169, with a background of 0.001
    # [{'radius': 60.0, 'thickness': 10.0, 'sld_core': 1.0, 'sld_shell': 2.0,
      # 'sld_solvent': 3.0, 'background': 0.0}, 0.4, 0.000698838],
# ]
