r"""
Definition
----------
This model is a micromagnetic approach to analyse the SANS that arises from
nanoscale variations in the magnitude and orientation of the magnetization in
bulk ferromagnets in the approach to magnetic saturation (single domain state).
Typical materials are cold-worked elemental magnets, hard and soft magnetic 
nanocomposites, amorphous alloys and precipitates in magnetic steel [#Michels2014]_.
The magnetic SANS depends on the magnetic interactions, the magnetic microstructure
(defect/particle size, magnetocrystalline anisotropy, saturation magnetisation)
and on the applied magnetic field. As shown in [#Michels2016]_ near magnetic 
saturation the scattering cross-section can be evaluated by means of micromagnetic theory

.. math::
    I(\mathbf{Q}) = I_{nuc} + I_{mag}(\mathbf{Q},H),

with the field-independent nuclear and magnetic SANS cross section (due
to nanoscale spatial variations of the magnetisation).

.. math::
    I_{mag}(\mathbf{Q},H)= S_K(Q) R_K(\mathbf{Q}, H_i) + S_M(Q) R_M(\mathbf{Q}, H_i),

with $H_i$ the internal field, i.e. the external magnetic field corrected for
demagnetizing effects and the influence of the magnetodipolar field and of the
magnetic anisotropy [#Bick2013]_. This magnetic field dependence of the scattering 
reflects the increasing magnetisation misalignment with decreasing
externally applied magnetic field with a contribution $S_K \times R_K$ due to
perturbations around magnetic anisotropy fields and a term $S_M \times R_M$
related to magnetostatic fields. The magnetic moments decorate perturbations in the
microstructure (precipitates, grain boundaries etc).
The anisotropy-field function $S_K$ depends on the Fourier transform of the magnetic
anisotropy distribution (strength and orientation) in the material, and the
scattering function of the longitudinal magnetisation $S_M$ reflects the
variations of the saturation magnetisation, e.g. jumps at the particle-matrix
interface. $R_K$ and $R_M$ denote the micromagnetic response functions that
describe the magnetisation distribution around a perturbation in magnetic
anisotropy and flucutations in the saturation magnetisation value.

.. figure:: img/micromagnetic_FF.png

    Magnetisation distribution around (left) a particle with magnetic easy axis
    in the vertical direction and (right) a precipitation with a magnetisation
    that is higher than the matrix phase.  

The micromagnetic response functions depend on magnetic material parameters $M_S$:
average saturation magnetisation of the material, $H_i$: the internal magnetic
field, $A$ the average exchange-stiffness constant. In the vicinity of lattice 
imperfection in ferromagnetic materials, antisymmetric Dzyaloshinskii–Moriya
interaction (DMI) can occur due to the local structural inversion symmetry
breaking [#Arrott1963]_. DMI with strength $D$ can give rise to nonuniform spin
textures resulting in a polarization-dependent asymmetric scattering term for 
polycrystalline ferromagnetic with a centrosymmetric crystal structure [#Michels2016]_.
We assume (for simplicity) an isotropic microstructure (for $S_M$) and random 
orientation of magnetic easy axes (additionally for $S_K$) such that the 
contributions of the magnetic microstructure only depend on the magnitude of $q$.
Considerations for a microstructure with a prefered orientation (texture) can be
found in [#Weissmueller2001]_. In the code the averaging procedure over the random
anisotropy is explicitely performed. A specific orientation distribution can be
implemented by rewriting the model. 

The magnetic field is oriented with an inclination of $\alpha$ to the neutron beam
and rotated by $\beta$. The model for the nuclear scattering amplitude, saturation
magnetisation is based on spherical particles with a core shell structure. For
simplicity, only the core has an effective anisotropy, that is varying randomly
in direction from particle to particle. The effect of different, more complex 
spatial profiles of the anisotropy can be seen in [#Michels2010]_.
The magnetic scattering length density (SLD) is defined as 
$\rho_{\mathrm{mag}}=b_H M_S$, where $b_H= 2.91*10^{8}A^{-1}m^{-1}$ and $M_S$
is the saturation magnetisation (in $A/m$).

The fraction of "upward" neutrons before ('up_frac_i') and after the sample 
('up_frac_f') must range between 0 to 1, with 0.5 denoting an unpolarised beam.
Note that a fit may result in a negative magnetic SLD, and hence magnetisation, 
when the polarisation state is inverted, i.e. if you have analysed for a $I_{00}$
state wheras your data are $I_{11}$. The model allows to construct the 4 
spin-resolved cross sections (non-spin-flip $I_{00}$, $I_{11}$ and spin-flip, here
$I_{01}=I_{10}$), half-polarised SANS (SANSpol, incoming polarised beam $I_0$ and
$I_1$, no analysis after sample 'up_frac_f'$=0.5$), and unpolarised beam 
('up_frac_i'$=$'up_frac_f'$=0.5$). Differences and other combinations between
polarised scattering cross section, e.g. to obtain the nuclear-magnetic 
interference scattering, or subtraction of the residual scattering of the high
field reference state can be constructed with a custom model (Fitting>
Add/Multiply Model) and using approbriate scales. For dense systems, special
care has to be taken as the nculear structure factor (arrangement of particles)
does not need to be identical with the magnetic microstructure e.g. local 
textures and correlations between easy axes (see [#Honecker2020]_ for further
details). The use of structure model is therefore strongly discouraged. Better
$I_{nuc}$, $S_K$ and $S_M$ are fit independent from each other in a model-free way. 



References
----------

.. [#Arrott1963] A. Arrott, J. Appl. Phys. 34, 1108 (1963).
.. [#Weissmueller2001] J. Weissmueller et al., *Phys. Rev. B* 63, 214414 (2001).
.. [#Bick2013] J.-P. Bick et al., *Appl. Phys. Lett.* 102, 022415 (2013).
.. [#Michels2010] A. Michels et al., *Phys. Rev. B* 82, 024433 (2010).
.. [#Michels2014] A. Michels, *J. Phys.: Condens. Matter* 26, 383201 (2014).
.. [#Michels2016] A. Michels et al., *Phys. Rev. B* 94, 054424 (2016).
.. [#Honecker2020] D. Honecker, L. Fernandez Barguin, and P. Bender, *Phys. Rev. B* 101, 134401 (2020).



Authorship and Verification
----------------------------

* **Author:** Dirk Honecker **Date:** January 14, 2021
* **Last Modified by:** Dirk Honecker **Date:** September 23, 2024
* **Last Reviewed by:**

"""

import numpy as np
from numpy import pi, inf

name = "micromagnetic_FF_3D"
title = "Field-dependent magnetic microstructure around imperfections in bulk ferromagnets"
description = """
    I(q) = A (F_N^2(q)+ C F_N F_M + D F_M^2) +B(H) I_mag(q,H) 
            A: weighting function =1 for unpolarised beam and non-neutron-spin-flip scattering, zero for spin-flip scattering. The terms in the bracket are the residual scattering at perfectly saturating magnetic field.
            B(H): weighting function for purely magnetic scattering I_mag(q,H)  due to misaligned magnetic moments, different for the various possible spin-resolved scattering cross sections
            C: weighting function for nuclear-magnetic interference scattering
            F_N: nuclear form factor
            F_M: magnetic form factor
The underlying defect can have a core-shell structure.
"""
category = "shape:sphere"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["radius", "Ang", 50., [0, inf], "volume", "Structural radius of the core"],
              ["thickness", "Ang", 40., [0, inf], "volume", "Structural thickness of shell"],            
              ["nuc_sld_core", "1e-6/Ang^2", 1.0, [-inf, inf], "", "Core scattering length density"],
              ["nuc_sld_shell", "1e-6/Ang^2", 1.7, [-inf, inf], "", "Scattering length density of shell"],
              ["nuc_sld_solvent", "1e-6/Ang^2", 6.4, [-inf, inf], "", "Solvent scattering length density"],
              ["mag_sld_core",    "1e-6/Ang^2", 1.0,  [-inf, inf], "",    "Magnetic scattering length density of core"],
              ["mag_sld_shell", "1e-6/Ang^2", 1.7, [-inf, inf], "", "Magnetic scattering length density of shell"],
              ["mag_sld_solvent", "1e-6/Ang^2", 3.0,  [-inf, inf], "",  "Magnetic scattering length density of solvent"],
              ["hk_sld_core", "1e-6/Ang^2", 1.0,  [0, inf], "",    "Anisotropy field of defect"],
              ["Hi", "T", 2.0,  [0, inf], "",    "Effective field inside the material"],
              ["Ms", "T", 1.0,  [0, inf], "",    "Volume averaged saturation magnetisation"], 
              ["A", "pJ/m", 10.0,  [0, inf], "",    "Average exchange stiffness constant"],  
              ["D", "mJ/m^2", 0.0,  [0, inf], "",    "Average DMI constant"],                          
              ["up_i",          "None",  0.5, [0, 1], "", "Polarisation incoming beam"],
              ["up_f",          "None",  0.5, [0, 1], "", "Polarisation outgoing beam"],
              ["alpha",          "None",  90, [0, 180], "", "Inclination of field to neutron beam"],
              ["beta",          "None",  0, [0, 360], "", "Rotation of field around neutron beam"],             
             ]
# pylint: enable=bad-whitespace, line-too-long




source = ["lib/sas_3j1x_x.c", "lib/core_shell.c", "lib/gauss76.c", "lib/magnetic_functions.c", "micromagnetic_FF_3D.c"]
structure_factor = False
have_Fq = False
single=False
opencl = False


def random():
    """Return a random parameter set for the model."""
    outer_radius = 10**np.random.uniform(1.3, 4.3)
    # Use a distribution with a preference for thin shell or thin core
    # Avoid core,shell radii < 1
    radius = np.random.beta(0.5, 0.5)*(outer_radius-2) + 1
    thickness = outer_radius - radius
    pars = dict(
        radius=radius,
        thickness=thickness,
        
    )
    return pars



tests = [
     [{},1.002266990452620e-03, 7.461046163627724e+03],
     [{},(0.0688124,  -0.0261013),  22.024],
]
