#fcc paracrystal model
#note model title and parameter table are automatically inserted
#note - calculation requires double precision
r"""
Calculates the scattering from a **face-centered cubic lattice** with
paracrystalline distortion. Thermal vibrations are considered to be
negligible, and the size of the paracrystal is infinitely large.
Paracrystalline distortion is assumed to be isotropic and characterized by
a Gaussian distribution.

The returned value is scaled to units of |cm^-1|\ |sr^-1|, absolute scale.

Definition
----------

The scattering intensity *I(q)* is calculated as

.. image:: img/image158.jpg

where *scale* is the volume fraction of spheres, *Vp* is the volume of
the primary particle, *V(lattice)* is a volume correction for the crystal
structure, *P(q)* is the form factor of the sphere (normalized), and *Z(q)*
is the paracrystalline structure factor for a face-centered cubic structure.

Equation (1) of the 1990 reference is used to calculate *Z(q)*, using
equations (23)-(25) from the 1987 paper for *Z1*\ , *Z2*\ , and *Z3*\ .

The lattice correction (the occupied volume of the lattice) for a
face-centered cubic structure of particles of radius *R* and nearest
neighbor separation *D* is

.. image:: img/image159.jpg

The distortion factor (one standard deviation) of the paracrystal is
included in the calculation of *Z(q)*

.. image:: img/image160.jpg

where *g* is a fractional distortion based on the nearest neighbor distance.

The face-centered cubic lattice is

.. image:: img/image161.jpg

For a crystal, diffraction peaks appear at reduced q-values given by

.. image:: img/image162.jpg

where for a face-centered cubic lattice *h*\ , *k*\ , *l* all odd or all
even are allowed and reflections where *h*\ , *k*\ , *l* are mixed odd/even
are forbidden. Thus the peak positions correspond to (just the first 5)

.. image:: img/image163.jpg

**NB: The calculation of** *Z(q)* **is a double numerical integral that
must be carried out with a high density of** **points to properly capture
the sharp peaks of the paracrystalline scattering.** So be warned that the
calculation is SLOW. Go get some coffee. Fitting of any experimental data
must be resolution smeared for any meaningful fit. This makes a triple
integral. Very, very slow. Go get lunch!

This example dataset is produced using 200 data points, *qmin* = 0.01 |Ang^-1|,
*qmax* = 0.1 |Ang^-1| and the above default values.

.. image:: img/image164.jpg

*Figure. 1D plot in the linear scale using the default values (w/200 data point).*

The 2D (Anisotropic model) is based on the reference below where *I(q)* is
approximated for 1d scattering. Thus the scattering pattern for 2D may not
be accurate. Note that we are not responsible for any incorrectness of the
2D model computation.

.. image:: img/image165.gif

.. image:: img/image166.jpg

*Figure. 2D plot using the default values (w/200X200 pixels).*

REFERENCE

Hideki Matsuoka et. al. *Physical Review B*, 36 (1987) 1754-1765
(Original Paper)

Hideki Matsuoka et. al. *Physical Review B*, 41 (1990) 3854 -3856
(Corrections to FCC and BCC lattice structure calculation)
"""

from numpy import inf

name = "fcc_paracrystal"
title = "Face-centred cubic lattic with paracrystalline distortion"
description = """
    Calculates the scattering from a **face-centered cubic lattice** with paracrystalline distortion. Thermal vibrations
    are considered to be negligible, and the size of the paracrystal is infinitely large. Paracrystalline distortion is
    assumed to be isotropic and characterized by a Gaussian distribution.
    """
category = "shape:paracrystal"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["dnn", "Ang", 220, [-inf, inf], "", "Nearest neighbour distance"],
              ["d_factor", "", 0.06, [-inf, inf], "", "Paracrystal distortion factor"],
              ["radius", "Ang", 40, [0, inf], "volume", "Particle radius"],
              ["sld", "1e-6/Ang^2", 4, [-inf, inf], "", "Particle scattering length density"],
              ["solvent_sld", "1e-6/Ang^2", 1, [-inf, inf], "", "Solvent scattering length density"],
              ["theta", "degrees", 60, [-inf, inf], "orientation", "In plane angle"],
              ["phi", "degrees", 60, [-inf, inf], "orientation", "Out of plane angle"],
              ["psi", "degrees", 60, [-inf, inf], "orientation", "Out of plane angle"]
             ]

source = ["lib/J1.c", "lib/gauss150.c", "fcc.c"]

# parameters for demo
demo = dict(scale=1, background=0,
            dnn=220, d_factor=0.06, sld=4, solvent_sld=1,
            radius=40,
            theta=60, phi=60, psi=60,
            radius_pd=.2, radius_pd_n=0.2,
            theta_pd=15, theta_pd_n=0,
            phi_pd=15, phi_pd_n=0,
            psi_pd=15, psi_pd_n=0,
           )

# For testing against the old sasview models, include the converted parameter
# names and the target sasview model name.
oldname = 'FCCrystalModel'
oldpars = dict(sld='sldSph', solvent_sld='sldSolv')
