r"""

This model provides the form factor, $P(q)$, for a micelle with a spherical
core and Gaussian polymer chains attached to the surface, thus may be applied
to block copolymer micelles. To work well the Gaussian chains must be much
smaller than the core, which is often not the case.  Please study the
reference carefully.

Definition
----------

The 1D scattering intensity for this model is calculated according to
the equations given by Pedersen (Pedersen, 2000), summarised briefly here.

The micelle core is imagined as $N\_aggreg$ polymer heads, each of volume $v\_core$,
which then defines a micelle core of $radius\_core$, which is a separate parameter
even though it could be directly determined.
The Gaussian random coil tails, of gyration radius $rg$, are imagined uniformly 
distributed around the spherical core, centred at a distance $radius\_core + d\_penetration.rg$
from the micelle centre, where $d\_penetration$ is of order unity.
A volume $v\_corona$ is defined for each coil.
The model in detail seems to separately parametrise the terms for the shape of I(Q) and the
relative intensity of each term, so use with caution and check parameters for consistency.
The spherical core is monodisperse, so it's intensity and the cross terms may have sharp
oscillations (use q resolution smearing if needs be to help remove them).

.. math::
    P(q) = N^2\beta^2_s\Phi(qR)^2+N\beta^2_cP_c(q)+2N^2\beta_s\beta_cS_{sc}s_c(q)+N(N-1)\beta_c^2S_{cc}(q)
    
    \beta_s = v\_core(sld\_core - sld\_solvent)
    
    \beta_c = v\_corona(sld\_corona - sld\_solvent)

where $N = n\_aggreg$, and for the spherical core of radius $R$ 

.. math::   
   \Phi(qR)= \frac{\sin(qr) - qr\cos(qr)}{(qr)^3}

whilst for the Gaussian coils

.. math::

   P_c(q) &= 2 [\exp(-Z) + Z - 1] / Z^2

   Z &= (q R_g)^2

The sphere to coil ( core to corona) and coil to coil (corona to corona) cross terms are
approximated by:

.. math::
   
   S_{sc}(q)=\Phi(qR)\psi(Z)\frac{sin(q(R+d.R_g))}{q(R+d.R_g)}
   
   S_{cc}(q)=\psi(Z)^2\left[\frac{sin(q(R+d.R_g))}{q(R+d.R_g)} \right ]^2
   
   \psi(Z)=\frac{[1-exp^{-Z}]}{Z}

Validation
----------

$P(q)$ above is multiplied by $ndensity$, and a units conversion of 10^{-13}, so $scale$
is likely 1.0 if the scattering data is in absolute units. This model has not yet been 
independently validated.


References
----------

J Pedersen, *J. Appl. Cryst.*, 33 (2000) 637-640

"""

from numpy import inf

name = "polymer_micelle"
title = "Polymer micelle model"
description = """
This model provides the form factor, $P(q)$, for a micelle with a spherical
core and Gaussian polymer chains attached to the surface, thus may be applied
to block copolymer micelles. To work well the Gaussian chains must be much
smaller than the core, which is often not the case.  Please study the
reference to Pedersen and full documentation carefully. 
    """


category = "shape:sphere"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["ndensity",      "1e15/cm^3",  8.94, [0.0, inf], "", "Number density of micelles"],
    ["v_core",        "Ang^3",  62624.0,  [0.0, inf], "", "Core volume "],
    ["v_corona",      "Ang^3",  61940.0,  [0.0, inf], "", "Corona volume"],
    ["sld_solvent",   "1e-6/Ang^2", 6.4,  [0.0, inf], "sld", "Solvent scattering length density"],
    ["sld_core",      "1e-6/Ang^2", 0.34, [0.0, inf], "sld", "Core scattering length density"],
    ["sld_corona",    "1e-6/Ang^2", 0.8,  [0.0, inf], "sld", "Corona scattering length density"],
    ["radius_core",   "Ang",       45.0,  [0.0, inf], "", "Radius of core ( must be >> rg )"],
    ["rg",    "Ang",       20.0,  [0.0, inf], "", "Radius of gyration of chains in corona"],
    ["d_penetration", "",           1.0,  [-inf, inf], "", "Factor to mimic non-penetration of Gaussian chains"],
    ["n_aggreg",      "",           6.0,  [-inf, inf], "", "Aggregation number of the micelle"],
    ]
# pylint: enable=bad-whitespace, line-too-long

single = False

source = ["lib/sph_j1c.c", "polymer_micelle.c"]

demo = dict(scale=1, background=0,
            ndensity=8.94,
            v_core=62624.0,
            v_corona=61940.0,
            sld_solvent=6.4,
            sld_core=0.34,
            sld_corona=0.8,
            radius_core=45.0,
            rg=20.0,
            d_penetration=1.0,
            n_aggreg=6.0)


tests = [
    [{}, 0.01, 15.3532],
    ]
# RKH 20Mar2016 - need to check whether the core & corona volumes are per
#                 monomer ??? and how aggregation number works!
# renamed from micelle_spherical_core to polymer_micelle,
# moved from shape-independent to spheres section.
# Ought to be able to add polydisp to core? And add ability to x by S(Q) ?
