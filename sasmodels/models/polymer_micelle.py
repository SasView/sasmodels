r"""

This model provides the form factor, $P(q)$, for a micelle with a spherical
core and Gaussian polymer chains attached to the surface, thus may be applied
to block copolymer micelles. To work well the Gaussian chains must be much
smaller than the core, which is often not the case.  Please study the
reference carefully.

Definition
----------

The 1D scattering intensity for this model is calculated according to
the equations given by Pedersen (Pedersen, 2000).

Validation
----------

This model has not yet been validated. Feb2015


References
----------

J Pedersen, *J. Appl. Cryst.*, 33 (2000) 637-640

"""

from numpy import inf

name = "polymer_micelle"
title = "Polymer micelle model"
description = """
    This model provides an approximate form factor, P(q), for a micelle with
    a spherical core with Gaussian polymer chains attached to the surface.
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
