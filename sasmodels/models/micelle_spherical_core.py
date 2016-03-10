r"""

This model provides the form factor, $P(q)$, for a micelle with a spherical
core and Gaussian polymer chains attached to the surface.

Definition
----------

The 1D scattering intensity for this model is calculated according to
the equations given by Pedersen (Pedersen, 2000).

Validation
----------

This model has not yet been validated. Feb2015


Reference
---------

J Pedersen, *J. Appl. Cryst.*, 33 (2000) 637-640

"""

from numpy import inf

name = "micelle_spherical_core"
title = "Micelle Spherical Core model"
description = """
    This model provides the form factor, $P(q)$, for a micelle with
    a spherical core and Gaussian polymer chains attached to the surface.
    """
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["ndensity",      "1e15/cm^3",  8.94, [0.0, inf], "", "Number density of micelles"],
    ["v_core",        "Ang^3",  62624.0,  [0.0, inf], "", "Core volume"],
    ["v_corona",      "Ang^3",  61940.0,  [0.0, inf], "", "Corona volume"],
    ["solvent_sld",   "1e-6/Ang^2", 6.4,  [0.0, inf], "", "Solvent scattering length density"],
    ["core_sld",      "1e-6/Ang^2", 0.34, [0.0, inf], "", "Core scattering length density"],
    ["corona_sld",    "1e-6/Ang^2", 0.8,  [0.0, inf], "", "Corona scattering length density"],
    ["radius_core",   "Ang",       45.0,  [0.0, inf], "", "Radius of core"],
    ["radius_gyr",    "Ang",       20.0,  [0.0, inf], "", "Radius of gyration of chains in corona"],
    ["d_penetration", "",           1.0,  [-inf, inf], "", "Factor to mimic non-penetration of Gaussian chains"],
    ["n_aggreg",      "",           6.0,  [-inf, inf], "", "Aggregation number of the micelle"],
    ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sph_j1c.c", "micelle_spherical_core.c"]

demo = dict(scale=1, background=0,
            ndensity=8.94,
            v_core=62624.0,
            v_corona=61940.0,
            solvent_sld=6.4,
            core_sld=0.34,
            corona_sld=0.8,
            radius_core=45.0,
            radius_gyr=20.0,
            d_penetration=1.0,
            n_aggreg=6.0)


oldname = 'MicelleSphCoreModel'
oldpars = dict(solvent_sld='rho_solv',
               core_sld='rho_core',
               corona_sld='rho_corona')

tests = [
    [{}, 0.01, 15.3532],
    ]
