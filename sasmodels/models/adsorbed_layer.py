#adsorbed_layer model
#conversion of Core2ndMomentModel.py
#converted by Steve King, Mar 2016


r"""
This model describes the scattering from a layer of surfactant or polymer adsorbed on spherical particles under the conditions that (i) the particles (cores) are contrast-matched to the dispersion medium, (ii) *S(Q)* ~ 1 (ie, the particle volume fraction is dilute), (iii) the particle radius is >> layer thickness (ie, the interface is locally flat), and (iv) scattering from excess unadsorbed adsorbate in the bulk medium is absent or has been corrected for.

Unlike many other core-shell models, this model does not assume any form for the density distribution of the adsorbed species normal to the interface (cf, a core-shell model normally assumes the density distribution to be a homogeneous step-function). For comparison, if the thickness of a (traditional core-shell like) step function distribution is *t*, the second moment about the mean of the density distribution (ie, the distance of the centre-of-mass of the distribution from the interface), |sigma| =  sqrt((*t* :sup:`2` )/12).

Definition
----------

.. math::

     I(q) = \text{scale} \cdot(\rho_\text{poly}-\rho_\text{solvent})^2    \left[\frac{6\pi\phi_\text{core}}{Q^2}\frac{\Gamma^2}{\delta_\text{poly}^2R_\text{core}} \exp(-Q^2\sigma^2)\right] + \text{background}

where *scale* is a scale factor, |rho|\ :sub:`poly` is the sld of the polymer (or surfactant) layer, |rho|\ :sub:`solv` is the sld of the solvent/medium and cores, |phi|\ :sub:`core` is the volume fraction of the core paraticles, |delta|\ :sub:`poly` is the bulk density of the polymer, |biggamma| is the adsorbed amount, and |sigma| is the second moment of the thickness distribution.

Note that all parameters except the |sigma| are correlated so fitting more than one of these parameters will generally fail. Also note that unlike other shape models, no volume normalization is applied to this model (the calculation is exact).

.. figure:: img/adsorbed_layer_1d.jpg

    1D plot using the default values.

References
----------

S King, P Griffiths, J. Hone, and T Cosgrove, *SANS from Adsorbed Polymer Layers*,
*Macromol. Symp.*, 190 (2002) 33-42.
"""

from numpy import inf, sqrt, pi, exp

name =  "adsorbed_layer"
title =  "Scattering from an adsorbed layer on particles"

description =  """
    Evaluates the scattering from particles
    with an adsorbed layer of surfactant or
    polymer, independent of the form of the
    density distribution.
    """
category =  "shape:sphere"

#             ["name", "units", default, [lower, upper], "type", "description"],
parameters =  [["second_moment", "Ang", 23.0, [0.0, inf], "", "Second moment"],
              ["adsorbed_amount", "mg/m2", 1.9, [0.0, inf], "", "Adsorbed amount"],
              ["density_poly", "g/cm3", 0.7, [0.0, inf], "", "Polymer density"],
              ["radius", "Ang", 500.0, [0.0, inf], "", "Particle radius"],
              ["vol_frac", "none", 0.14, [0.0, inf], "", "Particle vol fraction"],
              ["polymer_sld", "1/Ang^2", 1.5e-06, [-inf, inf], "", "Polymer SLD"],
              ["solvent_sld", "1/Ang^2", 6.3e-06, [-inf, inf], "", "Solvent SLD"]]

# NB: Scale and Background are implicit parameters on every model
def Iq(q, second_moment, adsorbed_amount, density_poly, radius, 
        vol_frac, polymer_sld, solvent_sld):
    # pylint: disable = missing-docstring
    deltarhosqrd =  (polymer_sld - solvent_sld) * (polymer_sld - solvent_sld)
    numerator =  6.0 * pi * vol_frac * (adsorbed_amount * adsorbed_amount)
    denominator =  (q * q) * (density_poly * density_poly) * radius
    eterm =  exp(-1.0 * (q * q) * (second_moment * second_moment))
    #scale by 10^10 for units conversion to cm^-1
    inten =  1.0e+10 * deltarhosqrd * ((numerator / denominator) * eterm)
    return inten * 9.4e-13
Iq.vectorized =  True  # Iq accepts an array of q values

def Iqxy(qx, qy, *args):
    # pylint: disable = missing-docstring
    return Iq(sqrt(qx ** 2 + qy ** 2), *args)
Iqxy.vectorized =  True # Iqxy accepts an array of qx, qy values

demo =  dict(scale = 1.0,
            second_moment = 23.0,
            adsorbed_amount = 1.9,
            density_poly = 0.7,
            radius = 500.0,
            vol_frac = 0.14,
            polymer_sld = 1.5e-06,
            solvent_sld = 6.3e-06,
            background = 0.0)

oldname =  "Core2ndMomentModel"
oldpars =  dict(scale = 'scale',
               second_moment = 'second_moment',
               adsorbed_amount = 'ads_amount',
               density_poly = 'density_poly',
               radius = 'radius_core',
               vol_frac = 'volf_cores',
               polymer_sld = 'sld_poly',
               solvent_sld = 'sld_solv',
               background = 'background')

tests =  [
    [{'scale': 1.0, 'second_moment': 23.0, 'adsorbed_amount': 1.9, 
     'density_poly': 0.7, 'radius': 500.0, 'vol_frac': 0.14, 
     'polymer_sld': 1.5e-06, 'solvent_sld': 6.3e-06, 'background': 0.0},
     [0.0106939, 0.469418], [73.741, 9.65391e-53]],
    ]
