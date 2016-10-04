r"""
This model calculates the scattering from fractal-like aggregates based
on the Mildner reference.

Definition
----------

The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = scale \times P(q)S(q) + background

.. math::

    P(q) = F(qR)^2

.. math::

    F(x) = \frac{3\left[sin(x)-xcos(x)\right]}{x^3}

.. math::

    S(q) = \frac{\Gamma(5-D_S)\zeta^{5-D_S}}{\left[1+(q\zeta)^2
    \right]^{(5-D_S)/2}}
    \frac{sin\left[(D_S - 5) tan^{-1}(q\zeta) \right]}{q}

.. math::

    scale = scale\_factor \times NV^2(\rho_{particle} - \rho_{solvent})^2

.. math::

    V = \frac{4}{3}\pi R^3

where $R$ is the radius of the building block, $D_S$ is the **surface** fractal
dimension,| \zeta\|  is the cut-off length, $\rho_{solvent}$ is the scattering
length density of the solvent,
and $\rho_{particle}$ is the scattering length density of particles.

.. note::
    The surface fractal dimension $D_s$ is only valid if $1<surface\_dim<3$.
    It is also only valid over a limited $q$ range (see the reference for
    details)


References
----------

D Mildner and P Hall, *J. Phys. D: Appl. Phys.*, 19 (1986) 1535-1545

"""

from numpy import inf

name = "surface_fractal"
title = "Fractal-like aggregates based on the Mildner reference"
description = """\
    [The scattering intensity  I(x) = scale*P(x)*S(x) + background, where
        scale = scale_factor  * V * delta^(2)
        p(x) = F(x*radius)^(2)
        F(x) = 3*[sin(x)-x cos(x)]/x**3
        S(x) = [(gamma(5-Ds)*colength^(5-Ds)*[1+(x^2*colength^2)]^((Ds-5)/2)
             * sin[(Ds-5)*arctan(x*colength)])/x]
        where
        delta        =  sldParticle -sldSolv.
        radius       =  Particle radius
        fractal_dim_surf  =  Surface fractal dimension (Ds)
        co_length    =  Cut-off length
        background   =  background

        Ref.   :Mildner, Hall,J Phys D Appl Phys(1986), 19, 1535-1545
        Note I : This model is valid for 1<fractal_dim_surf<3 with limited q range.
        Note II: This model is not in absolute scale.
"""
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["radius",        "Ang", 10.0, [0, inf],   "",
               "Particle radius"],
              ["fractal_dim_surf",   "",    2.0,  [1, 3],   "",
               "Surface fractal dimension"],
              ["cutoff_length", "Ang", 500., [0.0, inf], "",
               "Cut-off Length"],
             ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/sph_j1c.c", "lib/sas_gamma.c", "surface_fractal.c"]

demo = dict(scale=1, background=0,
            radius=10, fractal_dim_surf=2.0, cutoff_length=500)

tests = [
    # Accuracy tests based on content in test/utest_other_models.py
    [{'radius': 10.0,
      'fractal_dim_surf': 2.0,
      'cutoff_length': 500.0,
     }, 0.05, 301428.66016],

    # Additional tests with larger range of parameters
    [{'radius': 1.0,
      'fractal_dim_surf': 1.0,
      'cutoff_length': 10.0,
     }, 0.332070182643, 1125.00421004],

    [{'radius': 3.5,
      'fractal_dim_surf': 0.1,
      'cutoff_length': 30.0,
      'background': 0.01,
     }, 5.0, 0.00999998891322],

    [{'radius': 3.0,
      'fractal_dim_surf': 1.0,
      'cutoff_length': 33.0,
      'scale': 0.1,
     }, 0.51, 2.50120147004],
    ]
