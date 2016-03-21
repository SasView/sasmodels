r"""

A number of natural and commercial processes form high-surface area materials
as a result of the vapour-phase aggregation of primary particles.
Examples of such materials include soots, aerosols, and fume or pyrogenic
silicas. These are all characterised by cluster mass distributions (sometimes
also cluster size distributions) and internal surfaces that are fractal in
nature. The scattering from such materials displays two distinct breaks in
log-log representation, corresponding to the radius-of-gyration of the primary
particles, $rg$, and the radius-of-gyration of the clusters (aggregates),
$Rg$. Between these boundaries the scattering follows a power law related to
the mass fractal dimension, $Dm$, whilst above the high-Q boundary the
scattering follows a power law related to the surface fractal dimension of
the primary particles, $Ds$.

Definition
----------

The scattered intensity I(q) is calculated using a modified
Ornstein-Zernicke equation

.. math::

    I(q) = scale \times P(q) + background

    P(q) = \left\{ \left[ 1+(q^2a)\right]^{D_m/2} \times
                   \left[ 1+(q^2b)\right]^{(6-D_s-D_m)/2}
           \right\}^{-1}

    a = R_{g}^2/(3D_m/2)

    b = r_{g}^2/[-3(D_s+D_m-6)/2]

    scale = scale\_factor \times NV^2 (\rho_{particle} - \rho_{solvent})^2

where $R_g$ is the size of the cluster, $r_g$ is the size of the primary
particle, $D_s$ is the surface fractal dimension, $D_m$ is the mass fractal
dimension, $\rho_{solvent}$ is the scattering length density of the solvent,
and $\rho_{particle}$ is the scattering length density of particles.

.. note::

    The surface ( $D_s$ ) and mass ( $D_m$ ) fractal dimensions are only
    valid if $0 < surface\_dim < 6$ , $0 < mass\_dim < 6$ , and
    $(surface\_dim + mass\_dim ) < 6$ .


Reference
---------

P Schmidt, *J Appl. Cryst.*, 24 (1991) 414-435 Equation(19)

A J Hurd, D W Schaefer, J E Martin, *Phys. Rev. A*,
35 (1987) 2361-2364 Equation(2)

"""

from numpy import inf

name = "mass_surface_fractal"
title = "Mass Surface Fractal model"
description = """
        The scattering intensity  I(x) = scale*P(x)*S(x) + background, where
        p(x)= {[1+(x^2*a)]^(Dm/2) * [1+(x^2*b)]^(6-Ds-Dm)/2}^(-1)
        a = Rg^2/(3*Dm/2)
        b = rg^2/(3*(6-Ds-Dm)/2)
        scale        =  scale factor * N*Volume^2*contrast^2
        mass_dim       =  Dm (mass fractal dimension)
        surface_dim  =  Ds
        cluster_rg  =  Rg
        primary_rg    =  rg
        background   =  background
        Ref: Schmidt, J Appl Cryst, eq(19), (1991), 24, 414-435
        Hurd, Schaefer, Martin, Phys Rev A, eq(2),(1987),35, 2361-2364
        Note that 0 < Ds< 6 and 0 < Dm < 6.
        """
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["mass_dim",      "",    1.8, [1e-16, 6.0], "",
               "Mass fractal dimension"],
              ["surface_dim",   "",    2.3, [1e-16, 6.0], "",
               "Surface fractal dimension"],
              ["cluster_rg", "Ang",   86.7, [0.0, inf], "",
               "Cluster radius of gyration"],
              ["primary_rg", "Ang", 4000.,  [0.0, inf], "",
               "Primary particle radius of gyration"],
             ]
# pylint: enable=bad-whitespace, line-too-long

source = ["mass_surface_fractal.c"]

demo = dict(scale=1, background=0,
            mass_dim=1.8,
            surface_dim=2.3,
            cluster_rg=86.7,
            primary_rg=4000.0)

oldname = 'MassSurfaceFractal'
oldpars = dict(radius='radius',
               mass_dim='mass_dim',
               surface_dim='surface_dim',
               cluster_rg='cluster_rg',
               primary_rg='primary_rg')

tests = [

    # Accuracy tests based on content in test/utest_other_models.py
    [{'mass_dim':      1.8,
      'surface_dim':   2.3,
      'cluster_rg':   86.7,
      'primary_rg': 4000.0,
      'background':    0.0,
     }, 0.05, 1.77537e-05],

    # Additional tests with larger range of parameters
    [{'mass_dim':      3.3,
      'surface_dim':   1.0,
      'cluster_rg':   90.0,
      'primary_rg': 4000.0,
     }, 0.001, 0.18562699016],

    [{'mass_dim':      1.3,
      'surface_dim':   1.0,
      'cluster_rg':   90.0,
      'primary_rg': 2000.0,
      'background':    0.8,
     }, 0.001, 1.16539753641],

    [{'mass_dim':      2.3,
      'surface_dim':   1.0,
      'cluster_rg':   90.0,
      'primary_rg': 1000.0,
      'scale':        10.0,
      'background':    0.0,
     }, 0.051, 0.000169548800377],
    ]
