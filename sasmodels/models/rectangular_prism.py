# rectangular_prism model
# Note: model title and parameter table are inserted automatically
r"""

This model provides the form factor, $P(q)$, for a rectangular prism.

Note that this model is almost totally equivalent to the existing
:ref:`parallelepiped` model.
The only difference is that the way the relevant
parameters are defined here ($a$, $b/a$, $c/a$ instead of $a$, $b$, $c$)
which allows use of polydispersity with this model while keeping the shape of
the prism (e.g. setting $b/a = 1$ and $c/a = 1$ and applying polydispersity
to *a* will generate a distribution of cubes of different sizes).
Note also that, contrary to :ref:`parallelepiped`, it does not compute
the 2D scattering.


Definition
----------

The 1D scattering intensity for this model was calculated by Mittelbach and
Porod (Mittelbach, 1961), but the implementation here is closer to the
equations given by Nayuk and Huber (Nayuk, 2012).
Note also that the angle definitions used in the code and the present
documentation correspond to those used in (Nayuk, 2012) (see Fig. 1 of
that reference), with $\theta$ corresponding to $\alpha$ in that paper,
and not to the usual convention used for example in the
:ref:`parallelepiped` model. As the present model does not compute
the 2D scattering, this has no further consequences.

In this model the scattering from a massive parallelepiped with an
orientation with respect to the scattering vector given by $\theta$
and $\phi$

.. math::

  A_P\,(q) =
      \frac{\sin \left( \tfrac{1}{2}qC \cos\theta \right) }{\tfrac{1}{2} qC \cos\theta}
      \,\times\,
      \frac{\sin \left( \tfrac{1}{2}qA \cos\theta \right) }{\tfrac{1}{2} qA \cos\theta}
      \,\times\ ,
      \frac{\sin \left( \tfrac{1}{2}qB \cos\theta \right) }{\tfrac{1}{2} qB \cos\theta}

where $A$, $B$ and $C$ are the sides of the parallelepiped and must fulfill
$A \le B \le C$, $\theta$ is the angle between the $z$ axis and the
longest axis of the parallelepiped $C$, and $\phi$ is the angle between the
scattering vector (lying in the $xy$ plane) and the $y$ axis.

The normalized form factor in 1D is obtained averaging over all possible
orientations

.. math::
  P(q) =  \frac{2}{\pi} \int_0^{\frac{\pi}{2}} \,
  \int_0^{\frac{\pi}{2}} A_P^2(q) \, \sin\theta \, d\theta \, d\phi

And the 1D scattering intensity is calculated as

.. math::
  I(q) = \text{scale} \times V \times (\rho_\text{p} -
  \rho_\text{solvent})^2 \times P(q)

where $V$ is the volume of the rectangular prism, $\rho_\text{p}$
is the scattering length of the parallelepiped, $\rho_\text{solvent}$
is the scattering length of the solvent, and (if the data are in absolute
units) *scale* represents the volume fraction (which is unitless).

**The 2D scattering intensity is not computed by this model.**


Validation
----------

Validation of the code was conducted by comparing the output of the 1D model
to the output of the existing :ref:`parallelepiped` model.


References
----------

P Mittelbach and G Porod, *Acta Physica Austriaca*, 14 (1961) 185-211

R Nayuk and K Huber, *Z. Phys. Chem.*, 226 (2012) 837-854

"""

from numpy import pi, inf, sqrt

name = "rectangular_prism"
title = "Rectangular parallelepiped with uniform scattering length density."
description = """
    I(q)= scale*V*(sld - sld_solvent)^2*P(q,theta,phi)+background
        P(q,theta,phi) = (2/pi) * double integral from 0 to pi/2 of ...
           AP^2(q)*sin(theta)*dtheta*dphi
        AP = S(q*C*cos(theta)/2) * S(q*A*sin(theta)*sin(phi)/2) * S(q*B*sin(theta)*cos(phi)/2)
        S(x) = sin(x)/x
"""
category = "shape:parallelepiped"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 6.3, [-inf, inf], "sld",
               "Parallelepiped scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["length_a", "Ang", 35, [0, inf], "volume",
               "Shorter side of the parallelepiped"],
              ["b2a_ratio", "", 1, [0, inf], "volume",
               "Ratio sides b/a"],
              ["c2a_ratio", "", 1, [0, inf], "volume",
               "Ratio sides c/a"],
             ]

source = ["lib/gauss76.c", "rectangular_prism.c"]

def ER(length_a, b2a_ratio, c2a_ratio):
    """
        Return equivalent radius (ER)
    """
    b_side = length_a * b2a_ratio
    c_side = length_a * c2a_ratio

    # surface average radius (rough approximation)
    surf_rad = sqrt(length_a * b_side / pi)

    ddd = 0.75 * surf_rad * (2 * surf_rad * c_side + (c_side + surf_rad) * (c_side + pi * surf_rad))
    return 0.5 * (ddd) ** (1. / 3.)

def random():
    import numpy as np
    a, b, c = 10**np.random.uniform(1, 4.7, size=3)
    pars = dict(
        length_a=a,
        b2a_ratio=b/a,
        c2a_ratio=c/a,
    )
    return pars

# parameters for demo
demo = dict(scale=1, background=0,
            sld=6.3, sld_solvent=1.0,
            length_a=35, b2a_ratio=1, c2a_ratio=1,
            length_a_pd=0.1, length_a_pd_n=10,
            b2a_ratio_pd=0.1, b2a_ratio_pd_n=1,
            c2a_ratio_pd=0.1, c2a_ratio_pd_n=1)

tests = [[{}, 0.2, 0.375248406825],
         [{}, [0.2], [0.375248406825]],
        ]
