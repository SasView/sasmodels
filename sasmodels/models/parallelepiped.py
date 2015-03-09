# parallelepiped model
# Note: model title and parameter table are inserted automatically
r"""
The form factor is normalized by the particle volume.

For information about polarised and magnetic scattering, click here_.

Definition
----------

This model provides the form factor, *P(q)*, for a rectangular parallelepiped
(below) where the form factor is normalized by the volume of the
parallelepiped. If you need to apply polydispersity, see also the
RectangularPrismModel_.

The calculated form factor is:

.. math::

    P(Q) = {\text{scale} \over V} F^2(Q) + \text{background}

where the volume *V* = *A B C* and the averaging < > is applied over all
orientations for 1D.

.. image:: img/parallelepiped.jpg

*Figure. Parallelepiped with the corresponding Definition of sides.

The edge of the solid must satisfy the condition that** *A* < *B* < *C*.
Then, assuming *a* = *A* / *B* < 1, *b* = *B* / *B* = 1, and
*c* = *C* / *B* > 1, the form factor is

.. math::

    P(q) = \frac{\textstyle{scale}}{V}\int_0^1 \phi(\mu \sqrt{1-\sigma^2},a)
    [S(\mu c \sigma/2)]^2 d\sigma

with

.. math::

    \phi(\mu,a) = \int_0^1 \{S[\frac{\mu}{2}\cos(\frac{\pi}{2}u)]
    S[\frac{\mu a}{2}\sin(\frac{\pi}{2}u)]\}^2 du

    S(x) = \frac{\sin x}{x}

    \mu = qB

and the contrast is defined as

.. math::

    \Delta\rho = \rho_{\textstyle p} - \rho_{\textstyle solvent}

The scattering intensity per unit volume is returned in units of |cm^-1|;
ie, *I(q)* = |phi| *P(q)*\ .

NB: The 2nd virial coefficient of the parallelpiped is calculated based on
the averaged effective radius (= sqrt(*short_a* \* *short_b* / |pi|)) and
length(= *long_c*) values, and used as the effective radius for
*S(Q)* when *P(Q)* \* *S(Q)* is applied.

To provide easy access to the orientation of the parallelepiped, we define
three angles |theta|, |phi| and |bigpsi|. The definition of |theta| and |phi|
is the same as for the cylinder model (see also figures below).
The angle |bigpsi| is the rotational angle around the *long_c* axis against
the *q* plane. For example, |bigpsi| = 0 when the *short_b* axis is parallel
to the *x*-axis of the detector.


.. _parallelepiped-orientation:

.. figure:: img/orientation.jpg

    Definition of the angles for oriented parallelepipeds.

.. figure:: img/orientation2.jpg

    Examples of the angles for oriented parallelepipeds against the detector plane.


Validation
----------

Validation of the code was done by comparing the output of the 1D calculation
to the angular average of the output of a 2D calculation over all possible
angles. The Figure below shows the comparison where the solid dot refers to
averaged 2D while the line represents the result of the 1D calculation (for
the averaging, 76, 180, 76 points are taken for the angles of |theta|, |phi|,
and |psi| respectively).

.. _parallelepiped-compare:

.. figure:: img/parallelepiped_compare.jpg

*Figure. Comparison between 1D and averaged 2D.*

This model reimplements the form factor calculations implemented in a c-library
provided by the NIST Center for Neutron Research (Kline, 2006).

"""

from numpy import pi, inf, sqrt

name = "parallelepiped"
title = "Rectangular parallelepiped with uniform scattering length density."
description = """
     P(q)= scale/V*integral from 0 to 1 of ...
           phi(mu*sqrt(1-sigma^2),a) * S(mu*c*sigma/2)^2 * dsigma

            phi(mu,a) = integral from 0 to 1 of ..
        (S((mu/2)*cos(pi*u/2))*S((mu*a/2)*sin(pi*u/2)))^2 * du
            S(x) = sin(x)/x
        mu = q*B
"""
category = "shape:parallelpiped"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 4, [-inf, inf], "",
               "Parallelepiped scattering length density"],
              ["solvent_sld", "1e-6/Ang^2", 1, [-inf, inf], "",
               "Solvent scattering length density"],
              ["a_side", "Ang", 35, [0, inf], "volume",
               "Shorter side of the parallelepiped"],
              ["b_side", "Ang", 75, [0, inf], "volume",
               "Second side of the parallelepiped"],
              ["c_side", "Ang", 400, [0, inf], "volume",
               "Larger side of the parallelepiped"],
              ["theta", "degrees", 60, [-inf, inf], "orientation",
               "In plane angle"],
              ["phi", "degrees", 60, [-inf, inf], "orientation",
               "Out of plane angle"],
              ["psi", "degrees", 60, [-inf, inf], "orientation",
               "Rotation angle around its own c axis against q plane"],
             ]

source = ["lib/J1.c", "lib/gauss76.c", "parallelepiped.c"]

def ER(a_side, b_side, c_side):

    # surface average radius (rough approximation)
    surf_rad = sqrt(a_side * b_side / pi)

    # DiamCyl recoded here (to check and possibly put in a library?)
    a = surf_rad
    b = 0.5 * c_side
    t1 = a * a * b
    t2 = 1.0 + (b / a) * (1.0 + a / b / 2.0) * (1.0 + pi * a / b / 2.0)
    ddd = 3.0 * t1 * t2

    return 0.5 * (ddd) ** (1. / 3.)

# parameters for demo
demo = dict(scale=1, background=0,
            sld=6.3e-6, solvent_sld=1.0e-6,
            a_side=35, b_side=75, c_side=400,
            theta=45, phi=30, psi=15,
            a_side_pd=0.1, a_side_pd_n=10,
            b_side_pd=0.1, b_side_pd_n=1,
            c_side_pd=0.1, c_side_pd_n=10,
            theta_pd=10, theta_pd_n=1,
            phi_pd=10, phi_pd_n=1,
            psi_pd=10, psi_pd_n=10)

# For testing against the old sasview models, include the converted parameter
# names and the target sasview model name.
oldname = 'ParallelepipedModel'
oldpars = dict(theta='parallel_theta', phi='parallel_phi', psi='parallel_psi',
               a_side='short_a', b_side='short_b', c_side='long_c',
               sld='sldPipe', solvent_sld='sldSolv')

