# precessing ellipsoid model
# Note: model title and parameter table are inserted automatically
r"""

Definition
----------

The model describes the scattering of a single-domain haematite nanospindles oriented in a magnetic field [#Zakutna2019]_.
The magnetic easy axis of haematite is fixed in the basal plane coinciding with the equatorial plane of the particle. 
This leads to the peculiar behaviour that the particle aligns with its short axis along the external field and the long axis is free to rotate around
the magnetic magnetic moment. In zero field the nanospindles in dispersion are randomly oriented, whereas for high magnetic field the particle long axis rotates in the plane perpendicular to the field due to Brownian motion. 

The 2D scattering intensity is modeled with the form factor of an ellipsoid 
given by [#Feigin1987]_

.. math::

    P(q,\alpha) = \frac{\text{scale}}{V} F^2(q,\alpha) + \text{background}

where

.. math::

    F(q,\alpha) = 3\Delta \rho V \frac{(\sin qr  - qr \cos qr)}{(qr)^3}

for

.. math::

    r = \left[ R_e^2 \sin^2 \alpha + R_p^2 \cos^2 \alpha \right]^{1/2}


$\alpha$ is the angle between the polar axis of the ellipsoid and $\mathbf{q}$, 
$V = (4/3)\pi R_pR_e^2$ is the volume of the ellipsoid, $R_p$ is the polar
radius along the rotational axis of the ellipsoid, $R_e$ is the equatorial
radius and
$\Delta \rho$ is the scattering length density difference between
the scatterer and the solvent.

The field-induced orientation will be derived for anisometric particles with the polar axis 
n.

The easy axis 

.. math::

    \mathbf{e}_A = \left( \cos \psi, \sin \psi \sin \gamma_1, \sin \psi \cos \gamma_1 \right)

is oriented with the polar angle $\psi$ towards the magnetic field and precessing under an angle $\gamma_1$.

The particle polar axis (perpendicular to the magnetic easy axis) is found with

.. math::

    \mathbf{n} = \mathbf{n}_0 \cos \gamma_2 + \mathbf{n}_0 \times \mathbf{e}_A \sin gamma_2

and an arbitrary orientation vector perpendicular to the easy axis

.. math::

    \mathbf{n}_0 = \left( -\sin \psi, \cos \psi \sin \gamma_1, \cos \psi \cos \gamma_1 \right)

The projection angle $\alpha$ is then obtained with

.. math::

    \cos \alpha = \mathbf{n} \cdot \mathbf{q}

The particle is free to precess over the angles $\gamma_1$ and $\gamma_2$. The distribution $P(\psi)$ of the easy axis to the magnetic field is given by Boltzmann statistics

.. math::

    P(\psi)=\xi \exp(\xi \cos \psi)/exp(\xi -1)
    
and the Langevin parameter $\xi=\mu \mu_0 H/ k_B T$ with $\mu_0$ the permeability of free space, $\mu$ the integral particle moment, $k_B$ the Boltzmann constant and $T$ the absolute temperature.
    
    
The angle $\beta$ allows to choose the direction of the applied magnetic field, with $\beta=90$ Degrees the field aligned along the horizontal axis, and for $\beta=0$ Degrees the field directed along the beam direction.


The model is computationally demanding with 3 averages over orientational distributions, make sure to have a performant system with GPU and patience.



Validation
----------



References
----------

.. [#Feigin1987] L. A. Feigin and D. I. Svergun. *Structure Analysis by Small-Angle X-Ray and Neutron Scattering*, Plenum Press, New York, 1987
.. [#Zakutna2019]  D. Zakutna, Y. Falke, D. Dresen, S. Prevost, P. Bender, D. Honecker and S. Disch, *Nanoscale*, 11 (2019), 7149-7156, DOI: 10.1039/c8nr09583c
 

Authorship and Verification
----------------------------

* **Converted to sasmodels by:** Dominique Dresen **Date:** March 25, 2019
* **Last Modified by:** Dirk Honecker **Date:** September 23, 2020"""

import numpy as np
from numpy import inf, sin, cos, pi

name = "morp_ellipsoid"
title = "Magnetically oriented, rotating and precessing (MORP) ellipsoid with uniform scattering length density."

description = """\

P(q.alpha)= scale*f(q)^2 + background, where f(q)= 3*(sld
        - sld_solvent)*V*[sin(q*r(Rp,Re,alpha))
        -q*r*cos(qr(Rp,Re,alpha))]
        /[qr(Rp,Re,alpha)]^3"

     r(Rp,Re,alpha)= [Re^(2)*(sin(alpha))^2
        + Rp^(2)*(cos(alpha))^2]^(1/2)

        sld: SLD of the ellipsoid
        sld_solvent: SLD of the solvent
        V: volume of the ellipsoid
        Rp: polar radius of the ellipsoid
        Re: equatorial radius of the ellipsoid
"""
category = "shape:ellipsoid"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 40, [-inf, inf], "sld",
               "Ellipsoid scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 8, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["radius_polar", "Ang", 1630, [0, inf], "volume",
               "Polar radius"],
              ["radius_equatorial", "Ang", 270, [0, inf], "volume",
               "Equatorial radius"],
              ["xi", "", 1, [0, inf], "",	
	               "Langevin parameter"]]             


source = ["lib/sas_3j1x_x.c", "lib/gauss76.c", "morp_ellipsoid.c"]
have_Fq = True
effective_radius_type = [
    "average curvature", "equivalent volume sphere", "min radius", "max radius",
    ]

def random():
    """Return a random parameter set for the model."""
    volume = 10**np.random.uniform(5, 12)
    radius_polar = 10**np.random.uniform(1.3, 4)
    radius_equatorial = np.sqrt(volume/radius_polar) # ignore 4/3 pi
    pars = dict(

        radius_polar=radius_polar,
        radius_equatorial=radius_equatorial,
    )
    return pars

demo = dict(scale=1, background=0,
            sld=40, sld_solvent=8,
            radius_polar=1630, radius_equatorial=270,

            radius_equatorial_pd=.2, radius_equatorial_pd_n=15,

            )
            
# Test 1-D and 2-D models            
# q = 0.1
# # april 6 2017, rkh add unit tests, NOT compared with any other calc method, assume correct!
# qx = q*cos(pi/6.0)
# qy = q*sin(pi/6.0)
# tests = [
#     [{}, 0.05, 54.8525847025],
#     [{'theta':80., 'phi':10.}, (qx, qy), 1.74134670026],
# ]
# del qx, qy  # not necessary to delete, but cleaner
