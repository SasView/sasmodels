# Note: model title and parameter table are inserted automatically
r"""
This model calculates the scattering from a stack of repeating lamellar
structures. The stacks of lamellae (infinite in lateral dimension) are
treated as a paracrystal to account for the repeating spacing. The repeat
distance is further characterized by a Gaussian polydispersity. **This model
can be used for large multilamellar vesicles.**

Definition
----------

In the equations below,

- *scale* is used instead of the mass per area of the bilayer $\Gamma_m$
  (this corresponds to the volume fraction of the material in the bilayer,
  *not* the total excluded volume of the paracrystal),

- *sld* $-$ *solvent_sld* is the contrast $\Delta \rho$,

- *thickness* is the layer thickness $t$,

- *Nlayers* is the number of layers $N$,

- *spacing* is the average distance between adjacent layers
  $\langle D \rangle$, and

- *spacing_polydisp* is the relative standard deviation of the Gaussian
  layer distance distribution $\sigma_D / \langle D \rangle$.


The scattering intensity $I(q)$ is calculated as

.. math::

    I(q) = 2\pi\Delta\rho^2\Gamma_m\frac{P_\text{bil}(q)}{q^2} Z_N(q)

The form factor of the bilayer is approximated as the cross section of an
infinite, planar bilayer of thickness $t$

.. math::

    P_\text{bil}(q) = \left(\frac{\sin(qt/2)}{qt/2}\right)^2

$Z_N(q)$ describes the interference effects for aggregates
consisting of more than one bilayer. The equations used are (3-5)
from the Bergstrom reference:

.. math::


    Z_N(q) = \frac{1 - w^2}{1 + w^2 - 2w \cos(q \langle D \rangle)}
        + x_N S_N + (1 - x_N) S_{N+1}

where

.. math::

    S_N(q) = \frac{a_N}{N}[1 + w^2 - 2 w \cos(q \langle D \rangle)]^2

and

.. math::

    a_N &= 4w^2 - 2(w^3 + w) \cos(q \langle D \rangle) \\
        &\quad - 4w^{N+2}\cos(Nq \langle D \rangle)
        + 2 w^{N+3}\cos[(N-1)q \langle D \rangle]
        + 2w^{N+1}\cos[(N+1)q \langle D \rangle]

for the layer spacing distribution $w = \exp(-\sigma_D^2 q^2/2)$.

Non-integer numbers of stacks are calculated as a linear combination of
the lower and higher values

.. math::

    N_L = x_N N + (1 - x_N)(N+1)

The 2D scattering intensity is the same as 1D, regardless of the orientation
of the $q$ vector which is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


.. figure:: img/lamellarPC_1d.jpg

    1D plot using the default values above (w/20000 data point).

Reference
---------

M Bergstrom, J S Pedersen, P Schurtenberger, S U Egelhaaf,
*J. Phys. Chem. B*, 103 (1999) 9888-9897

"""

from numpy import inf

name = "lamellarPC"
title = "Random lamellar sheet with paracrystal structure factor"
description = """\
    [Random lamellar phase with paracrystal structure factor]
        randomly oriented stacks of infinite sheets
        with paracrytal S(Q), having polydisperse spacing.
        sld = sheet scattering length density
        sld_solvent = solvent scattering length density
        background = incoherent background
        scale = scale factor
"""
category = "shape:lamellae"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["thickness", "Ang", 33.0, [0, inf], "volume",
               "sheet thickness"],
              ["Nlayers", "", 20, [0, inf], "",
               "Number of layers"],
              ["spacing", "Ang", 250., [0.0, inf], "",
               "d-spacing of paracrystal stack"],
              ["spacing_polydisp", "Ang", 0.0, [0.0, inf], "",
               "d-spacing polydispersity"],
              ["sld", "1e-6/Ang^2", 1.0, [-inf, inf], "",
               "layer scattering length density"],
              ["solvent_sld", "1e-6/Ang^2", 6.34, [-inf, inf], "",
               "Solvent scattering length density"],
             ]


source = ["lamellarPC_kernel.c"]

form_volume = """
    return 1.0;
    """

Iqxy = """
    return Iq(sqrt(qx*qx+qy*qy), IQ_PARAMETERS);
    """

# ER defaults to 0.0
# VR defaults to 1.0

demo = dict(scale=1, background=0,
            thickness=33, Nlayers=20, spacing=250, spacing_polydisp=0.2,
            sld=1.0, solvent_sld=6.34,
            thickness_pd=0.2, thickness_pd_n=40)

oldname = 'LamellarPCrystalModel'
oldpars = dict(spacing_polydisp='pd_spacing', sld='sld_layer',
               solvent_sld='sld_solvent')
