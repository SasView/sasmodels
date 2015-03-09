# Note: model title and parameter table are inserted automatically
r"""
This model provides the scattering intensity, $I(q) = P(q) S(q)$, for a
lamellar phase where a random distribution in solution are assumed.
Here a Caille $S(Q)$ is used for the lamellar stacks.

The scattering intensity $I(q)$ is

.. math:

    I(q) = 2\pi \frac{P(q)S(q)}{\delta q^2}

The form factor is

.. math:

    P(q) = \frac{2\Delta\rho^2}{q^2}\left(1-\cos q\delta \right)

and the structure factor is

.. math:

    S(q) = 1 + 2 \sum_1^{N-1}\left(1-\frac{n}{N}\right)
           \cos(qdn)\exp\left(-\frac{2q^2d^2\alpha(n)}{2}\right)

where

.. math:

    \begin{eqnarray}
    \alpha(n) &=& \frac{\eta_{cp}}{4\pi^2} \left(\ln(\pi n)+\gamma_E\right)  \\
    \gamma_E &=& 0.5772156649 && \text{Euler's constant} \\
    \eta_{cp} &=& \frac{q_o^2k_B T}{8\pi\sqrt{K\overline{B}}} && \text{Caille constant}
    \end{eqnarray}

Here $d$ = (repeat) spacing, $\delta$ = bilayer thickness,
the contrast $\Delta\rho$ = SLD(headgroup) - SLD(solvent),
$K$ = smectic bending elasticity, $B$ = compression modulus, and
$N$ = number of lamellar plates (*n_plates*).

NB: **When the Caille parameter is greater than approximately 0.8 to 1.0, the
assumptions of the model are incorrect.** And due to a complication of the
model function, users are responsible for making sure that all the assumptions
are handled accurately (see the original reference below for more details).

Non-integer numbers of stacks are calculated as a linear combination of
results for the next lower and higher values.

The 2D scattering intensity is calculated in the same way as 1D, where the
$q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}

The returned value is in units of |cm^-1|, on absolute scale.

.. image:: img/lamellarCaille_1d.jpg

*Figure. 1D plot using the default values (w/6000 data point).*

Our model uses the form factor calculations as implemented in a c library
provided by the NIST Center for Neutron Research (Kline, 2006).

REFERENCE
---------

F Nallet, R Laversanne, and D Roux, J. Phys. II France, 3, (1993) 487-502

also in J. Phys. Chem. B, 105, (2001) 11081-11088
"""
from numpy import inf

name = "lamellarPS"
title = "Random lamellar sheet with Caille structure factor"
description = """\
	[Random lamellar phase with Caille  structure factor]
        randomly oriented stacks of infinite sheets
		with Caille S(Q), having polydisperse spacing.
    	sld = sheet scattering length density
		sld_solvent = solvent scattering length density
		background = incoherent background
		scale = scale factor
"""
category = "shape:lamellae"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["thickness", "Ang",  30.0, [0, inf], "volume", "sheet thickness"],
              ["Nlayers", "",  20, [0, inf], "", "Number of layers"],
              ["spacing", "Ang", 400., [0.0,inf], "volume", "d-spacing of Caille S(Q)"],
              ["Caille_parameter", "1/Ang^2", 0.1, [0.0,0.8], "", "Caille parameter"],
              ["sld", "1e-6/Ang^2", 6.3, [-inf,inf], "",
               "layer scattering length density"],
              ["solvent_sld", "1e-6/Ang^2", 1.0, [-inf,inf], "",
               "Solvent scattering length density"],
             ]

source = ["lamellarCaille_kernel.c"]

# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return 1.0;
    """

Iqxy = """
    return Iq(sqrt(qx*qx+qy*qy), IQ_PARAMETERS);
    """

# ER defaults to 0.0
# VR defaults to 1.0

demo = dict(scale=1, background=0,
            thickness=67.,Nlayers=3.75,spacing=200.,
            Caille_parameter=0.268,sld=1.0, solvent_sld=6.34,
            thickness_pd= 0.1, thickness_pd_n=100,
            spacing_pd= 0.05, spacing_pd_n=40)

oldname = 'LamellarPSModel'
oldpars = dict(thickness='delta', Nlayers='N_plates', Caille_parameter='caille',
               sld='sld_bi',solvent_sld='sld_sol')
