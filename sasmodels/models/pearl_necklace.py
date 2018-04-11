r"""
This model provides the form factor for a pearl necklace composed of two
elements: *N* pearls (homogeneous spheres of radius *R*) freely jointed by *M*
rods (like strings - with a total mass *Mw* = *M* \* *m*\ :sub:`r` + *N* \* *m*\
:sub:`s`, and the string segment length (or edge separation) *l*
(= *A* - 2\ *R*)). *A* is the center-to-center pearl separation distance.

.. figure:: img/pearl_necklace_geometry.jpg

    Pearl Necklace schematic

Definition
----------

The output of the scattering intensity function for the pearl_necklace is
given by (Schweins, 2004)

.. math::

    I(q)=\frac{ \text{scale} }{V} \cdot \frac{(S_{ss}(q)+S_{ff}(q)+S_{fs}(q))}
        {(M \cdot m_f + N \cdot m_s)^2} + \text{bkg}

where

.. math::

    S_{ss}(q) &= sm_s^2\psi^2(q)[\frac{N}{1-sin(qA)/qA}-\frac{N}{2}-
        \frac{1-(sin(qA)/qA)^N}{(1-sin(qA)/qA)^2}\cdot\frac{sin(qA)}{qA}] \\
    S_{ff}(q) &= sm_r^2[M\{2\Lambda(q)-(\frac{sin(ql/2)}{ql/2})\}+
        \frac{2M\beta^2(q)}{1-sin(qA)/qA}-2\beta^2(q)\cdot
        \frac{1-(sin(qA)/qA)^M}{(1-sin(qA)/qA)^2}] \\
    S_{fs}(q) &= m_r \beta (q) \cdot m_s \psi (q) \cdot 4[
        \frac{N-1}{1-sin(qA)/qA}-\frac{1-(sin(qA)/qA)^{N-1}}{(1-sin(qA)/qA)^2}
        \cdot \frac{sin(qA)}{qA}] \\
    \psi(q) &= 3 \cdot \frac{sin(qR)-(qR)\cdot cos(qR)}{(qR)^3} \\
    \Lambda(q) &= \frac{\int_0^{ql}\frac{sin(t)}{t}dt}{ql} \\
    \beta(q) &= \frac{\int_{qR}^{q(A-R)}\frac{sin(t)}{t}dt}{ql}

where the mass *m*\ :sub:`i` is (SLD\ :sub:`i` - SLD\ :sub:`solvent`) \*
(volume of the *N* pearls/rods). *V* is the total volume of the necklace.

The 2D scattering intensity is the same as $P(q)$ above, regardless of the
orientation of the *q* vector.

The returned value is scaled to units of |cm^-1| and the parameters of the
pearl_necklace model are the following

NB: *num_pearls* must be an integer.

References
----------

R Schweins and K Huber, *Particle Scattering Factor of Pearl Necklace Chains*,
*Macromol. Symp.* 211 (2004) 25-42 2004
"""

import numpy as np
from numpy import inf, pi

name = "pearl_necklace"
title = "Colloidal spheres chained together with no preferential orientation"
description = """
Calculate form factor for Pearl Necklace Model
[Macromol. Symp. 2004, 211, 25-42]
Parameters:
background:background
scale: scale factor
sld: the SLD of the pearl spheres
sld_string: the SLD of the strings
sld_solvent: the SLD of the solvent
num_pearls: number of the pearls
radius: the radius of a pearl
edge_sep: the length of string segment; surface to surface
thick_string: thickness (ie, diameter) of the string
"""
category = "shape:cylinder"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["radius", "Ang", 80.0, [0, inf], "volume",
               "Mean radius of the chained spheres"],
              ["edge_sep", "Ang", 350.0, [0, inf], "volume",
               "Mean separation of chained particles"],
              ["thick_string", "Ang", 2.5, [0, inf], "volume",
               "Thickness of the chain linkage"],
              ["num_pearls", "none", 3, [1, inf], "volume",
               "Number of pearls in the necklace (must be integer)"],
              ["sld", "1e-6/Ang^2", 1.0, [-inf, inf], "sld",
               "Scattering length density of the chained spheres"],
              ["sld_string", "1e-6/Ang^2", 1.0, [-inf, inf], "sld",
               "Scattering length density of the chain linkage"],
              ["sld_solvent", "1e-6/Ang^2", 6.3, [-inf, inf], "sld",
               "Scattering length density of the solvent"],
             ]

source = ["lib/sas_Si.c", "lib/sas_3j1x_x.c", "pearl_necklace.c"]
single = False  # use double precision unless told otherwise

def volume(radius, edge_sep, thick_string, num_pearls):
    """
    Calculates the total particle volume of the necklace.
    Redundant with form_volume.
    """
    num_pearls = int(num_pearls + 0.5)
    number_of_strings = num_pearls - 1.0
    string_vol = edge_sep * pi * pow((thick_string / 2.0), 2.0)
    pearl_vol = 4.0 /3.0 * pi * pow(radius, 3.0)
    total_vol = number_of_strings * string_vol
    total_vol += num_pearls * pearl_vol
    return total_vol

def ER(radius, edge_sep, thick_string, num_pearls):
    """
    Calculation for effective radius.
    """
    num_pearls = int(num_pearls + 0.5)
    tot_vol = volume(radius, edge_sep, thick_string, num_pearls)
    rad_out = (tot_vol/(4.0/3.0*pi)) ** (1./3.)
    return rad_out

def random():
    radius = 10**np.random.uniform(1, 3) # 1 - 1000
    thick_string = 10**np.random.uniform(0, np.log10(radius)-1) # 1 - radius/10
    edge_sep = 10**np.random.uniform(0, 3)  # 1 - 1000
    num_pearls = np.round(10**np.random.uniform(0.3, 3)) # 2 - 1000
    pars = dict(
        radius=radius,
        edge_sep=edge_sep,
        thick_string=thick_string,
        num_pearls=num_pearls,
    )
    return pars

# parameters for demo
demo = dict(scale=1, background=0, radius=80.0, edge_sep=350.0,
            num_pearls=3, sld=1, sld_solvent=6.3, sld_string=1,
            thick_string=2.5,
            radius_pd=.2, radius_pd_n=5,
            edge_sep_pd=25.0, edge_sep_pd_n=5,
            num_pearls_pd=0, num_pearls_pd_n=0,
            thick_string_pd=0.2, thick_string_pd_n=5,
           )

tests = [[{}, 0.001, 17380.245], [{}, 'ER', 115.39502]]
