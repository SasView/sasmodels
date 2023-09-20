r"""
Definition
----------

Core multi-shell cylinders or discs.

There must be a minumum of ONE shell (May set its sld to match solvent, 
thickness1 and face1 to zero or very small so the shell does not contribute 
to the normalisation volume.)

There may be numerical integration issues at extremes of Q and/or with 
extreme aspect ratios - particles 
which are very large in some dimensions compared to others.

2d scattering is so far only minimally tested.

Use of S(Q) with this I(Q) currently not giving correct I(Q) in sasview, 
(though passes unit tests), due to a more general sasview v5 bug.

The sld profile plot show profiles along both radius and half length 
simultaneously! (A simple edit of the py code will change which is 
displayed.)

Scattering is normalised to $V_{total} = \pi.Router^2Louter$ the outer 
volume of the particle.

$Louter = length\_core + 2*(face1 + face2 + face3 + ...)$

$Router = radius\_core + thickness1 + thickness2 + thickness3 + ...$

.. math::

Scattered intensity is calculated by a numerical integration over angle 
$\alpha$ between the axis of the cylinder and $\vec q$,

.. math::

    I(q,\alpha) = \frac{\text{scale}}{V_{total}} F^2(q,\alpha).sin(\alpha) 
                   + \text{background}

where

.. math::

    F(q,\alpha) =  \sum_{k=0}^{n} \left[ (\rho_k - \rho_{k+1}) V_k
           \frac{\sin \left( q.h_k\cos\alpha \right)}
                {q.h_k\cos\alpha}
           \frac{2 J_1 \left( q.r_k\sin\alpha \right)}
                {q.r_k\sin\alpha} 
           \right]
           exp\left \{ -\frac{1}{2}q^2\sigma^2 \right \}

and

    $\rho_0 = sld\_core$

    $\rho_k = sld\_shell_k\quad\text{for  k = 1 to n, the number of shells}$

    $\rho_{n+1} = sld\_solvent$

    $r_0 = radius\_core$

    $r_1 = r_0 + thickness_1$

    $r_k = r_{k-1} + thickness_k\quad\text{for  k = 1 to n}$

    $h_0 = \tfrac12 length\_core$

    $h_1 = h_0 + face_1$

    $h_k = h_{k-1} + face_k\quad\text{for  k = 1 to n}$

    $V_k = 2\pi.h_k. R_k^2$

    $J_1$ is the first order Bessel function.


.. _core-multi_shell-cylinder-geometry:

.. figure:: img/core_multi_shell_cylinder_geometry.png

    Core multi_shell cylinder or disc - schematic cross section with two shells.

An approximation for the effects of "Gaussian interfacial roughness" $\sigma$
is included, by multiplying $I(Q)$ by
$\exp\left \{ -\frac{1}{2}Q^2\sigma^2 \right \}$. This applies, in some way, to
all interfaces in the model not just the external ones. (Note that for a one
dimensional system convolution of the scattering length density profile with
a Gaussian of standard deviation $\sigma$ does exactly this multiplication.)
Leave $\sigma$ set to zero for the usual sharp interfaces.
There is some debate as to whether the factor of 1/2 is needed or not.

Since the number of parameters may become too large to fit well, fix or 
constrain as many as possible. e.g. to constrain face1 to thickness1, 
click on the first, ctrl/click the second, right click and select 
"mutual constraint". Note that such constraints are NOT applied inside
integrations over polydispersity, they are only applied to starting values.
To do this the constraints would need to be included inside a customized model.


To provide easy access to the orientation of the core-shell cylinder, we
define the axis of the cylinder using two angles $\theta$ and $\phi$.
(see :ref:`cylinder model <cylinder-angle-definition>`)


The $\theta$ and $\phi$ parameters are not used for the 1D output.

Reference
---------

See also Livsey [#Livsey1987]_ and Onsager [#Onsager1949]_.

.. [#Livsey1987] I Livsey, *J. Chem. Soc., Faraday Trans. 2*, 83 (1987) 1445-1452

.. [#Kline2006] S R Kline, *J Appl. Cryst.*, 39 (2006) 895

.. [#Onsager1949] L. Onsager, *Ann. New York Acad. Sci.*, 51 (1949) 627-659

Authorship and Verification
----------------------------

* **Author:** Richard Heenan **Date:** April 2021
* **Last Modified by:** Richard Heenan **Date:** April 2021
* **Last Reviewed by:** Richard Heenan **Date:** April 2021
"""

from __future__ import division

import numpy as np
from numpy import inf, nan, pi, sin, cos

name = "core_multi_shell_cylinder"
title = "core multi-shell cylinder or disc"

description = """\
Core multi_shell cylinder or disc"""

category = "shape:cylinder"

# TODO: n is a volume parameter that is not polydisperse

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["sld_core", "1e-6/Ang^2", 5.0, [-inf, inf], "sld", "Core scattering length density"],
    ["radius_core", "Ang", 35., [0, inf], "volume", "Radius of the core"],
    ["length_core", "Ang", 80., [0, inf], "volume", "Radius of the core"],
    ["sld_solvent", "1e-6/Ang^2", 6.4, [-inf, inf], "sld", "Solvent scattering length density"],
    ["sigma",       "Ang",        0,    [0, inf],    "",            "Interfacial roughness"],
    ["n_shells", "", 1, [1, 10], "volume", "number of shells (must be integer)"],
    ["sld_shell[n_shells]", "1e-6/Ang^2", 1.0, [-inf, inf], "sld", "scattering length density of shell k"],
    ["thickness[n_shells]", "Ang", 20., [0, inf], "volume", "Thickness of cylinder shell k"],
    ["face[n_shells]", "Ang", 10., [0, inf], "volume", "Thickness of disc (or cylinder end) shell k"],
 # orientation parameters 
    ["theta", "degrees", 60, [-360, 360], "orientation",
                             "cylinder axis to beam angle"],
    ["phi", "degrees", 60, [-360, 360], "orientation",
                            "rotation about beam"],
                            
    ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c", "core_multi_shell_cylinder.c"]
have_Fq = True
radius_effective_modes = [
    "excluded volume", "equivalent volume sphere", "outer radius", "half outer length",
    "half min outer dimension", "half max outer dimension", "half outer diagonal",
    ]
    
# TODO cpoied code as from core_multi_shell, WHY are the sld's not being set????
def random():
    """Return a random parameter set for the model."""
    num_shells = np.minimum(np.random.poisson(3)+1, 10)
    total_radius = 10**np.random.uniform(1.7, 4)
    thickness = np.random.exponential(size=num_shells+1)
    thickness *= total_radius/np.sum(thickness)
    
    total_length = 10**np.random.uniform(1.7, 4)
    face = np.random.exponential(size=num_shells+1)
    face *= total_length/np.sum(thickness)
    pars = dict(
        #background=0,
        n=num_shells,
        radius_core=thickness[0],
        length_core=face[0],
    )
    for k, v in enumerate(thickness[1:]):
        pars['thickness%d'%(k+1)] = v
    for k, v in enumerate(face[1:]):
        pars['face%d'%(k+1)] = v
    return pars

profile_axes = ['Radius (A)', 'SLD (1e-6/A^2)']
# This is called from sasmodels\sasmodels\sasview_model.py around line 485
# TODO have superimposed axial & radial profiles here, really need a second 
#   curve in different colour
def profile(sld_core, radius_core, length_core, sld_solvent, sigma, n_shells,
            sld_shell, thickness, face, theta, phi):
    """
    Needs to pass ALL the params through, else get an eror on building docs etc
    Returns shape profile with x=radius, y=SLD.
    Do radius profile for now, can we prodcue a 2nd plot for the face profile ?
    """
    n_shells = int(n_shells+0.5)
    z = []
    rho = []
    z_next = 0
    # radius profile,  two sld points for core
    z.append(z_next)
    rho.append(sld_core)
    z_next = radius_core
    z.append(z_next)
    rho.append(sld_core)

    for i in range(0, n_shells):
        z.append(z_next)
        rho.append(sld_shell[i])
        z_next += thickness[i]
        z.append(z_next)
        rho.append(sld_shell[i])
    z.append(z_next)
    rho.append(sld_solvent)
    z.append(z_next*1.5)
    rho.append(sld_solvent)
    
    # now the profile in length
    z2 = []
    # create new list of rho values
    rho2 = list(rho)
    z_next = 0
    # two sld points for core
    z2.append(z_next)
    z_next = 0.5*length_core
    z2.append(z_next)

    for i in range(0, n_shells):
        z2.append(z_next)
        z_next += face[i]
        z2.append(z_next)
    z2.append(z_next)
    z2.append(z_next*1.5)
    # now reverse it then append the radial list, note both profiles start
    # at (zero,sld_core), so get a continuous but possibly confusing trace
    # on the plot!
    z2.reverse()
    rho2.reverse()
    for i in range(0,len(z)):
        z2.append(z[i])
        rho2.append(rho[i])
    # return just z and rho for radial plot only, or comment out the final
    # section above, below "now reverse" for the original (z2, rho2) half 
    # length plot
    return np.asarray(z2), np.asarray(rho2)

# unit tests, need tests for at least 1, 2 and 3 shells to check all
# routes through the c++ code.
# Some cases of cylinder, core_shell_cylinder and a 2 shell approximation from 
# core_shell_bicelle_elliptical_belt_rough worked OK in this new model.
# TODO: 2d calc has NOT been checked !
q = 0.1
qx = q*cos(pi/6.0)
qy = q*sin(pi/6.0)
qxten, qyten = 0.2 * np.cos(2.5), 0.2 * np.sin(2.5)
theta, phi = 80.1534480601659, 10.1510817110481  # (10, 10) in sasview 3.x

# numerical values for I(Q) here not checked independently, but similar test 
# cases were working for 1d data.
tests = [
    # 
    [{"sld_core": 4.0, "radius_core": 10.0, "length_core": 100.0, 
      "sld_solvent": 6.0, "sigma": 0.0, "n_shells": 1,
      "sld_shell": [8.0],
      "thickness": [30.0],
      "face": [5.0]
     }, [0.001,0.01,0.05],[173.64805993875723, 162.92094280325867, 30.73988178927162]],
    [{"sld_core": 4.0, "radius_core": 10.0, "length_core": 100.0, 
      "sld_solvent": 6.0, "sigma": 0.0, "n_shells": 2,
      "sld_shell": [8.0, 0.5],
      "thickness": [30.0, 5.0],
      "face": [5.0, 10.0]
     }, [0.001,0.01,0.05],[33.55611252965908, 26.607752107659834, 14.388355760013978]],
    [{"sld_core": 4.0, "radius_core": 10.0, "length_core": 100.0, 
      "sld_solvent": 6.0, "sigma": 0.0, "n_shells": 3,
      "sld_shell": [8.0, 0.5,-0.5],
      "thickness": [30.0, 5.0,7.5],
      "face": [5.0, 10.0, 3.5]
     }, [0.001,0.01,0.05],[689.6953333069359, 590.4682661780487, 14.469981695284062]],
     # 2d test using default params
     [{}, (qx, qy), 0.39415111587961554],
     #
     # tests copied from cylinder.py, where R=20 L=400, del_rho =3
     # set up 2 shell cylinder to look like uniform cylinder
     # can only do the monodisperse cases.
    [{"sld_core": 3.0, "radius_core": 7.0, "length_core": 250.0, 
      "sld_solvent": 6.0, "sigma": 0.0, "n_shells": 2,
      "sld_shell": [3.0, 3.0],
      "thickness": [5.0, 8.0],
      "face": [50.0, 25.0]
     }, [0.2], [0.042761386790780453]],
    [{"scale": 1., "background": 0.,"sld_core": 3.0, "radius_core": 7.0, 
      "length_core": 250.0, 
      "sld_solvent": 6.0, "sigma": 0.0, "n_shells": 2,
      "sld_shell": [3.0, 3.0],
      "thickness": [5.0, 8.0],
      "face": [50.0, 25.0]}, [0.01, 0.05, 0.2],
     # these numerical results from cylinder NOT independently verified
     [3.01823887e+02, 5.36661653e+01, 4.17613868e-02]],
    [{"scale": 1., "background": 0.,"sld_core": 3.0, "radius_core": 7.0, 
      "length_core": 250.0, 
      "sld_solvent": 6.0, "sigma": 0.0, "n_shells": 2,
      "sld_shell": [3.0, 3.0],
      "thickness": [5.0, 8.0],
      "face": [50.0, 25.0]},
     # the longer list here checks  F1, F2=I(Q)*V, R_eff, volume, volume_ratio
     # F1 changed sign compared to cylinder - does that matter ?? - suspect not.
     0.05, -2214.9614083046904, 26975556.88749548, 73.34013315261608,
     502654.8245743669, 1.0],
     # 2d test from uniform cylinder
    [{"scale": 1., "background": 0.001,"sld_core": 3.0, "radius_core": 7.0, 
      "length_core": 250.0, 
      "sld_solvent": 6.0, "sigma": 0.0, "n_shells": 2,
      "sld_shell": [3.0, 3.0],
      "thickness": [5.0, 8.0],
      "face": [50.0, 25.0],
      'theta': theta, 'phi': phi}, [(qxten, qyten)], [0.03514647218513852]],
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
    [{"@S": "hardsphere",         # MONODISPERSE
      "scale": 5., "background": 0., "volfraction": 0.2,
      "structure_factor_mode": 0,  # normal decoupling approx
      "radius_effective_mode": 1,  # Reff "excluded volume"
      "sld_core": 3.0, "radius_core": 7.0, 
      "length_core": 250.0, 
      "sld_solvent": 6.0, "sigma": 0.0, "n_shells": 2,
      "sld_shell": [3.0, 3.0],
      "thickness": [5.0, 8.0],
      "face": [50.0, 25.0]
     }, [0.01, 0.05, 0.2], [7.35571916e+01, 5.78147797e+01, 4.15623248e-2]
     ],
    [{"@S": "hardsphere",
      "scale": 5., "background": 0., "volfraction": 0.2,
      "structure_factor_mode": 1,  # beta approx
      "radius_effective_mode": 1,  # Reff "excluded volume"
      "sld_core": 3.0, "radius_core": 7.0, 
      "length_core": 250.0, 
      "sld_solvent": 6.0, "sigma": 0.0, "n_shells": 2,
      "sld_shell": [3.0, 3.0],
      "thickness": [5.0, 8.0],
      "face": [50.0, 25.0]
     }, [0.01, 0.05, 0.2], [8.29729770e+01, 5.44206752e+01, 4.17598382e-2]
     ],
     ]
     
del qx, qy, qxten, qyten  # not necessary to delete, but cleaner
