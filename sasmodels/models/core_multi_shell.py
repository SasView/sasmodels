r"""
Definition
----------

This model is a trivial extension of the CoreShell function to a larger number
of shells. The scattering length density profile for the default sld values 
(w/ 4 shells).

.. figure:: img/core_multi_shell_sld_default_profile.jpg

    SLD profile of the core_multi_shell object from the center of sphere out
    for the default SLDs.*

The 2D scattering intensity is the same as $P(q)$ above, regardless of the
orientation of the $q$ vector which is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}

.. note:: **Be careful!** The SLDs and scale can be highly correlated. Hold as
         many of these parameters fixed as possible.

.. note:: The outer most radius (= *radius* + *thickness*) is used as the
          effective radius for $S(Q)$ when $P(Q)*S(Q)$ is applied.

For information about polarised and magnetic scattering, see 
the :doc:`magnetic help <../sasgui/perspectives/fitting/mag_help>` documentation.

Our model uses the form factor calculations implemented in a c-library provided
by the NIST Center for Neutron Research (Kline, 2006).

References
----------
See the :ref:`core_shell_sphere <core_shell_sphere>` model documentation.

L A Feigin and D I Svergun,
*Structure Analysis by Small-Angle X-Ray and Neutron Scattering*,
Plenum Press, New York, 1987.

**Author:** NIST IGOR/DANSE **on:** pre 2010

**Last Modified by:** in progress **on:** March 20, 2016

**Last Reviewed by:** in progress **on:** March 20, 2016
"""



from __future__ import division

import numpy as np
from numpy import inf, nan
from math import fabs, exp, expm1

name = "core_multi_shell"
title = "This model provides the scattering from a spherical core with 1 to 4 \
 concentric shell structures. The SLDs of the core and each shell are \
 individually specified."

description = """\
Form factor for a core muti-shell (up to 4) sphere normalized by the volume.
Each shell can have a unique thickness and sld.

	background:background,
	rad_core0: radius of sphere(core)
	thick_shell#:the thickness of the shell#
	sld_core0: the SLD of the sphere
	sld_solv: the SLD of the solvent
	sld_shell: the SLD of the shell#
	A_shell#: the coefficient in the exponential function
	
	
    scale: 1.0 if data is on absolute scale
    volfraction: volume fraction of spheres
    radius: the radius of the core
    sld: the SLD of the core
    thick_shelli: the thickness of the i'th shell from the core
    sld_shelli: the SLD of the i'th shell from the core
    sld_solvent: the SLD of the solvent
    background: incoherent background

"""

category = "shape:sphere"


#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld_core", "1e-6/Ang^2", 1.0, [-inf, inf], "",
               "Core scattering length density"],
              ["radius", "Ang", 200., [0, inf], "volume",
               "Radius of the core"],
              ["sld_solvent", "1e-6/Ang^2", 6.4, [-inf, inf], "",
               "Solvent scattering length density"],
              ["n", "", 1, [0, 10], "volume",
               "number of shells"],
              ["sld[n]", "1e-6/Ang^2", 1.7, [-inf, inf], "",
               "scattering length density of shell k"],
              ["thickness[n]", "Ang", 40., [0, inf], "volume",
               "Thickness of shell k"],
              ]

source = ["lib/sph_j1c.c", "core_multi_shell.c"]

def profile(sld_core, radius, sld_solvent, n, sld, thickness):
    """
    SLD profile
    
    :return: (r, beta) where r is a list of radius of the transition points\
      and beta is a list of the corresponding SLD values.

    """
    total_radius = 1.25*(sum(thickness[:n]) + radius + 1)

    r = []
    beta = []

    # add in the core
    r.append(0)
    beta.append(sld)
    r.append(radius)
    beta.append(sld)

    # add in the shells
    for k in range(n):
        # Left side of each shells
        r0 = r[-1]
        r.append(r0)
        beta.append(sld[k])
        r.append(r0 + thickness[k])
        beta.append(sld[k])
    # add in the solvent
    r.append(r[-1])
    beta.append(sld_solvent)
    r.append(total_radius)
    beta.append(sld_solvent)

    return np.asarray(r), np.asarray(beta)

def ER(radius, n, thickness):
    n = n[0]  # n cannot be polydisperse
    return np.sum(thickness[:n], axis=0) + radius

def VR(radius, n, thickness):
    return 1.0, 1.0

demo = dict(sld_core = 6.4,
            radius = 60,
            sld_solvent = 6.4,
            n = 2,
            sld = [2.0, 3.0],
            thickness = 20,
            thickness1_pd = 0.3,
            thickness2_pd = 0.3,
            thickness1_pd_n = 10,
            thickness2_pd_n = 10,
            )
