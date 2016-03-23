r"""
Definition
----------

This model is a trivial extension of the CoreShell function to a larger number
of shells. The scattering length density profile for the default sld values 
(w/ 4 shells).

.. figure:: img/core_multi_shell_sld_default_profile.jpg

    SLD profile of the core_multi_shell object from the center of sphere out
    for the default SLDs.*

The 2D scattering intensity is the same as *P(q)* above, regardless of the
orientation of the *q* vector which is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}

.. note:: **Be careful!** The SLDs and scale can be highly correlated. Hold as
         many of these parameters fixed as possible.

.. note:: The outer most radius (= *radius* + *thickness*) is used as the
          effective radius for *S(Q)* when *P(Q)* \* *S(Q)* is applied.

For information about polarised and magnetic scattering, see 
the :ref:`mag_help` documentation.

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
parameters = [["volfraction", "", 0.05, [0,1],"",
               "volume fraction of spheres"],
              ["sld", "1e-6/Ang^2", 1.0, [-inf, inf], "",
               "Core scattering length density"],
              ["radius", "Ang", 200., [0, inf], "",
               "Radius of the core"],
              ["sld_solvent", "1e-6/Ang^2", 6.4, [-inf, inf], "",
               "Solvent scattering length density"],
              ["n", "", 1, [0, 10], "volume",
               "number of shells"],
              ["sld_shell[n]", "1e-6/Ang^2", 1.7, [-inf, inf], "",
               "scattering length density of shell k"],
              ["thick_shell[n]", "Ang", 40., [0, inf], "volume",
               "Thickness of shell k"],
              ]

#source = ["lib/sph_j1c.c", "onion.c"]

def Iq(q, *args, **kw):
    return q

def Iqxy(qx, *args, **kw):
    return qx


def shape(core_sld, core_radius, solvent_sld, n, in_sld, out_sld, thickness, A):
    """
    SLD profile
    
    :return: (r, beta) where r is a list of radius of the transition points\
      and beta is a list of the corresponding SLD values.

    """
#        r = []
#        beta = []
#        # for core at r=0
#        r.append(0)
#        beta.append(self.params['sld_core0'])
#        # for core at r=rad_core
#        r.append(self.params['rad_core0'])
#        beta.append(self.params['sld_core0'])
#
#        # for shells
#        for n in range(1, self.n_shells+1):
#            # Left side of each shells
#            r0 = r[len(r)-1]
#            r.append(r0)
#            exec "beta.append(self.params['sld_shell%s'% str(n)])"
#
#            # Right side of each shells
#            exec "r0 += self.params['thick_shell%s'% str(n)]"
#            r.append(r0)
#            exec "beta.append(self.params['sld_shell%s'% str(n)])"
#
#        # for solvent
#        r0 = r[len(r)-1]            
#        r.append(r0)
#        beta.append(self.params['sld_solv'])
#        r_solv = 5*r0/4
#        r.append(r_solv)
#        beta.append(self.params['sld_solv'])
#
#        return r, beta
# above is directly from old code -- below is alotered from Woitek's first
# cut an the onion.  Seems like we should be consistant?

    total_radius = 1.25*(sum(thickness[:n]) + core_radius + 1)
    dr = total_radius/400  # 400 points for a smooth plot

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
        beta.append(sld_shell[k])
        r.append(r0 + thickness[k])
        beta.append(sld_shell[k])
   # add in the solvent
    r.append(r[-1])
    beta.append(solvent_sld)
    r.append(total_radius)
    beta.append(solvent_sld)

    return np.asarray(r), np.asarray(beta)

def ER(radius, n, thick_shell):
    return np.sum(thick_shell[:n], axis=0) + core_radius

def VR(radius, n, thick_shell):
    return 1.0

parameters = [["volfraction", "", 0.05, [0,1],"",
               "volume fraction of spheres"],
              ["sld", "1e-6/Ang^2", 1.0, [-inf, inf], "",
               "Core scattering length density"],
              ["radius", "Ang", 200., [0, inf], "",
               "Radius of the core"],
              ["sld_solvent", "1e-6/Ang^2", 6.4, [-inf, inf], "",
               "Solvent scattering length density"],
              ["n", "", 1, [0, 10], "volume",
               "number of shells"],
              ["sld_shell[n]", "1e-6/Ang^2", 1.7, [-inf, inf], "",
               "scattering length density of shell k"],
              ["thick_shell[n]", "Ang", 40., [0, inf], "volume",
               "Thickness of shell k"],
               ]

demo = dict(volfraction = 1.0,
            sld = 6.4,
            radius = 60,
            sld_solvent = 6.4,
            n = 1,
            sld_shell = [2.0],
            thick_shell = [10])

oldname = "CoreMultiShellModel"
oldpars = dict(
    sld="sld_core0",
    radius="rad_core0",
    sld_solvent="sld_solv",
    n="n_shells",
    sld_shell="sld_in_shell",
    thick_shell="thick_shell",
    # func_shell is always 2 in the user interface, so isn't included
    )
