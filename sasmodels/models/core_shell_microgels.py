r"""
Polydisperse spherical core-shell microgel particles with fuzzy interfaces.
====================

This model provides the form factor, $P(q)$, for a spherical particle with
core-shell structure and fuzzy interfaces, as follows. A Lorentzian peak is
used to account for intra-particle composition fluctuations.

.. math::
    P(q) = \bigl[(\rho _{shell} - \rho _{solv})
    \Phi _{{shell}_{out}}(q, R_{{shell}_{out}}, s_{{shell}_{out}})
    + (\rho _{core} - \rho _{solv})\Phi _{core}(q, R_{core}, s_{core})
    - (\rho _{shell} - \rho _{solv})\Phi _{{shell}_{in}}(q, R_{{shell}_{in}}, s_{{shell}_{in}})\bigr]^2
    + \left(\frac{L_0}{1 + q^2\xi _{intra}^2}\right)

where $\Phi _k$ is the radial density profile of component $k$, given by:

.. math::
\begin{align*}
  \Phi _k(q,R_k,s_k) = \frac{4\pi }{q^4s_k^2} \biggl((R_k+s_k)cos[q(R_k+s_k)]
  + (R_k-s_k)cos[q(R_k-s_k)] - 3sin[q(R_k+s_k)]&  \\
  - 3sin[q(R_k-s_k)] + \frac{2cos(qR_k)}{q} + \frac{6sin(qR_k)}{q}&\biggr)
\end{align*}

and $k = core, shell_{in}, shell_{out}$.

The core-shell interface can be either symmetric or asymmetric, and arises
from interpenetration of the core polymer and the shell polymer. The model
assumes a box density profile for each component, with a parabolic gradation
at the interfaces. The model also includes a Lorentzian peak to account for
composition fluctuations within the particle. A radial density profile
$\Phi _{overall}(R)$ may be constructed from the fitted SANS curves. Polydispersity
is accounted for in the overall particle radius, $R_tot$ (i.e. core, shell,
inner interfaces, and outer surface/interface).

Add Figure 1 from Angew. paper

Parameters:\\
$\rho _{core}$: scattering length density of the core\\
$\rho _{shell}$: scattering length density of the shell\\
$\rho _{solv}$: scattering length density of the solvent\\
$R_{tot}$: radius of the full core-shell particle (polydisperse)\\
$R_{core}$: half-width radius of the core interface\\
$R_{{shell}_{in}}$: half-width radius of the shell inner interface\\
$s_{core}$: half of the thickness of the core interface\\
$s_{{shell}_{in}}$: half of the thickness of the shell inner interface\\
$s_{{shell}_{out}}$: half of the thickness of the shell outer interface\\
$L_0$: Lorentzian scale factor\\
$\xi _{intra}$: characteristic length of intra-particle fluctuations\\


Relationships to other physical parameters:\\
$W_{shell} = R_{tot} -  R_{{shell}_{in}} - s_{{shell}_{in}} - 2*s_{{shell}_{out}}$ = width/thickness of the shell \\
$W_{core} = R_{core} - s_{core}$ = width/radius of the core \\
$R_{{shell}_{out}} = R_{tot} - s_{{shell}_{out}}$ = half-width radius of the shell outer interface\\

Additionally, one can calculate the following volumes:\\
.. math::
\begin{align*}
  V_{core} &= 4\pi \left(\frac{R_{core}^3}{3} +
  \frac{R_{core}s_{core}^2}{6}\right) \\
  V_{shell} &= 4\pi \left[\frac{R_{{shell}_{out}}^3 - R_{{shell}_{in}}^3}{3} +
  \frac{R_{{shell}_{out}}s_{{shell}_{out}}^2 - R_{{shell}_{in}}s_{{shell}_{in}}^2}{6}\right]
\end{align*}

Special case for a symmetric core-shell interface:\\
$R_{{shell}_{in}} = R_{core}$ and $s_{{shell}_{in}} = s_{core}$

====================
References
-----------

Berndt, I., Pedersen, J. S., Richtering, W. J. Am. Chem. Soc. 2005, 127,
pp. 9372 - 9373.

Berndt, I., Pedersen, J. S., Lindner, P., Richtering, W. Langmuir. 2006, 22,
pp. 459 - 468.

Berndt, I., Pedersen, J. S., Richtering, W. Angew. Chem. 2006, 118,
pp. 1769 - 1773.

====================
Author: Rachel Ford

"""

import numpy as np

#name = """core_shell_microgels"""
#title = """Form factor for fuzzy spherical core-shell microgels
#    with a Lorentzian peak for intra-particle fluctuations."""
#description = """
#    P(q) = [(rho_shell - rho_solv)*Phi_shell_out(q, R_shell_out, s_shell_out)
#    + (rho_core - rho_solv)*Phi_core(q, R_core, s_core)
#    - (rho_shell - rho_solv)*Phi_shell_in(q, R_shell_in, s_shell_in)]**2
#    + L_0/(1 + q**2*xi_intra**2)
#
#    where density profile Phi(q, R, s) is given by:
#    Phi_k(q,R_k,s_k) = (4*pi/q**4*s_k**2)*((R_k+s_k)cos[q(R_k+s_k)] +
#    (R_k-s_k)cos[q(R_k-s_k)] - 3sin[q(R_k+s_k)] - 3sin[q(R_k-s_k)] +
#    2cos(qR_k)/q + 6sin(qR_k)/q)
#
#    and k = core, shell_in, shell_out.
#
#    List of parameters:
#    rho_core: SLD of the core
#    rho_shell: SLD of the shell
#    rho_solv: SLD of the solvent
#    R_tot: radius of the full particle (polydisperse)
#    R_core: half-width radius of the core
#    R_shell_in: half-width radius of the shell inner interface
#    s_core: half the thickness of the core interface
#    s_shell_in: half the thickness of the shell inner interface
#    s_shell_out: half of the thickness of the outer shell interface
#    L_0: Lorentzian scale factor
#    xi_intra: characteristic length of intra-particle fluctuations
#    """
#category = "shape:sphere"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
#parameters = [["rho_core",    "1e-06 Ang^-2",   2.0, [-inf, inf], "",       "SLD of core"],
#              ["rho_shell",   "1e-06 Ang^-2",   3.0, [-inf, inf], "",       "SLD of shell"],
#              ["rho_solv",    "1e-06 Ang^-2",   1.0, [-inf, inf], "",       "SLD of solvent"],
#              ["R_tot",       "Ang",          700.0,    [0, inf], "volume", "Radius of full particle"],
#              ["R_core",      "Ang",          500.0,    [0, inf], "",       "Radius of core"],
#              ["R_shell_in",  "Ang",          550.0,    [0, inf], "",       "Radius of shell inner interface"],
#              ["s_core",      "Ang",           20.0,    [0, inf], "",       "Half-width of core interface"],
#              ["s_shell_in",  "Ang",           20.0,    [0, inf], "",       "Half-width of shell inner interface"],
#              ["s_shell_out", "Ang",           20.0,    [0, inf], "",       "Half-width of shell outer interface"],
#              ["L_0",         "",              10.0,    [0, inf], "",       "Lorentzian scale factor"],
#              ["xi_intra",    "Ang",            1.0,    [0, inf], "",       "Intra-particle characteristic length"]
#             ]
# pylint: enable=bad-whitespace, line-too-long
#
#
def R_shell_out(R, s):
    """Return outer radius given total radius and surface thickenss."""
    return R - s

def phi_full(q, R, s):
    """Return the Fourier Transform of the particle density profile.
    Box profile with parabolic gradation of interfaces."""
    return (4*np.pi/((q**4)*(s**2)))*((R+s)*np.cos(q*(R+s)) + (R-s)*np.cos(q*(R-s))
           - (3/q)*np.sin(q*(R+s)) - (3/q)*np.sin(q*(R-s)) + (2/q)*np.cos(q*R)
           + (6/q)*np.sin(q*R))

def Fq(q, rho_core, rho_shell, rho_solv, R_tot, R_core, R_shell_in, s_core,
       s_shell_in, s_shell_out):
    """Return the scattering amplitude for core-shell microgels."""
    return (rho_shell - rho_solv)*phi_full(q=q, R=R_shell_out(R=R_tot, s=s_shell_out), s=s_shell_out)
    + (rho_core - rho_solv)*phi_full(q=q, R=R_core, s=s_core)
    - (rho_shell - rho_solv)*phi_full(q=q, R=R_shell_in, s=s_shell_in)

def Lor(q, L_0, xi):
    """Return a Lorenztian peak with correlation length xi."""
    return L_0/(1 + q**2*xi**2)

def Iq(q,
       rho_core=2.0,
       rho_shell=3.0,
       rho_solv=1.0,
       R_tot=700.0,
       R_core=500.0,
       R_shell_in=550.0,
       s_core=20.0,
       s_shell_in=20.0,
       s_shell_out=20.0,
       L_0=1.0,
       xi_intra=1.0):
    """Return the form factor of core-shell microgels.
    :param q:           Input q-value
    :param rho_core:    SLD of core
    :param rho_shell:   SLD of shell
    :param rho_solv:    SLD of solvent
    :param R_tot:       Radius of full particle
    :param R_core:      Radius of core
    :param R_shell_in:  Radius of shell inner interface
    :param s_core:      Half-width of core interface
    :param s_shell_in:  Half-width of shell inner interface
    :param s_shell_out: Half-width of shell outer interface
    :param L_0:         Lorentzian scale factor
    :param xi_intra:    Intra-particle characteristic length
    :return:            Calculated intensity
    """
    A = Fq(q=q, rho_core=rho_core, rho_shell=rho_shell, rho_solv=rho_solv,
    R_tot=R_tot, R_core=R_core, R_shell_in=R_shell_in, s_core=s_core,
    s_shell_in=s_shell_in, s_shell_out=s_shell_out)
    return A**2 + Lor(q=q, L_0=L_0, xi=xi_intra)

def norm(Iq, n, u):
    """Normalize Iq to number density of particles and convert to cm^-1"""
    return Iq*n*u

print(norm(Iq=Iq(q=0.1, rho_core=0.62, rho_shell=0.74, rho_solv=6.4, R_tot=1170, R_core=965, R_shell_in=965, s_core=15, s_shell_in=15, s_shell_out=20, L_0=1, xi_intra=5), n=10**-10, u=10**-8))
