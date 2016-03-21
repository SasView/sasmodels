r"""
This model calculates an empirical functional form for SAS data using SpericalSLD profile

Similarly to the OnionExpShellModel, this model provides the form factor, P(q), for a multi-shell sphere,
where the interface between the each neighboring shells can be described by one of a number of functions
including error, power-law, and exponential functions. This model is to calculate the scattering intensity
by building a continuous custom SLD profile against the radius of the particle. The SLD profile is composed
of a flat core, a flat solvent, a number (up to 9 ) flat shells, and the interfacial layers between
the adjacent flat shells (or core, and solvent) (see below).

.. figure:: img/spherical_sld_profile.gif

    Exemplary SLD profile

Unlike the OnionExpShellModel (using an analytical integration),
the interfacial layers here are sub-divided and numerically integrated assuming each of the sub-layers are described
by a line function. The number of the sub-layer can be given by users by setting the integer values of npts_inter.
The form factor is normalized by the total volume of the sphere.

Definition
----------

The form factor $P(q)$ in 1D is calculated by:

.. math::

    P(q) = \frac{f^2}{V_\text{particle}} \text{ where }
    f = f_\text{core} + \sum_{\text{inter}_i=0}^N f_{\text{inter}_i} +
    \sum_{\text{flat}_i=0}^N f_{\text{flat}_i} +f_\text{solvent}

For a spherically symmetric particle with a particle density $\rho_x(r)$ the sld function can be defined as:

.. math::

    f_x = 4 \pi \int_{0}^{\infty} \rho_x(r)  \frac{\sin(qr)} {qr^2} r^2 dr


so that individual terms can be calcualted as follows:

.. math::
    f_\text{core} = 4 \pi \int_{0}^{r_\text{core}} \rho_\text{core} \frac{\sin(qr)} {qr} r^2 dr =
    3 \rho_\text{core} V(r_\text{core})
    \Big[ \frac{\sin(qr_\text{core}) - qr_\text{core} \cos(qr_\text{core})} {qr_\text{core}^3} \Big]

    f_{\text{inter}_i} = 4 \pi \int_{\Delta t_{ \text{inter}_i } } \rho_{ \text{inter}_i } \frac{\sin(qr)} {qr} r^2 dr

    f_{\text{shell}_i} = 4 \pi \int_{\Delta t_{ \text{inter}_i } } \rho_{ \text{flat}_i } \frac{\sin(qr)} {qr} r^2 dr =
    3 \rho_{ \text{flat}_i } V ( r_{ \text{inter}_i } + \Delta t_{ \text{inter}_i } )
    \Big[ \frac{\sin(qr_{\text{inter}_i} + \Delta t_{ \text{inter}_i } ) - q (r_{\text{inter}_i} +
    \Delta t_{ \text{inter}_i }) \cos(q( r_{\text{inter}_i} + \Delta t_{ \text{inter}_i } ) ) }
    {q ( r_{\text{inter}_i} + \Delta t_{ \text{inter}_i } )^3 }  \Big]
    -3 \rho_{ \text{flat}_i } V(r_{ \text{inter}_i })
    \Big[ \frac{\sin(qr_{\text{inter}_i}) - qr_{\text{flat}_i} \cos(qr_{\text{inter}_i}) } {qr_{\text{inter}_i}^3} \Big]

    f_\text{solvent} = 4 \pi \int_{r_N}^{\infty} \rho_\text{solvent} \frac{\sin(qr)} {qr} r^2 dr =
    3 \rho_\text{solvent} V(r_N)
    \Big[ \frac{\sin(qr_N) - qr_N \cos(qr_N)} {qr_N^3} \Big]


Here we assumed that the SLDs of the core and solvent are constant against $r$.
The SLD at the interface between shells, $\rho_{\text {inter}_i}$
is calculated with a function chosen by an user, where the functions are

Exp:

.. math::
    \rho_{{inter}_i} (r) = \begin{cases}
    B \exp\Big( \frac {\pm A(r - r_{\text{flat}_i})} {\Delta t_{ \text{inter}_i }} \Big) +C  & \text{for} A \neq 0 \\
    B \Big( \frac {(r - r_{\text{flat}_i})} {\Delta t_{ \text{inter}_i }} \Big) +C  & \text{for} A = 0 \\
    \end{cases}

Power-Law

.. math::
    \rho_{{inter}_i} (r) = \begin{cases}
    \pm B \Big( \frac {(r - r_{\text{flat}_i} )} {\Delta t_{ \text{inter}_i }} \Big) ^A  +C  & \text{for} A \neq 0 \\
    \rho_{\text{flat}_{i+1}}  & \text{for} A = 0 \\
    \end{cases}

Erf:

.. math::
    \rho_{{inter}_i} (r) = \begin{cases}
    B \text{erf} \Big( \frac { A(r - r_{\text{flat}_i})} {\sqrt{2} \Delta t_{ \text{inter}_i }} \Big) +C  & \text{for} A \neq 0 \\
    B \Big( \frac {(r - r_{\text{flat}_i} )} {\Delta t_{ \text{inter}_i }} \Big)  +C  & \text{for} A = 0 \\
    \end{cases}

The functions are normalized so that they vary between 0 and 1, and they are constrained such that the SLD
is continuous at the boundaries of the interface as well as each sub-layers. Thus B and C are determined.

Once $\rho_{\text{inter}_i}$ is found at the boundary of the sub-layer of the interface, we can find its contribution
to the form factor $P(q)$

.. math::
    f_{\text{inter}_i} = 4 \pi \int_{\Delta t_{ \text{inter}_i } } \rho_{ \text{inter}_i } \frac{\sin(qr)} {qr} r^2 dr =
    4 \pi \sum_{j=0}^{npts_{\text{inter}_i} -1 }
    \int_{r_j}^{r_{j+1}} \rho_{ \text{inter}_i } (r_j) \frac{\sin(qr)} {qr} r^2 dr \approx

    4 \pi \sum_{j=0}^{npts_{\text{inter}_i} -1 } \Big[
    3 ( \rho_{ \text{inter}_i } ( r_{j+1} ) - \rho_{ \text{inter}_i } ( r_{j} ) V ( r_{ \text{sublayer}_j } )
    \Big[ \frac {r_j^2 \beta_\text{out}^2 \sin(\beta_\text{out}) - (\beta_\text{out}^2-2) \cos(\beta_\text{out}) }
    {\beta_\text{out}^4 } \Big]

    - 3 ( \rho_{ \text{inter}_i } ( r_{j+1} ) - \rho_{ \text{inter}_i } ( r_{j} ) V ( r_{ \text{sublayer}_j-1 } )
    \Big[ \frac {r_{j-1}^2 \sin(\beta_\text{in}) - (\beta_\text{in}^2-2) \cos(\beta_\text{in}) }
    {\beta_\text{in}^4 } \Big]

    + 3 \rho_{ \text{inter}_i } ( r_{j+1} )  V ( r_j )
    \Big[ \frac {\sin(\beta_\text{out}) - \cos(\beta_\text{out}) }
    {\beta_\text{out}^4 } \Big]

    - 3 \rho_{ \text{inter}_i } ( r_{j} )  V ( r_j )
    \Big[ \frac {\sin(\beta_\text{in}) - \cos(\beta_\text{in}) }
    {\beta_\text{in}^4 } \Big]
    \Big]

where

.. math::
    V(a) = \frac {4\pi}{3}a^3

    a_\text{in} ~ \frac{r_j}{r_{j+1} -r_j} \text{, } a_\text{out} ~ \frac{r_{j+1}}{r_{j+1} -r_j}

    \beta_\text{in} = qr_j \text{, } \beta_\text{out} = qr_{j+1}


We assume the $\rho_{\text{inter}_i} (r)$ can be approximately linear within a sub-layer $j$

Finally form factor can be calculated by

.. math::

    P(q) = \frac{[f]^2} {V_\text{particle}} \text{where} V_\text{particle} = V(r_{\text{shell}_N})

For 2D data the scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


.. figure:: img/spherical_sld_1d.jpg

    1D plot using the default values (w/400 data point).

.. figure:: img/spherical_sld_default_profile.jpg

    SLD profile from the default values.

.. note::
    The outer most radius is used as the effective radius for S(Q) when $P(Q) * S(Q)$ is applied.

References
----------
L A Feigin and D I Svergun, Structure Analysis by Small-Angle X-Ray and Neutron Scattering, Plenum Press, New York, (1987)

"""

from numpy import inf

name = "spherical_sld"
title = "Sperical SLD intensity calculation"
description = """
            I(q) =
               background = Incoherent background [1/cm]
        """
category = "sphere-based"

# pylint: disable=bad-whitespace, line-too-long
#            ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["n_shells",         "",               1,      [0, 9],         "", "number of shells"],
              ["thick_inter[n]",   "Ang",            50,     [-inf, inf],    "", "intern layer thickness"],
              ["func_inter[n]",    "",               0,      [0, 4],         "", "Erf:0, RPower:1, LPower:2, RExp:3, LExp:4"],
              ["sld_core",         "1e-6/Ang^2",     2.07,   [-inf, inf],    "", "sld function flat"],
              ["sld_solvent",      "1e-6/Ang^2",     1.0,    [-inf, inf],    "", "sld function solvent"],
              ["sld_flat[n]",      "1e-6/Ang^2",     4.06,   [-inf, inf],    "", "sld function flat"],
              ["thick_inter[n]",   "Ang",            50.0,   [0, inf],    "", "intern layer thickness"],
              ["thick_flat[n]",    "Ang",            100.0,  [0, inf],    "", "flat layer_thickness"],
              ["inter_nu[n]",      "",               2.5,    [-inf, inf],    "", "steepness parameter"],
              ["npts_inter",       "",               35,     [0, 35],        "", "number of points in each sublayer Must be odd number"],
              ["core_rad",         "Ang",            50.0,   [0, inf],    "", "intern layer thickness"],
              ]
# pylint: enable=bad-whitespace, line-too-long
#source = ["lib/librefl.c",  "lib/sph_j1c.c", "spherical_sld.c"]

Iq = """
    return q;
    """

Iqxy = """
    // never called since no orientation or magnetic parameters.
    //return -1.0;
    """

demo = dict(
        n_shells=4,
        scale=1.0,
        solvent_sld=1.0,
        background=0.0,
        npts_inter=35.0,
        func_inter_0=0,
        nu_inter_0=2.5,
        rad_core_0=50.0,
        core0_sld=2.07,
        thick_inter_0=50.0,
        func_inter_1=0,
        nu_inter_1=2.5,
        thick_inter_1=50.0,
        flat1_sld=4.0,
        thick_flat_1=100.0,
        func_inter_2=0,
        nu_inter_2=2.5,
        thick_inter_2=50.0,
        flat2_sld=3.5,
        thick_flat_2=100.0,
        func_inter_3=0,
        nu_inter_3=2.5,
        thick_inter_3=50.0,
        flat3_sld=4.0,
        thick_flat_3=100.0,
        func_inter_4=0,
        nu_inter_4=2.5,
        thick_inter_4=50.0,
        flat4_sld=3.5,
        thick_flat_4=100.0,
        func_inter_5=0,
        nu_inter_5=2.5,
        thick_inter_5=50.0,
        flat5_sld=4.0,
        thick_flat_5=100.0,
        func_inter_6=0,
        nu_inter_6=2.5,
        thick_inter_6=50.0,
        flat6_sld=3.5,
        thick_flat_6=100.0,
        func_inter_7=0,
        nu_inter_7=2.5,
        thick_inter_7=50.0,
        flat7_sld=4.0,
        thick_flat_7=100.0,
        func_inter_8=0,
        nu_inter_8=2.5,
        thick_inter_8=50.0,
        flat8_sld=3.5,
        thick_flat_8=100.0,
        func_inter_9=0,
        nu_inter_9=2.5,
        thick_inter_9=50.0,
        flat9_sld=4.0,
        thick_flat_9=100.0,
        func_inter_10=0,
        nu_inter_10=2.5,
        thick_inter_10=50.0,
        flat10_sld=3.5,
        thick_flat_10=100.0
        )

oldname = "SphereSLDModel"
oldpars = dict(
        n_shells="n_shells",
        scale="scale",
        npts_inter='npts_inter',
        solvent_sld='sld_solv',
        func_inter_0='func_inter0',
        nu_inter_0='nu_inter0',
        background='background',
        rad_core_0='rad_core0',
        core0_sld='sld_core0',
        thick_inter_0='thick_inter0',
        func_inter_1='func_inter1',
        nu_inter_1='nu_inter1',
        thick_inter_1='thick_inter1',
        flat1_sld='sld_flat1',
        thick_flat_1='thick_flat1',
        func_inter_2='func_inter2',
        nu_inter_2='nu_inter2',
        thick_inter_2='thick_inter2',
        flat2_sld='sld_flat2',
        thick_flat_2='thick_flat2',
        func_inter_3='func_inter3',
        nu_inter_3='nu_inter3',
        thick_inter_3='thick_inter3',
        flat3_sld='sld_flat3',
        thick_flat_3='thick_flat3',
        func_inter_4='func_inter4',
        nu_inter_4='nu_inter4',
        thick_inter_4='thick_inter4',
        flat4_sld='sld_flat4',
        thick_flat_4='thick_flat4',
        func_inter_5='func_inter5',
        nu_inter_5='nu_inter5',
        thick_inter_5='thick_inter5',
        flat5_sld='sld_flat5',
        thick_flat_5='thick_flat5',
        func_inter_6='func_inter6',
        nu_inter_6='nu_inter6',
        thick_inter_6='thick_inter6',
        flat6_sld='sld_flat6',
        thick_flat_6='thick_flat6',
        func_inter_7='func_inter7',
        nu_inter_7='nu_inter7',
        thick_inter_7='thick_inter7',
        flat7_sld='sld_flat7',
        thick_flat_7='thick_flat7',
        func_inter_8='func_inter8',
        nu_inter_8='nu_inter8',
        thick_inter_8='thick_inter8',
        flat8_sld='sld_flat8',
        thick_flat_8='thick_flat8',
        func_inter_9='func_inter9',
        nu_inter_9='nu_inter9',
        thick_inter_9='thick_inter9',
        flat9_sld='sld_flat9',
        thick_flat_9='thick_flat9',
        func_inter_10='func_inter10',
        nu_inter_10='nu_inter10',
        thick_inter_10='thick_inter10',
        flat10_sld='sld_flat10',
        thick_flat_10='thick_flat10')

#TODO: Not working yet
tests = [
    # Accuracy tests based on content in test/utest_extra_models.py
    [{'npts_iter':35,
        'sld_solv':1,
        'func_inter_1':0.0,
        'nu_inter':2.5,
        'thick_inter_1':50,
        'sld_flat_1':4,
        'thick_flat_1':100,
        'func_inter_0':0.0,
        'nu_inter_0':2.5,
        'rad_core_0':50.0,
        'sld_core_0':2.07,
        'thick_inter_0':50,
        'background': 0.0,
    }, 0.001, 1000],
]
