r"""
This model calculates an empirical functional form for SAS data using SpericalSLD profile

Similarly to the OnionExpShellModel, this model provides the form factor, P(q), for a multi-shell sphere,
where the interface between the each neighboring shells can be described by one of a number of functions
including error, power-law, and exponential functions. This model is to calculate the scattering intensity
by building a continuous custom SLD profile against the radius of the particle. The SLD profile is composed
of a flat core, a flat solvent, a number (up to 9 ) flat shells, and the interfacial layers between
the adjacent flat shells (or core, and solvent) (see below).
Unlike the OnionExpShellModel (using an analytical integration),
the interfacial layers here are sub-divided and numerically integrated assuming each of the sub-layers are described
by a line function. The number of the sub-layer can be given by users by setting the integer values of npts_inter
in the GUI. The form factor is normalized by the total volume of the sphere.

Definition
----------

The scattering intensity $I(q)$ in 1D is calculated as:

.. math::

    P(q) = \frac{f^2}{V_\text{particle}}
    f = f_\text{core} + \sum_{\text{inter}_i=0}^N f_text{inter}_i +
    \sum_{\text{flat}_i=0}^N f_text{flat}_i +f_\text{solvent}

where, for a spherically symmetric particle with a particle density \rho(r)


The scaling of the second power law region (coefficent C) is then automatically scaled to
match the first by following formula

.. math::
    C = \frac{A}{qc^{-m1} qc^{-m2}}

.. note::
    Be sure to enter the power law exponents as positive values!

For 2D data the scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


.. figure:: img/two_power_law_1d.jpg

    1D plot using the default values (w/500 data point).

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
category = "shere-based"


# pylint: disable=bad-whitespace, line-too-long
#            ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["n_shells",       "",         1, [0, 9], "", "number of shells"],
              ["thick_inter_0",       "Ang",         50, [-inf, inf], "", "intern layer thickness"],
              ["func_inter_0",       "",         0, [0, 4], "", "'Erf(|nu|*z)':0, 'RPower(z^|nu|)':1, 'LPower(z^|nu|)':2, 'RExp(-|nu|*z)':3, 'LExp(-|nu|*z)':4"],
              ["sld_core_0",       "1/Ang^2",         2.07E-6, [-inf, inf],"", "sld function flat"],
              ["sld_solv",     "1/Ang^2",    1E-6,[-inf, inf], "","sld function solvent"],
              ["sld_flat_1",       "1/Ang^2",         4.06E-6, [-inf, inf],"", "sld function flat"],
              ["sld_flat_2",       "1/Ang^2",         3.5E-6, [-inf, inf],"", "sld function flat"],
              ["sld_flat_3",       "1/Ang^2",         4.06E-6, [-inf, inf],"", "sld function flat"],
              ["sld_flat_4",       "1/Ang^2",         3.5E-6, [-inf, inf],"", "sld function flat"],
              ["sld_flat_5",       "1/Ang^2",         4.06E-6, [-inf, inf],"", "sld function flat"],
              ["sld_flat_6",       "1/Ang^2",         3.5E-6, [-inf, inf],"", "sld function flat"],
              ["sld_flat_7",       "1/Ang^2",         4.06E-6, [-inf, inf],"", "sld function flat"],
              ["sld_flat_8",       "1/Ang^2",         3.5E-6, [-inf, inf],"", "sld function flat"],
              ["sld_flat_9",       "1/Ang^2",         4.06E-6, [-inf, inf],"", "sld function flat"],
              ["sld_flat_10",      "1/Ang^2",         3.5E-6, [-inf, inf],"", "sld function flat"],
              ["thick_inter_1",       "Ang",         50.0, [-inf, inf], "", "intern layer thickness"],
              ["thick_inter_2",       "Ang",         50.0, [-inf, inf], "", "intern layer thickness"],
              ["thick_inter_3",       "Ang",         50.0, [-inf, inf], "", "intern layer thickness"],
              ["thick_inter_4",       "Ang",         50.0, [-inf, inf], "", "intern layer thickness"],
              ["thick_inter_5",       "Ang",         50.0, [-inf, inf], "", "intern layer thickness"],
              ["thick_inter_6",       "Ang",         50.0, [-inf, inf], "", "intern layer thickness"],
              ["thick_inter_7",       "Ang",         50.0, [-inf, inf], "", "intern layer thickness"],
              ["thick_inter_8",       "Ang",         50.0, [-inf, inf], "", "intern layer thickness"],
              ["thick_inter_9",       "Ang",         50.0, [-inf, inf], "", "intern layer thickness"],
              ["thick_inter_10",       "Ang",         50.0, [-inf, inf], "", "intern layer thickness"],
              ["thick_flat_1",       "Ang",         100.0, [-inf, inf], "", "flat layer_thickness"],
              ["thick_flat_2",       "Ang",         100.0, [-inf, inf], "", "flat layer_thickness"],
              ["thick_flat_3",       "Ang",         100.0, [-inf, inf], "", "flat layer_thickness"],
              ["thick_flat_4",       "Ang",         100.0, [-inf, inf], "", "flat layer_thickness"],
              ["thick_flat_5",       "Ang",         100.0, [-inf, inf], "", "flat layer_thickness"],
              ["thick_flat_6",       "Ang",         100.0, [-inf, inf], "", "flat layer_thickness"],
              ["thick_flat_7",       "Ang",         100.0, [-inf, inf], "", "flat layer_thickness"],
              ["thick_flat_8",       "Ang",         100.0, [-inf, inf], "", "flat layer_thickness"],
              ["thick_flat_9",       "Ang",         100.0, [-inf, inf], "", "flat layer_thickness"],
              ["thick_flat_10",       "Ang",         100.0, [-inf, inf], "", "flat layer_thickness"],
              ["func_inter_1",       "",         0, [0, 4], "", "'Erf(|nu|*z)':0, 'RPower(z^|nu|)':1, 'LPower(z^|nu|)':2, 'RExp(-|nu|*z)':3, 'LExp(-|nu|*z)':4"],
              ["func_inter_2",       "",         0, [0, 4], "", "'Erf(|nu|*z)':0, 'RPower(z^|nu|)':1, 'LPower(z^|nu|)':2, 'RExp(-|nu|*z)':3, 'LExp(-|nu|*z)':4"],
              ["func_inter_3",       "",         0, [0, 4], "", "'Erf(|nu|*z)':0, 'RPower(z^|nu|)':1, 'LPower(z^|nu|)':2, 'RExp(-|nu|*z)':3, 'LExp(-|nu|*z)':4"],
              ["func_inter_4",       "",         0, [0, 4], "", "'Erf(|nu|*z)':0, 'RPower(z^|nu|)':1, 'LPower(z^|nu|)':2, 'RExp(-|nu|*z)':3, 'LExp(-|nu|*z)':4"],
              ["func_inter_5",       "",         0, [0, 4], "", "'Erf(|nu|*z)':0, 'RPower(z^|nu|)':1, 'LPower(z^|nu|)':2, 'RExp(-|nu|*z)':3, 'LExp(-|nu|*z)':4"],
              ["func_inter_6",       "",         0, [0, 4], "", "'Erf(|nu|*z)':0, 'RPower(z^|nu|)':1, 'LPower(z^|nu|)':2, 'RExp(-|nu|*z)':3, 'LExp(-|nu|*z)':4"],
              ["func_inter_7",       "",         0, [0, 4], "", "'Erf(|nu|*z)':0, 'RPower(z^|nu|)':1, 'LPower(z^|nu|)':2, 'RExp(-|nu|*z)':3, 'LExp(-|nu|*z)':4"],
              ["func_inter_8",       "",         0, [0, 4], "", "'Erf(|nu|*z)':0, 'RPower(z^|nu|)':1, 'LPower(z^|nu|)':2, 'RExp(-|nu|*z)':3, 'LExp(-|nu|*z)':4"],
              ["func_inter_9",       "",         0, [0, 4], "", "'Erf(|nu|*z)':0, 'RPower(z^|nu|)':1, 'LPower(z^|nu|)':2, 'RExp(-|nu|*z)':3, 'LExp(-|nu|*z)':4"],
              ["func_inter_10",       "",        0, [0, 4], "", "'Erf(|nu|*z)':0, 'RPower(z^|nu|)':1, 'LPower(z^|nu|)':2, 'RExp(-|nu|*z)':3, 'LExp(-|nu|*z)':4"],
              ["nu_inter_1",       "",         2.5, [-inf, inf], "", "steepness parameter"],
              ["nu_inter_2",       "",         2.5, [-inf, inf], "", "steepness parameter"],
              ["nu_inter_3",       "",         2.5, [-inf, inf], "", "steepness parameter"],
              ["nu_inter_4",       "",         2.5, [-inf, inf], "", "steepness parameter"],
              ["nu_inter_5",       "",         2.5, [-inf, inf], "", "steepness parameter"],
              ["nu_inter_6",       "",         2.5, [-inf, inf], "", "steepness parameter"],
              ["nu_inter_7",       "",         2.5, [-inf, inf], "", "steepness parameter"],
              ["nu_inter_8",       "",         2.5, [-inf, inf], "", "steepness parameter"],
              ["nu_inter_9",       "",         2.5, [-inf, inf], "", "steepness parameter"],
              ["nu_inter_10",       "",         2.5, [-inf, inf], "", "steepness parameter"],
              ["npts_inter",  "",         35, [0, inf], "", "number of points in each sublayer"],
              ["nu_inter_0",       "",         2.5, [-inf, inf], "", "steepness parameter"],
              ["rad_core_0",       "Ang",         50.0, [-inf, inf], "", "intern layer thickness"],
              ]
# pylint: enable=bad-whitespace, line-too-long
source = ["lib/librefl.c", "spherical_sld.c"]


#TODO: Not clear if dispersion function is needed
#def _set_dispersion(self):
#        """
#        Setting dispersion for all shells
#        """
#        ##set dispersion from model
#        for name , value in self.model.dispersion.iteritems():
#
#            nshell = -1
#            if name.split('_')[0] == 'thick':
#                while nshell < 1:
#                    nshell += 1
#                    if name.split('_')[1] == 'inter%s' % str(nshell):
#                        self.dispersion[name] = value
#                    else:
#                        continue
#            else:
#                self.dispersion[name] = value


def ER(radius):
        """
        Calculate the effective radius for P(q)*S(q)
        Tale different radius values from corresponding shells

        :return: the value of the effective radius
        """

        return radius

demo = dict(scale=1, background=0.0,
        n_shells=10,
        sld_solv=1E-6,
        npts_inter=35,
        func_inter_0=0,
        nu_inter_0=2.5,
        rad_core_0=50.0,
        sld_core_0=2.07E-6,
        thick_inter_0=50,
        func_inter_1=0,
        nu_inter_1=2.5,
        thick_inter_1=50,
        sld_flat_1=4E-6,
        thick_flat_1=100,
        func_inter_2=0,
        nu_inter_2=2.5,
        thick_inter_2=50,
        sld_flat_2=3.5E-6,
        thick_flat_2=100,
        func_inter_3=0,
        nu_inter_3=2.5,
        thick_inter_3=50,
        sld_flat_3=4E-6,
        thick_flat_3=100,
        func_inter_4=0,
        nu_inter_4=2.5,
        thick_inter_4=50,
        sld_flat_4=3.5E-6,
        thick_flat_4=100,
        func_inter_5=0,
        nu_inter_5=2.5,
        thick_inter_5=50,
        sld_flat_5=4E-6,
        thick_flat_5=100,
        func_inter_6=0,
        nu_inter_6=2.5,
        thick_inter_6=50,
        sld_flat_6=3.5E-6,
        thick_flat_6=100,
        func_inter_7=0,
        nu_inter_7=2.5,
        thick_inter_7=50,
        sld_flat_7=4E-6,
        thick_flat_7=100,
        func_inter_8=0,
        nu_inter_8=2.5,
        thick_inter_8=50,
        sld_flat_8=3.5E-6,
        thick_flat_8=100,
        func_inter_9=0,
        nu_inter_9=2.5,
        thick_inter_9=50,
        sld_flat_9=4E-6,
        thick_flat_9=100,
        func_inter_10=0,
        nu_inter_10=2.5,
        thick_inter_10=50,
        sld_flat_10=3.5E-6,
        thick_flat_10=100
        )

oldname = "SphereSLDModel"
oldpars = dict(
        scale="scale",
        background="background",
        n_shells="n_shells",
        npts_inter='npts_inter',
        sld_solv='sld_solv',
        func_inter_0='func_inter0',
        nu_inter_0='nu_inter0',
        rad_core_0='rad_core0',
        sld_core_0='sld_core0',
        thick_inter_0='thick_inter0',
        func_inter_1='func_inter1',
        nu_inter_1='nu_inter1',
        thick_inter_1='thick_inter1',
        sld_flat_1='sld_flat1',
        thick_flat_1='thick_flat1',
        func_inter_2='func_inter2',
        nu_inter_2='nu_inter2',
        thick_inter_2='thick_inter2',
        sld_flat_2='sld_flat2',
        thick_flat_2='thick_flat2',
        func_inter_3='func_inter3',
        nu_inter_3='nu_inter3',
        thick_inter_3='thick_inter3',
        sld_flat_3='sld_flat3',
        thick_flat_3='thick_flat3',
        func_inter_4='func_inter4',
        nu_inter_4='nu_inter4',
        thick_inter_4='thick_inter4',
        sld_flat_4='sld_flat4',
        thick_flat_4='thick_flat4',
        func_inter_5='func_inter5',
        nu_inter_5='nu_inter5',
        thick_inter_5='thick_inter5',
        sld_flat_5='sld_flat5',
        thick_flat_5='thick_flat5',
        func_inter_6='func_inter6',
        nu_inter_6='nu_inter6',
        thick_inter_6='thick_inter6',
        sld_flat_6='sld_flat6',
        thick_flat_6='thick_flat6',
        func_inter_7='func_inter7',
        nu_inter_7='nu_inter7',
        thick_inter_7='thick_inter7',
        sld_flat_7='sld_flat7',
        thick_flat_7='thick_flat7',
        func_inter_8='func_inter8',
        nu_inter_8='nu_inter8',
        thick_inter_8='thick_inter8',
        sld_flat_8='sld_flat8',
        thick_flat_8='thick_flat8',
        func_inter_9='func_inter9',
        nu_inter_9='nu_inter9',
        thick_inter_9='thick_inter9',
        sld_flat_9='sld_flat9',
        thick_flat_9='thick_flat9',
        func_inter_10='func_inter10',
        nu_inter_10='nu_inter10',
        thick_inter_10='thick_inter10',
        sld_flat_10='sld_flat10',
        thick_flat_10='thick_flat10')

tests = [
    # Accuracy tests based on content in test/utest_extra_models.py
    [{'npts_iter':35,
        'sld_solv':1E-6,
        'func_inter_1':0.0,
        'nu_inter':2.5,
        'thick_inter_1':50,
        'sld_flat_1':4E-6,
        'thick_flat_1':100,
        'func_inter_0':0.0,
        'nu_inter_0':2.5,
        'rad_core_0':50.0,
        'sld_core_0':2.07E-6,
        'thick_inter_0':50,
        'background': 0.0,
    }, 0.001, 1000],
]
