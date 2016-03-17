r"""
Polydispersity in the bilayer thickness can be applied from the GUI.

Definition
----------

The scattering intensity $I(q)$ is

.. math::

    I(q) = \frac{2\pi P(q)}{\delta q^2}


The form factor is

.. math::

    P(q) = \frac{2\Delta\rho^2}{q^2}(1-cos(q\delta))


where $\delta$ is the bilayer thickness.

The 2D scattering intensity is calculated in the same way as 1D, where
the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


References
----------

F Nallet, R Laversanne, and D Roux, J. Phys. II France, 3, (1993) 487-502

also in J. Phys. Chem. B, 105, (2001) 11081-11088


"""

from numpy import inf

name = "lamellar"
title = "Lyotropic lamellar phase with uniform SLD and random distribution"
description = """\
    [Dilute Lamellar Form Factor](from a lyotropic lamellar phase)
        I(q)= 2*pi*P(q)/(delta *q^(2)), where
        P(q)=2*(contrast/q)^(2)*(1-cos(q*delta))^(2))
        thickness = layer thickness
        sld = layer scattering length density
        sld_solvent = solvent scattering length density
        background = incoherent background
        scale = scale factor
"""
category = "shape:lamellae"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 1, [-inf, inf], "",
               "Layer scattering length density" ],
              ["solvent_sld", "1e-6/Ang^2", 6, [-inf, inf], "",
               "Solvent scattering length density" ],
              ["thickness", "Ang", 50, [0, inf], "volume","Bilayer thickness" ],
             ]


# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return 1.0;
    """

Iq = """
    const double sub = sld - solvent_sld;
    const double qsq = q*q;
    // Original expression
    //return 4.0e-4*M_PI*sub*sub/qsq * (1.0-cos(q*thickness)) / (thickness*qsq);
    // const double alpha = fmod(q*thickness+0.1, 2.0*M_PI)-0.1;
    // Use small angle fix 1-cos(theta) = 2 sin^2(theta/2)
    const double sinq2 = sin(0.5*q*thickness);
    return 4.0e-4*M_PI*sub*sub/qsq * 2.0*sinq2*sinq2 / (thickness*qsq);
    """

Iqxy = """
    return Iq(sqrt(qx*qx+qy*qy), IQ_PARAMETERS);
    """

# ER defaults to 0.0
# VR defaults to 1.0

demo = dict(scale=1, background=0,
            sld=6, solvent_sld=1,
            thickness=40,
            thickness_pd=0.2, thickness_pd_n=40)
oldname = 'LamellarModel'
oldpars = dict(sld='sld_bi', solvent_sld='sld_sol', thickness='bi_thick')
tests = [
        [ {'scale': 1.0, 'background' : 0.0, 'thickness' : 50.0, 'sld' : 1.0,'solvent_sld' : 6.3, 'thickness_pd' : 0.0, 
           }, [0.001], [882289.54309]]
        ]
# ADDED by: converted by PAK? (or RKH?)     ON: 16Mar2016 - RKH adding unit tests from sasview to early 2015 conversion
#  [(qx1, qy1), (qx2, qy2), ...], [I(qx1,qy1), I(qx2,qy2), ...]],