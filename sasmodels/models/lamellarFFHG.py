# Note: model title and parameter table are inserted automatically
r"""
This model provides the scattering intensity, *I(q)*, for a lyotropic lamellar
phase where a random distribution in solution are assumed. The SLD of the head
region is taken to be different from the SLD of the tail region.

*2.1.31.1. Definition*

The scattering intensity *I(q)* is

.. math::

    I(Q) = 2\pi{P(Q) \over (2(|delta|\ H +|delta|\ T) Q^2)

The form factor is

.. image:: img/lamellarFFHG_.jpg

where |delta|\ T = tail length (or *tail_length*), |delta|\ H = head thickness
(or *h_thickness*), |drho|\ H = SLD(headgroup) - SLD(solvent),
and |drho|\ T = SLD(tail) - SLD(solvent).

The 2D scattering intensity is calculated in the same way as 1D, where
the *q* vector is defined as

.. math::

    Q = \sqrt{Q_x^2 + Q_y^2}

The returned value is in units of |cm^-1|, on absolute scale. In the
parameters, *sld_tail* = SLD of the tail group, and *sld_head* = SLD
of the head group.

.. image:: img/lamellarFFHG_138.jpg

*Figure. 1D plot using the default values (w/1000 data point).*

Our model uses the form factor calculations implemented in a C library
provided by the NIST Center for Neutron Research (Kline, 2006).

REFERENCE

F Nallet, R Laversanne, and D Roux, J. Phys. II France, 3, (1993) 487-502

also in J. Phys. Chem. B, 105, (2001) 11081-11088

*2014/04/17 - Description reviewed by S King and P Butler.*
"""

from numpy import inf

name = "lamellar_FFHG"
title = "Random lamellar phase with Head Groups "
description = """\
    [Random lamellar phase with Head Groups]
        I(q)= 2*pi*P(q)/(2(H+T)*q^(2)), where
        P(q)= see manual
        layer thickness =(H+T+T+H) = 2(Head+Tail)
        sld = Tail scattering length density
        sld_head = Head scattering length density
        sld_solvent = solvent scattering length density
        background = incoherent background
        scale = scale factor
"""
category = "shape:lamellae"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["tail_length", "Ang",  15, [0, inf], "volume",
               "Tail thickness"],
              ["head_length", "Ang",  10, [0, inf], "volume",
               "head thickness"],
              ["sld", "1e-6/Ang^2", 0.4, [-inf,inf], "",
               "Tail scattering length density"],
              ["head_sld", "1e-6/Ang^2", 3.0, [-inf,inf], "",
               "Head scattering length density"],
              ["solvent_sld", "1e-6/Ang^2", 6, [-inf,inf], "",
               "Solvent scattering length density"],
             ]

# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return 1.0;
    """

Iq = """
    const double qsq = q*q;
    const double drh = head_sld - solvent_sld;
    const double drt = sld - solvent_sld;    //correction 13FEB06 by L.Porcar
    const double qT = q*tail_length;
    double Pq, inten;
    Pq = drh*(sin(q*(head_length+tail_length))-sin(qT)) + drt*sin(qT);
    Pq *= Pq;
    Pq *= 4.0/(qsq);

    inten = 2.0e-4*M_PI*Pq/qsq;

    // normalize by the bilayer thickness
    inten /= 2.0*(head_length+tail_length);

    return inten;
    """

Iqxy = """
    return Iq(sqrt(qx*qx+qy*qy), IQ_PARAMETERS);
    """

# ER defaults to 0.0
# VR defaults to 1.0

demo = dict(scale=1, background=0,
            tail_length=15,head_length=10,
            sld=0.4, head_sld=3.0, solvent_sld=6.0,
            tail_length_pd= 0.2, tail_length_pd_n=40,
            head_length_pd= 0.01, head_length_pd_n=40)

oldname = 'LamellarFFHGModel'
oldpars = dict(head_length='h_thickness', sld='sld_tail',
               head_sld='sld_head', solvent_sld='sld_solvent')

