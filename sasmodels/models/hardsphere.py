# Note: model title and parameter table are inserted automatically
r"""Calculate the interparticle structure factor for monodisperse
spherical particles interacting through hard sphere (excluded volume)
interactions.

The calculation uses the Percus-Yevick closure where the interparticle
potential is

.. math:

    U(r) = \begin{cases}
    \infty & r < 2R \\
    0 & r \geq 2R
    \end{cases}

where *r* is the distance from the center of the sphere of a radius *R*.

For a 2D plot, the wave transfer is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


.. image:: img/HardSphere_1d.jpg

*Figure. 1D plot using the default values (in linear scale).*

REFERENCE

J K Percus, J Yevick, *J. Phys. Rev.*, 110, (1958) 1
"""

from numpy import inf

name = "hardsphere"
title = "Hard sphere structure factor, with Percus-Yevick closure"
description = """\
    [Hard sphere structure factor, with Percus-Yevick closure]
        Interparticle S(Q) for random, non-interacting spheres.
    May be a reasonable approximation for other shapes of
    particles that freely rotate, and for moderately polydisperse
        systems. Though strictly the maths needs to be modified -
    which sasview does not do yet.
    effect_radius is the hard sphere radius
    volfraction is the volume fraction occupied by the spheres.
"""
category = "structure-factor"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["effect_radius", "Ang", 50.0, [0, inf], "volume",
               "effective radius of hard sphere"],
              ["volfraction", "", 0.2, [0, 0.74], "",
               "volume fraction of hard spheres"],
             ]

# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return 1.0;
    """

Iq = """
    double denom,dnum,alpha,beta,gamm,a,asq,ath,afor,rca,rsa;
    double calp,cbeta,cgam,prefac,c,vstruc;
    double struc;

    //  compute constants
    denom = pow((1.0-volfraction),4);
    dnum = pow((1.0 + 2.0*volfraction),2);
    alpha = dnum/denom;
    beta = -6.0*volfraction*pow((1.0 + volfraction/2.0),2)/denom;
    gamm = 0.50*volfraction*dnum/denom;
    //
    //  calculate the structure factor
    //
    a = 2.0*q*effect_radius;
    asq = a*a;
    ath = asq*a;
    afor = ath*a;
    SINCOS(a,rsa,rca);
    //rca = cos(a);
    //rsa = sin(a);
    calp = alpha*(rsa/asq - rca/a);
    cbeta = beta*(2.0*rsa/asq - (asq - 2.0)*rca/ath - 2.0/ath);
    cgam = gamm*(-rca/a + (4.0/a)*((3.0*asq - 6.0)*rca/afor + (asq - 6.0)*rsa/ath + 6.0/afor));
    prefac = -24.0*volfraction/a;
    c = prefac*(calp + cbeta + cgam);
    vstruc = 1.0/(1.0-c);
    struc = vstruc;

    return(struc);
   """

Iqxy = """
    // never called since no orientation or magnetic parameters.
    return Iq(sqrt(qx*qx+qy*qy), IQ_PARAMETERS);
    """

# ER defaults to 0.0
# VR defaults to 1.0

demo = dict(effect_radius=200, volfraction=0.2, effect_radius_pd=0.1, effect_radius_pd_n=40)
oldname = 'HardsphereStructure'
oldpars = dict()

