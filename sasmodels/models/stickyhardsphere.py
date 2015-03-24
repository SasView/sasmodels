# Note: model title and parameter table are inserted automatically
r"""
This calculates the interparticle structure factor for a hard sphere fluid
with a narrow attractive well. A perturbative solution of the Percus-Yevick
closure is used. The strength of the attractive well is described in terms
of "stickiness" as defined below. The returned value is a dimensionless
structure factor, *S(q)*.

The perturb (perturbation parameter), |epsilon|, should be held between 0.01
and 0.1. It is best to hold the perturbation parameter fixed and let
the "stickiness" vary to adjust the interaction strength. The stickiness,
|tau|, is defined in the equation below and is a function of both the
perturbation parameter and the interaction strength. |tau| and |epsilon|
are defined in terms of the hard sphere diameter (|sigma| = 2\*\ *R*\ ), the
width of the square well, |bigdelta| (same units as *R*), and the depth of
the well, *Uo*, in units of kT. From the definition, it is clear that
smaller |tau| means stronger attraction.

.. image:: img/stickyhardsphere_228.PNG

where the interaction potential is

.. image:: img/stickyhardsphere_229.PNG

The Percus-Yevick (PY) closure was used for this calculation, and is an
adequate closure for an attractive interparticle potential. This solution
has been compared to Monte Carlo simulations for a square well fluid, with
good agreement.

The true particle volume fraction, |phi|, is not equal to *h*, which appears
in most of the reference. The two are related in equation (24) of the
reference. The reference also describes the relationship between this
perturbation solution and the original sticky hard sphere (or adhesive
sphere) model by Baxter.

NB: The calculation can go haywire for certain combinations of the input
parameters, producing unphysical solutions - in this case errors are
reported to the command window and the *S(q)* is set to -1 (so it will
disappear on a log-log plot). Use tight bounds to keep the parameters to
values that you know are physical (test them) and keep nudging them until
the optimization does not hit the constraints.

In sasview the effective radius will be calculated from the parameters
used in the form factor P(Q) that this S(Q) is combined with.

For 2D data: The 2D scattering intensity is calculated in the same way
as 1D, where the *q* vector is defined as

.. math::

    Q = \sqrt{Q_x^2 + Q_y^2}

==============  ========  =============
Parameter name  Units     Default value
==============  ========  =============
effect_radius   |Ang|     50
perturb         None      0.05
volfraction     None      0.1
stickiness      K         0.2
==============  ========  =============

.. image:: img/stickyhardsphere_230.jpg

*Figure. 1D plot using the default values (in linear scale).*

REFERENCE

S V G Menon, C Manohar, and K S Rao, *J. Chem. Phys.*, 95(12) (1991) 9186-9190
"""

# TODO: refactor so that we pull in the old sansmodels.c_extensions

from numpy import inf

name = "stickyhardsphere"
title = "Sticky hard sphere structure factor, with Percus-Yevick closure"
description = """\
    [Sticky hard sphere structure factor, with Percus-Yevick closure]
        Interparticle structure factor S(Q)for a hard sphere fluid with
        a narrow attractive well. Fits are prone to deliver non-physical
        parameters, use with care and read the references in the full manual.
        In sasview the effective radius will be calculated from the
        parameters used in P(Q).
"""
category = "structure-factor"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    #   [ "name", "units", default, [lower, upper], "type",
    #     "description" ],
    ["effect_radius", "Ang", 50.0, [0, inf], "volume",
     "effective radius of hard sphere"],
    ["volfraction", "", 0.2, [0, 0.74], "",
     "volume fraction of hard spheres"],
    ["perturb", "", 0.05, [0.01, 0.1], "",
     "perturbation parameter, epsilon"],
    ["stickiness", "", 0.20, [-inf, inf], "",
     "stickiness, tau"],
    ]

# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return 1.0;
    """

Iq = """
    double onemineps,eta;
    double sig,aa,etam1,etam1sq,qa,qb,qc,radic;
    double lam,lam2,test,mu,alpha,beta;
    double kk,k2,k3,ds,dc,aq1,aq2,aq3,aq,bq1,bq2,bq3,bq,sq;

    onemineps = 1.0-perturb;
    eta = volfraction/onemineps/onemineps/onemineps;

    sig = 2.0 * effect_radius;
    aa = sig/onemineps;
    etam1 = 1.0 - eta;
    etam1sq=etam1*etam1;
    //C
    //C  SOLVE QUADRATIC FOR LAMBDA
    //C
    qa = eta/12.0;
    qb = -1.0*(stickiness + eta/etam1);
    qc = (1.0 + eta/2.0)/etam1sq;
    radic = qb*qb - 4.0*qa*qc;
    if(radic<0) {
        //if(x>0.01 && x<0.015)
        //    Print "Lambda unphysical - both roots imaginary"
        //endif
        return(-1.0);
    }
    //C   KEEP THE SMALLER ROOT, THE LARGER ONE IS UNPHYSICAL
    lam = (-1.0*qb-sqrt(radic))/(2.0*qa);
    lam2 = (-1.0*qb+sqrt(radic))/(2.0*qa);
    if(lam2<lam) {
        lam = lam2;
    }
    test = 1.0 + 2.0*eta;
    mu = lam*eta*etam1;
    if(mu>test) {
        //if(x>0.01 && x<0.015)
        // Print "Lambda unphysical mu>test"
        //endif
        return(-1.0);
    }
    alpha = (1.0 + 2.0*eta - mu)/etam1sq;
    beta = (mu - 3.0*eta)/(2.0*etam1sq);
    //C
    //C   CALCULATE THE STRUCTURE FACTOR
    //C
    kk = q*aa;
    k2 = kk*kk;
    k3 = kk*k2;
    SINCOS(kk,ds,dc);
    //ds = sin(kk);
    //dc = cos(kk);
    aq1 = ((ds - kk*dc)*alpha)/k3;
    aq2 = (beta*(1.0-dc))/k2;
    aq3 = (lam*ds)/(12.0*kk);
    aq = 1.0 + 12.0*eta*(aq1+aq2-aq3);
    //
    bq1 = alpha*(0.5/kk - ds/k2 + (1.0 - dc)/k3);
    bq2 = beta*(1.0/kk - ds/k2);
    bq3 = (lam/12.0)*((1.0 - dc)/kk);
    bq = 12.0*eta*(bq1+bq2-bq3);
    //
    sq = 1.0/(aq*aq +bq*bq);

    return(sq);
"""

Iqxy = """
    return Iq(sqrt(qx*qx+qy*qy), IQ_PARAMETERS);
    """

# ER defaults to 0.0
# VR defaults to 1.0

oldname = 'StickyHSStructure'
oldpars = dict()
demo = dict(effect_radius=200, volfraction=0.2, perturb=0.05,
            stickiness=0.2, effect_radius_pd=0.1, effect_radius_pd_n=40)

