# Note: model title and parameter table are inserted automatically
r"""
This calculates the interparticle structure factor for a square well fluid
spherical particles. The mean spherical approximation (MSA) closure was
used for this calculation, and is not the most appropriate closure for
an attractive interparticle potential. This solution has been compared
to Monte Carlo simulations for a square well fluid, showing this calculation
to be limited in applicability to well depths $\epsilon < 1.5$ kT and
volume fractions $\phi < 0.08$.

Positive well depths correspond to an attractive potential well. Negative
well depths correspond to a potential "shoulder", which may or may not be
physically reasonable. The stickyhardsphere model may be a better choice in
some circumstances. Computed values may behave badly at extremely small $qR$.

The well width $(\lambda)$ is defined as multiples of the particle diameter
$(2 R)$.

The interaction potential is:

  .. image:: img/squarewell.png

.. math::

    U(r) = \begin{cases}
    \infty & r < 2R \\
    -\epsilon & 2R \leq r < 2R\lambda \\
    0 & r \geq 2R\lambda
    \end{cases}

where $r$ is the distance from the center of the sphere of a radius $R$.

In sasview the effective radius may be calculated from the parameters
used in the form factor $P(q)$ that this $S(q)$ is combined with.

For 2D data: The 2D scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}

References
----------

R V Sharma, K C Sharma, *Physica*, 89A (1977) 213.
"""

import numpy as np
from numpy import inf

name = "squarewell"
title = "Square well structure factor, with MSA closure"
description = """\
    [Square well structure factor, with MSA closure]
        Interparticle structure factor S(Q)for a hard sphere fluid with
        a narrow attractive well. Fits are prone to deliver non-physical
        parameters, use with care and read the references in the full manual.
        In sasview the effective radius will be calculated from the
        parameters used in P(Q).
"""
category = "structure-factor"
structure_factor = True
single = False

#single = False
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    #   [ "name", "units", default, [lower, upper], "type",
    #     "description" ],
    ["radius_effective", "Ang", 50.0, [0, inf], "volume",
     "effective radius of hard sphere"],
    ["volfraction", "", 0.04, [0, 0.08], "",
     "volume fraction of spheres"],
    ["welldepth", "kT", 1.5, [0.0, 1.5], "",
     "depth of well, epsilon"],
    ["wellwidth", "diameters", 1.2, [1.0, inf], "",
     "width of well in diameters (=2R) units, must be > 1"],
    ]

# No volume normalization despite having a volume parameter
# This should perhaps be volume normalized?
form_volume = """
    return 1.0;
    """

Iq = """
// single precision is very poor at extreme small Q, would need a Taylor series
	double req,phis,edibkb,lambda,struc;
	double sigma,eta,eta2,eta3,eta4,etam1,etam14,alpha,beta,gamm;
	double x,sk,sk2,sk3,sk4,t1,t2,t3,t4,ck;
	double S,C,SL,CL;
	x= q;

	req = radius_effective;
	phis = volfraction;
	edibkb = welldepth;
	lambda = wellwidth;

	sigma = req*2.;
	eta = phis;
	eta2 = eta*eta;
	eta3 = eta*eta2;
	eta4 = eta*eta3;
	etam1 = 1. - eta;
	etam14 = etam1*etam1*etam1*etam1;
	// temp borrow sk for an intermediate calc
	sk = 1.0 +2.0*eta;
	sk *= sk;
	alpha = (  sk + eta3*( eta-4.0 )  )/etam14;
	beta = -(eta/3.0) * ( 18. + 20.*eta - 12.*eta2 + eta4 )/etam14;
	gamm = 0.5*eta*( sk + eta3*(eta-4.) )/etam14;

	//  calculate the structure factor

	sk = x*sigma;
	sk2 = sk*sk;
	sk3 = sk*sk2;
	sk4 = sk3*sk;
	SINCOS(sk,S,C);
	SINCOS(lambda*sk,SL,CL);
	t1 = alpha * sk3 * ( S - sk * C );
	t2 = beta * sk2 * 2.0*( sk*S - (0.5*sk2 - 1.)*C - 1.0 );
	t3 = gamm*( ( 4.0*sk3 - 24.*sk ) * S - ( sk4 - 12.0*sk2 + 24.0 )*C + 24.0 );
	t4 = -edibkb*sk3*(SL +sk*(C - lambda*CL) - S );
	ck =  -24.0*eta*( t1 + t2 + t3 + t4 )/sk3/sk3;
	struc  = 1./(1.-ck);

	return(struc);
"""

# ER defaults to 0.0
# VR defaults to 1.0

def random():
    pars = dict(
        scale=1, background=0,
        radius_effective=10**np.random.uniform(1, 4.7),
        volfraction=np.random.uniform(0.00001, 0.08),
        welldepth=np.random.uniform(0, 1.5),
        wellwidth=np.random.uniform(1, 1.2),
    )
    return pars

demo = dict(radius_effective=50, volfraction=0.04, welldepth=1.5,
            wellwidth=1.2, radius_effective_pd=0, radius_effective_pd_n=0)
#
tests = [
    [{'scale': 1.0, 'background': 0.0, 'radius_effective': 50.0,
      'volfraction': 0.04, 'welldepth': 1.5, 'wellwidth': 1.2,
      'radius_effective_pd': 0}, [0.001], [0.97665742]],
    ]
# ADDED by: converting from sasview RKH  ON: 16Mar2016
