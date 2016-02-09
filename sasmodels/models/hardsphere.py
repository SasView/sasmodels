# Note: model title and parameter table are inserted automatically
r"""Calculate the interparticle structure factor for monodisperse
spherical particles interacting through hard sphere (excluded volume)
interactions.

The calculation uses the Percus-Yevick closure where the interparticle
potential is

.. math::

    U(r) = \begin{cases}
    \infty & r < 2R \\
    0 & r \geq 2R
    \end{cases}

where $r$ is the distance from the center of the sphere of a radius $R$.

For a 2D plot, the wave transfer is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


.. figure:: img/hardSphere_1d.jpg

    1D plot using the default values (in linear scale).

References
----------

J K Percus, J Yevick, *J. Phys. Rev.*, 110, (1958) 1
"""

from numpy import inf

name = "hardsphere_fish"
title = "Hard sphere structure factor from FISH, with Percus-Yevick closure"
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
      double D,A,B,G,X,X2,X4,S,C,FF,HARDSPH;

      if(fabs(effect_radius) < 1.E-12) {
               HARDSPH=1.0;
               return(HARDSPH);
      }
      D=pow((1.-volfraction),2);
      A=pow((1.+2.*volfraction)/D, 2);
      X=fabs(q*effect_radius*2.0);

      if(X < 5.E-06) {
                 HARDSPH=1./A;
                 return(HARDSPH);
      }
      X2=pow(X,2);
      X4=pow(X2,2);
      B=-6.*volfraction* pow((1.+0.5*volfraction)/D ,2);
      G=0.5*volfraction*A;

      if(X < 0.2) {
      // use Taylor series expansion for small X, IT IS VERY PICKY ABOUT THE X CUT OFF VALUE, ought to be lower in double. 
      // No obvious way to rearrange the equations to avoid needing a very high number of significant figures. 
      // Series expansion found using Mathematica software. Numerical test in .xls showed terms to X^2 are sufficient 
      // for 5 or 6 significant figures but I put the X^4 one in anyway 
            FF = 8*A +6*B + 4*G - (0.8*A +2.0*B/3.0 +0.5*G)*X2 +(A/35. +B/40. +G/50.)*X4;
            // combining the terms makes things worse at smallest Q in single precision
            //FF = (8-0.8*X2)*A +(3.0-X2/3.)*2*B + (4+0.5*X2)*G +(A/35. +B/40. +G/50.)*X4;
            // note that G = -volfraction*A/2, combining this makes no further difference at smallest Q
            //FF = (8 +2.*volfraction + ( volfraction/4. -0.8 +(volfraction/100. -1./35.)*X2 )*X2 )*A  + (3.0 -X2/3. +X4/40)*2*B;
            HARDSPH= 1./(1. + volfraction*FF );
            return(HARDSPH);
      }
      SINCOS(X,S,C);

// RKH Feb 2016, use version from FISH code as it is better than original sasview one at small Q in single precision
      FF=A*(S-X*C)/X + B*(2.*X*S -(X2-2.)*C -2.)/X2 + G*( (4.*X2*X -24.*X)*S -(X4 -12.*X2 +24.)*C +24. )/X4;
      HARDSPH= 1./(1. + 24.*volfraction*FF/X2 );

// rearrange the terms, is now about same as sasmodels
//     FF=A*(S/X-C) + B*(2.*S/X - C +2.0*(C-1.0)/X2) + G*( (4./X -24./X3)*S -(1.0 -12./X2 +24./X4)*C +24./X4 );
//     HARDSPH= 1./(1. + 24.*volfraction*FF/X2 );
// remove 1/X2 from final line, take more powers of X inside the brackets, stil bad
//      FF=A*(S/X3-C/X2) + B*(2.*S/X3 - C/X2 +2.0*(C-1.0)/X4) + G*( (4./X -24./X3)*S -(1.0 -12./X2 +24./X4)*C +24./X4 )/X2;
//      HARDSPH= 1./(1. + 24.*volfraction*FF );
      return(HARDSPH);
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

tests = [
        [ {'scale': 1.0, 'background' : 0.0, 'effect_radius' : 50.0, 'volfraction' : 0.2,
           'effect_radius_pd' : 0}, [0.001], [0.209128]]
        ]

