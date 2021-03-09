/*!
 * \file
 * \brief Implementation of Airy function
 * \author Tony Ottosson
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 *
 * This is slightly modified routine from the Cephes library:
 * http://www.netlib.org/cephes/
 */




/*
 * Airy function
 *
 * double x, ai, aip, bi, bip;
 * int airy();
 *
 * airy( x, _&ai, _&aip, _&bi, _&bip );
 *
 * DESCRIPTION:
 *
 * Solution of the differential equation
 *
 * y"(x) = xy.
 *
 * The function returns the two independent solutions Ai, Bi
 * and their first derivatives Ai'(x), Bi'(x).
 *
 * Evaluation is by power series summation for small x,
 * by rational minimax approximations for large x.
 *
 * ACCURACY:
 * Error criterion is absolute when function <= 1, relative
 * when function > 1, except * denotes relative error criterion.
 * For large negative x, the absolute error increases as x^1.5.
 * For large positive x, the relative error increases as x^1.5.
 *
 * Arithmetic  domain   function  # trials      peak         rms
 * IEEE        -10, 0     Ai        10000       1.6e-15     2.7e-16
 * IEEE          0, 10    Ai        10000       2.3e-14*    1.8e-15*
 * IEEE        -10, 0     Ai'       10000       4.6e-15     7.6e-16
 * IEEE          0, 10    Ai'       10000       1.8e-14*    1.5e-15*
 * IEEE        -10, 10    Bi        30000       4.2e-15     5.3e-16
 * IEEE        -10, 10    Bi'       30000       4.9e-15     7.3e-16
 */

/*
  Cephes Math Library Release 2.8:  June, 2000
  Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

static double c1 = 0.35502805388781723926;
static double c2 = 0.258819403792806798405;
static double sqrt3 = 1.732050807568877293527;
static double sqpii = 5.64189583547756286948E-1;

#define MAXAIRY 25.77

#define MAXNUM 1.79769313486231570815E308    /* 2**1024*(1-MACHEP) */
#define MACHEP 1.11022302462515654042E-16   /* 2**-53 */

static double AN[8] = {
  3.46538101525629032477E-1,
  1.20075952739645805542E1,
  7.62796053615234516538E1,
  1.68089224934630576269E2,
  1.59756391350164413639E2,
  7.05360906840444183113E1,
  1.40264691163389668864E1,
  9.99999999999999995305E-1,
};
static double AD[8] = {
  5.67594532638770212846E-1,
  1.47562562584847203173E1,
  8.45138970141474626562E1,
  1.77318088145400459522E2,
  1.64234692871529701831E2,
  7.14778400825575695274E1,
  1.40959135607834029598E1,
  1.00000000000000000470E0,
};

static double APN[8] = {
  6.13759184814035759225E-1,
  1.47454670787755323881E1,
  8.20584123476060982430E1,
  1.71184781360976385540E2,
  1.59317847137141783523E2,
  6.99778599330103016170E1,
  1.39470856980481566958E1,
  1.00000000000000000550E0,
};
static double APD[8] = {
  3.34203677749736953049E-1,
  1.11810297306158156705E1,
  7.11727352147859965283E1,
  1.58778084372838313640E2,
  1.53206427475809220834E2,
  6.86752304592780337944E1,
  1.38498634758259442477E1,
  9.99999999999999994502E-1,
};


static double BN16[5] = {
  -2.53240795869364152689E-1,
  5.75285167332467384228E-1,
  -3.29907036873225371650E-1,
  6.44404068948199951727E-2,
  -3.82519546641336734394E-3,
};
static double BD16[5] = {
  /* 1.00000000000000000000E0,*/
  -7.15685095054035237902E0,
  1.06039580715664694291E1,
  -5.23246636471251500874E0,
  9.57395864378383833152E-1,
  -5.50828147163549611107E-2,
};


static double BPPN[5] = {
  4.65461162774651610328E-1,
  -1.08992173800493920734E0,
  6.38800117371827987759E-1,
  -1.26844349553102907034E-1,
  7.62487844342109852105E-3,
};
static double BPPD[5] = {
  /* 1.00000000000000000000E0,*/
  -8.70622787633159124240E0,
  1.38993162704553213172E1,
  -7.14116144616431159572E0,
  1.34008595960680518666E0,
  -7.84273211323341930448E-2,
};

static double AFN[9] = {
  -1.31696323418331795333E-1,
  -6.26456544431912369773E-1,
  -6.93158036036933542233E-1,
  -2.79779981545119124951E-1,
  -4.91900132609500318020E-2,
  -4.06265923594885404393E-3,
  -1.59276496239262096340E-4,
  -2.77649108155232920844E-6,
  -1.67787698489114633780E-8,
};
static double AFD[9] = {
  /* 1.00000000000000000000E0,*/
  1.33560420706553243746E1,
  3.26825032795224613948E1,
  2.67367040941499554804E1,
  9.18707402907259625840E0,
  1.47529146771666414581E0,
  1.15687173795188044134E-1,
  4.40291641615211203805E-3,
  7.54720348287414296618E-5,
  4.51850092970580378464E-7,
};

static double AGN[11] = {
  1.97339932091685679179E-2,
  3.91103029615688277255E-1,
  1.06579897599595591108E0,
  9.39169229816650230044E-1,
  3.51465656105547619242E-1,
  6.33888919628925490927E-2,
  5.85804113048388458567E-3,
  2.82851600836737019778E-4,
  6.98793669997260967291E-6,
  8.11789239554389293311E-8,
  3.41551784765923618484E-10,
};
static double AGD[10] = {
  /*  1.00000000000000000000E0,*/
  9.30892908077441974853E0,
  1.98352928718312140417E1,
  1.55646628932864612953E1,
  5.47686069422975497931E0,
  9.54293611618961883998E-1,
  8.64580826352392193095E-2,
  4.12656523824222607191E-3,
  1.01259085116509135510E-4,
  1.17166733214413521882E-6,
  4.91834570062930015649E-9,
};


static double APFN[9] = {
  1.85365624022535566142E-1,
  8.86712188052584095637E-1,
  9.87391981747398547272E-1,
  4.01241082318003734092E-1,
  7.10304926289631174579E-2,
  5.90618657995661810071E-3,
  2.33051409401776799569E-4,
  4.08718778289035454598E-6,
  2.48379932900442457853E-8,
};
static double APFD[9] = {
  /*  1.00000000000000000000E0,*/
  1.47345854687502542552E1,
  3.75423933435489594466E1,
  3.14657751203046424330E1,
  1.09969125207298778536E1,
  1.78885054766999417817E0,
  1.41733275753662636873E-1,
  5.44066067017226003627E-3,
  9.39421290654511171663E-5,
  5.65978713036027009243E-7,
};

static double APGN[11] = {
  -3.55615429033082288335E-2,
  -6.37311518129435504426E-1,
  -1.70856738884312371053E0,
  -1.50221872117316635393E0,
  -5.63606665822102676611E-1,
  -1.02101031120216891789E-1,
  -9.48396695961445269093E-3,
  -4.60325307486780994357E-4,
  -1.14300836484517375919E-5,
  -1.33415518685547420648E-7,
  -5.63803833958893494476E-10,
};
static double APGD[11] = {
  /*  1.00000000000000000000E0,*/
  9.85865801696130355144E0,
  2.16401867356585941885E1,
  1.73130776389749389525E1,
  6.17872175280828766327E0,
  1.08848694396321495475E0,
  9.95005543440888479402E-2,
  4.78468199683886610842E-3,
  1.18159633322838625562E-4,
  1.37480673554219441465E-6,
  5.79912514929147598821E-9,
};


int airy(double x, double *ai, double *aip, double *bi, double *bip)
{
  double z, zz, t, f, g, uf, ug, k, zeta, theta;
  int domflg;

  domflg = 0;
  if (x > MAXAIRY) {
    *ai = 0;
    *aip = 0;
    *bi = MAXNUM;
    *bip = MAXNUM;
    return(-1);
  }

  if (x < -2.09) {
    domflg = 15;
    t = sqrt(-x);
    zeta = -2.0 * x * t / 3.0;
    t = sqrt(t);
    k = sqpii / t;
    z = 1.0 / zeta;
    zz = z * z;
    uf = 1.0 + zz * polevl(zz, AFN, 8) / p1evl(zz, AFD, 9);
    ug = z * polevl(zz, AGN, 10) / p1evl(zz, AGD, 10);
    theta = zeta + 0.25 * M_PI;
    f = sin(theta);
    g = cos(theta);
    *ai = k * (f * uf - g * ug);
    *bi = k * (g * uf + f * ug);
    uf = 1.0 + zz * polevl(zz, APFN, 8) / p1evl(zz, APFD, 9);
    ug = z * polevl(zz, APGN, 10) / p1evl(zz, APGD, 10);
    k = sqpii * t;
    *aip = -k * (g * uf + f * ug);
    *bip = k * (f * uf - g * ug);
    return(0);
  }

  if (x >= 2.09) { /* cbrt(9) */
    domflg = 5;
    t = sqrt(x);
    zeta = 2.0 * x * t / 3.0;
    g = exp(zeta);
    t = sqrt(t);
    k = 2.0 * t * g;
    z = 1.0 / zeta;
    f = polevl(z, AN, 7) / polevl(z, AD, 7);
    *ai = sqpii * f / k;
    k = -0.5 * sqpii * t / g;
    f = polevl(z, APN, 7) / polevl(z, APD, 7);
    *aip = f * k;

    if (x > 8.3203353) { /* zeta > 16 */
      f = z * polevl(z, BN16, 4) / p1evl(z, BD16, 5);
      k = sqpii * g;
      *bi = k * (1.0 + f) / t;
      f = z * polevl(z, BPPN, 4) / p1evl(z, BPPD, 5);
      *bip = k * t * (1.0 + f);
      return(0);
    }
  }

  f = 1.0;
  g = x;
  t = 1.0;
  uf = 1.0;
  ug = x;
  k = 1.0;
  z = x * x * x;
  while (t > MACHEP) {
    uf *= z;
    k += 1.0;
    uf /= k;
    ug *= z;
    k += 1.0;
    ug /= k;
    uf /= k;
    f += uf;
    k += 1.0;
    ug /= k;
    g += ug;
    t = fabs(uf / f);
  }
  uf = c1 * f;
  ug = c2 * g;
  if ((domflg & 1) == 0)
    *ai = uf - ug;
  if ((domflg & 2) == 0)
    *bi = sqrt3 * (uf + ug);

  /* the deriviative of ai */
  k = 4.0;
  uf = x * x / 2.0;
  ug = z / 3.0;
  f = uf;
  g = 1.0 + ug;
  uf /= 3.0;
  t = 1.0;

  while (t > MACHEP) {
    uf *= z;
    ug /= k;
    k += 1.0;
    ug *= z;
    uf /= k;
    f += uf;
    k += 1.0;
    ug /= k;
    uf /= k;
    g += ug;
    k += 1.0;
    t = fabs(ug / g);
  }

  uf = c1 * f;
  ug = c2 * g;
  if ((domflg & 4) == 0)
    *aip = uf - ug;
  if ((domflg & 8) == 0)
    *bip = sqrt3 * (uf + ug);
  return(0);
}




/*!
 * \file
 * \brief Implementation of gamma functions
 * \author Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 *
 * This is slightly modified routine from the Cephes library:
 * http://www.netlib.org/cephes/
 */



/*
 * Gamma function
 *
 *
 * SYNOPSIS:
 *
 * double x, y, gam();
 * extern int sgngam;
 *
 * y = gam( x );
 *
 *
 * DESCRIPTION:
 *
 * Returns gamma function of the argument.  The result is
 * correctly signed, and the sign (+1 or -1) is also
 * returned in a global (extern) variable named sgngam.
 * This variable is also filled in by the logarithmic gamma
 * function lgam().
 *
 * Arguments |x| <= 34 are reduced by recurrence and the function
 * approximated by a rational function of degree 6/7 in the
 * interval (2,3).  Large arguments are handled by Stirling's
 * formula. Large negative arguments are made positive using
 * a reflection formula.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE    -170,-33      20000       2.3e-15     3.3e-16
 *    IEEE     -33,  33     20000       9.4e-16     2.2e-16
 *    IEEE      33, 171.6   20000       2.3e-15     3.2e-16
 *
 * Error for arguments outside the test range will be larger
 * owing to error amplification by the exponential function.
 */

/*
 * Natural logarithm of gamma function
 *
 *
 * SYNOPSIS:
 *
 * double x, y, lgam();
 * extern int sgngam;
 *
 * y = lgam( x );
 *
 *
 * DESCRIPTION:
 *
 * Returns the base e (2.718...) logarithm of the absolute
 * value of the gamma function of the argument.
 * The sign (+1 or -1) of the gamma function is returned in a
 * global (extern) variable named sgngam.
 *
 * For arguments greater than 13, the logarithm of the gamma
 * function is approximated by the logarithmic version of
 * Stirling's formula using a polynomial approximation of
 * degree 4. Arguments between -33 and +33 are reduced by
 * recurrence to the interval [2,3] of a rational approximation.
 * The cosecant reflection formula is employed for arguments
 * less than -33.
 *
 * Arguments greater than MAXLGM return INFINITY and an error
 * message.  MAXLGM = 2.556348e305 for IEEE arithmetic.
 *
 *
 * ACCURACY:
 *
 * arithmetic      domain        # trials     peak         rms
 *    IEEE    0, 3                 28000     5.4e-16     1.1e-16
 *    IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
 * The error criterion was relative when the function magnitude
 * was greater than one but absolute when it was less than one.
 *
 * The following test used the relative error criterion, though
 * at certain points the relative error could be much higher than
 * indicated.
 *    IEEE    -200, -4             10000     4.8e-16     1.3e-16
 */

/*
  Cephes Math Library Release 2.8:  June, 2000
  Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
*/

#ifndef INFINITY
#  define INFINITY 1.79769313486231570815E308 /* 2**1024*(1-MACHEP) */
#endif
#ifndef NAN
#  define NAN 0.0
#endif

static double P[] = {
  1.60119522476751861407E-4,
  1.19135147006586384913E-3,
  1.04213797561761569935E-2,
  4.76367800457137231464E-2,
  2.07448227648435975150E-1,
  4.94214826801497100753E-1,
  9.99999999999999996796E-1
};
static double Q[] = {
  -2.31581873324120129819E-5,
  5.39605580493303397842E-4,
  -4.45641913851797240494E-3,
  1.18139785222060435552E-2,
  3.58236398605498653373E-2,
  -2.34591795718243348568E-1,
  7.14304917030273074085E-2,
  1.00000000000000000320E0
};
static double LOGPI = 1.14472988584940017414;

/* Stirling's formula for the gamma function */
static double STIR[5] = {
  7.87311395793093628397E-4,
  -2.29549961613378126380E-4,
  -2.68132617805781232825E-3,
  3.47222221605458667310E-3,
  8.33333333333482257126E-2,
};
static double MAXSTIR = 143.01608;
static double SQTPI = 2.50662827463100050242E0;
static double MAXLGM = 2.556348e305;

int sgngam = 0;

/*!
 * \brief Gamma function computed by Stirling's formula.
 * The polynomial STIR is valid for 33 <= x <= 172.
 */
static double stirf(double x)
{
  double y, w, v;

  w = 1.0 / x;
  w = 1.0 + w * polevl(w, STIR, 4);
  y = exp(x);
  if (x > MAXSTIR) { /* Avoid overflow in pow() */
    v = pow(x, 0.5 * x - 0.25);
    y = v * (v / y);
  }
  else {
    y = pow(x, x - 0.5) / y;
  }
  y = SQTPI * y * w;
  return(y);
}



double gam(double x)
{
  double p, q, z;
  int i;

  sgngam = 1;
  if (isnan(x))
    return(x);

  if (isinf(x) == 1)
    return(x);
  if (isinf(x) == -1)
    return(NAN);

  q = fabs(x);

  if (q > 33.0) {
    if (x < 0.0) {
      p = floor(q);
      if (p == q) {
      gamnan:
        
        return (NAN);
      }
      i = (int)(p);
      if ((i & 1) == 0)
        sgngam = -1;
      z = q - p;
      if (z > 0.5) {
        p += 1.0;
        z = q - p;
      }
      z = q * sin(M_PI * z);
      if (z == 0.0) {
        return(sgngam * INFINITY);
      }
      z = fabs(z);
      z = M_PI/ (z * stirf(q));
    }
    else {
      z = stirf(x);
    }
    return(sgngam * z);
  }

  z = 1.0;
  while (x >= 3.0) {
    x -= 1.0;
    z *= x;
  }

  while (x < 0.0) {
    if (x > -1.E-9)
      goto small;
    z /= x;
    x += 1.0;
  }

  while (x < 2.0) {
    if (x < 1.e-9)
      goto small;
    z /= x;
    x += 1.0;
  }

  if (x == 2.0)
    return(z);

  x -= 2.0;
  p = polevl(x, P, 6);
  q = polevl(x, Q, 7);
  return(z * p / q);

small:
  if (x == 0.0) {
    goto gamnan;
  }
  else
    return(z / ((1.0 + 0.5772156649015329 * x) * x));
}



/* A[]: Stirling's formula expansion of log gamma
 * B[], C[]: log gamma function between 2 and 3
 */
static double A[] = {
  8.11614167470508450300E-4,
  -5.95061904284301438324E-4,
  7.93650340457716943945E-4,
  -2.77777777730099687205E-3,
  8.33333333333331927722E-2
};
static double B[] = {
  -1.37825152569120859100E3,
  -3.88016315134637840924E4,
  -3.31612992738871184744E5,
  -1.16237097492762307383E6,
  -1.72173700820839662146E6,
  -8.53555664245765465627E5
};
static double C[] = {
  /* 1.00000000000000000000E0, */
  -3.51815701436523470549E2,
  -1.70642106651881159223E4,
  -2.20528590553854454839E5,
  -1.13933444367982507207E6,
  -2.53252307177582951285E6,
  -2.01889141433532773231E6
};
/* log( sqrt( 2*pi ) ) */
static double LS2PI  =  0.91893853320467274178;


/*! Logarithm of gamma function */
double lgam(double x)
{
  double p, q, u, w, z;
  int i;

  sgngam = 1;
  if (isnan(x))
    return(x);

  if (isinf(x))
    return(INFINITY);

  if (x < -34.0) {
    q = -x;
    w = lgam(q); /* note this modifies sgngam! */
    p = floor(q);
    if (p == q) {
    lgsing:
      
      return (INFINITY);
    }
    i = (int)(p);
    if ((i & 1) == 0)
      sgngam = -1;
    else
      sgngam = 1;
    z = q - p;
    if (z > 0.5) {
      p += 1.0;
      z = p - q;
    }
    z = q * sin(M_PI * z);
    if (z == 0.0)
      goto lgsing;
    /*      z = log(PI) - log( z ) - w;*/
    z = LOGPI - log(z) - w;
    return(z);
  }

  if (x < 13.0) {
    z = 1.0;
    p = 0.0;
    u = x;
    while (u >= 3.0) {
      p -= 1.0;
      u = x + p;
      z *= u;
    }
    while (u < 2.0) {
      if (u == 0.0)
        goto lgsing;
      z /= u;
      p += 1.0;
      u = x + p;
    }
    if (z < 0.0) {
      sgngam = -1;
      z = -z;
    }
    else
      sgngam = 1;
    if (u == 2.0)
      return(log(z));
    p -= 2.0;
    x = x + p;
    p = x * polevl(x, B, 5) / p1evl(x, C, 6);
    return(log(z) + p);
  }

  if (x > MAXLGM) {
    return(sgngam * INFINITY);
  }

  q = (x - 0.5) * log(x) - x + LS2PI;
  if (x > 1.0e8)
    return(q);

  p = 1.0 / (x * x);
  if (x >= 1000.0)
    q += ((7.9365079365079365079365e-4 * p
           - 2.7777777777777777777778e-3) * p
          + 0.0833333333333333333333) / x;
  else
    q += polevl(p, A, 4) / x;
  return(q);
}






/*!
 * \file
 * \brief Implementation of Bessel functions of noninteager order
 * \author Tony Ottosson
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 *
 * This is slightly modified routine from the Cephes library:
 * http://www.netlib.org/cephes/
 */



/*
 * Bessel function of noninteger order
 *
 * double v, x, y, jv();
 *
 * y = jv( v, x );
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order v of the argument,
 * where v is real.  Negative x is allowed if v is an integer.
 *
 * Several expansions are included: the ascending power
 * series, the Hankel expansion, and two transitional
 * expansions for large v.  If v is not too large, it
 * is reduced by recurrence to a region of best accuracy.
 * The transitional expansions give 12D accuracy for v > 500.
 *
 * ACCURACY:
 * Results for integer v are indicated by *, where x and v
 * both vary from -125 to +125.  Otherwise,
 * x ranges from 0 to 125, v ranges as indicated by "domain."
 * Error criterion is absolute, except relative when |jv()| > 1.
 *
 * arithmetic  v domain  x domain    # trials      peak       rms
 *    IEEE      0,125     0,125      100000      4.6e-15    2.2e-16
 *    IEEE   -125,0       0,125       40000      5.4e-11    3.7e-13
 *    IEEE      0,500     0,500       20000      4.4e-15    4.0e-16
 * Integer v:
 *    IEEE   -125,125   -125,125      50000      3.5e-15*   1.9e-16*
 */

/*
  Cephes Math Library Release 2.8:  June, 2000
  Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
*/

static double recur(double *, double, double *, int);
static double jvs(double, double);
static double hankel(double, double);
static double jnx(double, double);
static double jnt(double, double);

#define MAXGAM 171.624376956302725


#define MAXLOG 7.08396418532264106224E2     /* log 2**1022 */
#define MINLOG -7.08396418532264106224E2    /* log 2**-1022 */


#define BIG  1.44115188075855872E+17



// ---------------------------- jv() -------------------------------------------------------
double jv(double n, double x)
{
  double k, q, t, y, an;
  int i, sign, nint;

  nint = 0; /* Flag for integer n */
  sign = 1; /* Flag for sign inversion */
  an = fabs(n);
  y = floor(an);
  if (y == an) {
    nint = 1;
    i = (int)(an - 16384.0 * floor(an / 16384.0));
    if (n < 0.0) {
      if (i & 1)
        sign = -sign;
      n = an;
    }
    if (x < 0.0) {
      if (i & 1)
        sign = -sign;
      x = -x;
    }
    if (n == 0.0)  // use 0th order bessel function
      return sas_J0(x);
    if (n == 1.0)  // use 1th order bessel function
      return (sign * sas_J1(x));
  }



  y = fabs(x);

  if (y < MACHEP)
    goto underf;

  k = 3.6 * sqrt(y);
  t = 3.6 * sqrt(an);
  if ((y < t) && (an > 21.0))
    return(sign * jvs(n, x));
  if ((an < k) && (y > 21.0))
    return(sign * hankel(n, x));

  if (an < 500.0) {
    /* Note: if x is too large, the continued
     * fraction will fail; but then the
     * Hankel expansion can be used.
     */
    if (nint != 0) {
      k = 0.0;
      q = recur(&n, x, &k, 1);
      if (k == 0.0) {
        y =sas_J0(x)/ q;
        goto done;
      }
      if (k == 1.0) {
        y = sas_J1(x) / q;
        goto done;
      }
    }

    if (an > 2.0 * y)
      goto rlarger;

    if ((n >= 0.0) && (n < 20.0)
        && (y > 6.0) && (y < 20.0)) {
      /* Recur backwards from a larger value of n
       */
    rlarger:
      k = n;

      y = y + an + 1.0;
      if (y < 30.0)
        y = 30.0;
      y = n + floor(y - n);
      q = recur(&y, x, &k, 0);
      y = jvs(y, x) * q;
      goto done;
    }

    if (k <= 30.0) {
      k = 2.0;
    }
    else if (k < 90.0) {
      k = (3 * k) / 4;
    }
    if (an > (k + 3.0)) {
      if (n < 0.0)
        k = -k;
      q = n - floor(n);
      k = floor(k) + q;
      if (n > 0.0)
        q = recur(&n, x, &k, 1);
      else {
        t = k;
        k = n;
        q = recur(&t, x, &k, 1);
        k = t;
      }
      if (q == 0.0) {
      underf:
        y = 0.0;
        goto done;
      }
    }
    else {
      k = n;
      q = 1.0;
    }

    /* boundary between convergence of
     * power series and Hankel expansion
     */
    y = fabs(k);
    if (y < 26.0)
      t = (0.0083 * y + 0.09) * y + 12.9;
    else
      t = 0.9 * y;

    if (x > t)
      y = hankel(k, x);
    else
      y = jvs(k, x);

    if (n > 0.0)
      y /= q;
    else
      y *= q;
  }

  else {
    /* For large n, use the uniform expansion
     * or the transitional expansion.
     * But if x is of the order of n**2,
     * these may blow up, whereas the
     * Hankel expansion will then work.
     */

    if (n < 0.0) {
      
      y = 0.0;
      goto done;
    }
    t = x / n;
    t /= n;
    if (t > 0.3)
      y = hankel(n, x);
    else
      y = jnx(n, x);
  }

done:
  return(sign * y);
}

/* Reduce the order by backward recurrence.
 * AMS55 #9.1.27 and 9.1.73.
 */

static double recur(double *n, double x, double *newn, int cancel)
{
  double pkm2, pkm1, pk, qkm2, qkm1;
  /* double pkp1; */
  double k, ans, qk, xk, yk, r, t, kf;
  static double big = BIG;
  int nflag, ctr;

  /* continued fraction for Jn(x)/Jn-1(x)  */
  if (*n < 0.0)
    nflag = 1;
  else
    nflag = 0;

fstart:
  pkm2 = 0.0;
  qkm2 = 1.0;
  pkm1 = x;
  qkm1 = *n + *n;
  xk = -x * x;
  yk = qkm1;
  ans = 1.0;
  ctr = 0;
  do {
    yk += 2.0;
    pk = pkm1 * yk +  pkm2 * xk;
    qk = qkm1 * yk +  qkm2 * xk;
    pkm2 = pkm1;
    pkm1 = pk;
    qkm2 = qkm1;
    qkm1 = qk;
    if (qk != 0)
      r = pk / qk;
    else
      r = 0.0;
    if (r != 0) {
      t = fabs((ans - r) / r);
      ans = r;
    }
    else
      t = 1.0;

    if (++ctr > 1000) {
      
      goto done;
    }

    if (t < MACHEP)
      goto done;

    if (fabs(pk) > big) {
      pkm2 /= big;
      pkm1 /= big;
      qkm2 /= big;
      qkm1 /= big;
    }
  }
  while (t > MACHEP);

done:

  /* Change n to n-1 if n < 0 and the continued fraction is small
   */
  if (nflag > 0) {
    if (fabs(ans) < 0.125) {
      nflag = -1;
      *n = *n - 1.0;
      goto fstart;
    }
  }


  kf = *newn;

  /* backward recurrence
   *              2k
   *  J   (x)  =  --- J (x)  -  J   (x)
   *   k-1         x   k         k+1
   */

  pk = 1.0;
  pkm1 = 1.0 / ans;
  k = *n - 1.0;
  r = 2 * k;
  do {
    pkm2 = (pkm1 * r  -  pk * x) / x;
    /* pkp1 = pk; */
    pk = pkm1;
    pkm1 = pkm2;
    r -= 2.0;
    /*
    t = fabs(pkp1) + fabs(pk);
    if( (k > (kf + 2.5)) && (fabs(pkm1) < 0.25*t) )
    {
    k -= 1.0;
    t = x*x;
    pkm2 = ( (r*(r+2.0)-t)*pk - r*x*pkp1 )/t;
    pkp1 = pk;
    pk = pkm1;
    pkm1 = pkm2;
    r -= 2.0;
    }
    */
    k -= 1.0;
  }
  while (k > (kf + 0.5));

  /* Take the larger of the last two iterates
   * on the theory that it may have less cancellation error.
   */

  if (cancel) {
    if ((kf >= 0.0) && (fabs(pk) > fabs(pkm1))) {
      k += 1.0;
      pkm2 = pk;
    }
  }
  *newn = k;

  return(pkm2);
}



/* Ascending power series for Jv(x).
 * AMS55 #9.1.10.
 */


static double jvs(double n, double x)
{
  double t, u, y, z, k;
  int ex;

  z = -x * x / 4.0;
  u = 1.0;
  y = u;
  k = 1.0;
  t = 1.0;

  while (t > MACHEP) {
    u *= z / (k * (n + k));
    y += u;
    k += 1.0;
    if (y != 0)
      t = fabs(u / y);
  }

  t = frexp(0.5 * x, &ex);
  ex = (int)(ex * n);
  if ((ex > -1023)
      && (ex < 1023)
      && (n > 0.0)
      && (n < (MAXGAM - 1.0))) {
    t = pow(0.5 * x, n) / gam(n + 1.0);

    y *= t;
  }
  else {
    t = n * log(0.5 * x) - lgam(n + 1.0);
    if (y < 0) {
      sgngam = -sgngam;
      y = -y;
    }
    t += log(y);

    if (t < -MAXLOG) {
      return(0.0);
    }

    if (t > MAXLOG) {
     
      return(MAXNUM);
    }
    y = sgngam * exp(t);
  }
  return(y);
}

/* Hankel's asymptotic expansion
 * for large x.
 * AMS55 #9.2.5.
 */
static double hankel(double n, double x)
{
  double t, u, z, k, sign, conv;
  double p, q, j, m, pp, qq;
  int flag;

  m = 4.0 * n * n;
  j = 1.0;
  z = 8.0 * x;
  k = 1.0;
  p = 1.0;
  u = (m - 1.0) / z;
  q = u;
  sign = 1.0;
  conv = 1.0;
  flag = 0;
  t = 1.0;
  pp = 1.0e38;
  qq = 1.0e38;

  while (t > MACHEP) {
    k += 2.0;
    j += 1.0;
    sign = -sign;
    u *= (m - k * k) / (j * z);
    p += sign * u;
    k += 2.0;
    j += 1.0;
    u *= (m - k * k) / (j * z);
    q += sign * u;
    t = fabs(u / p);
    if (t < conv) {
      conv = t;
      qq = q;
      pp = p;
      flag = 1;
    }
    /* stop if the terms start getting larger */
    if ((flag != 0) && (t > conv)) {
      goto hank1;
    }
  }

hank1:
  u = x - (0.5 * n + 0.25) * M_PI;
  t = sqrt(2.0 / (M_PI * x)) * (pp * cos(u) - qq * sin(u));

  return(t);
}


/* Asymptotic expansion for large n.
 * AMS55 #9.3.35.
 */

static double lambda[] = {
  1.0,
  1.041666666666666666666667E-1,
  8.355034722222222222222222E-2,
  1.282265745563271604938272E-1,
  2.918490264641404642489712E-1,
  8.816272674437576524187671E-1,
  3.321408281862767544702647E+0,
  1.499576298686255465867237E+1,
  7.892301301158651813848139E+1,
  4.744515388682643231611949E+2,
  3.207490090890661934704328E+3
};
static double mu[] = {
  1.0,
  -1.458333333333333333333333E-1,
  -9.874131944444444444444444E-2,
  -1.433120539158950617283951E-1,
  -3.172272026784135480967078E-1,
  -9.424291479571202491373028E-1,
  -3.511203040826354261542798E+0,
  -1.572726362036804512982712E+1,
  -8.228143909718594444224656E+1,
  -4.923553705236705240352022E+2,
  -3.316218568547972508762102E+3
};
static double P1[] = {
  -2.083333333333333333333333E-1,
  1.250000000000000000000000E-1
};
static double P2[] = {
  3.342013888888888888888889E-1,
  -4.010416666666666666666667E-1,
  7.031250000000000000000000E-2
};
static double P3[] = {
  -1.025812596450617283950617E+0,
  1.846462673611111111111111E+0,
  -8.912109375000000000000000E-1,
  7.324218750000000000000000E-2
};
static double P4[] = {
  4.669584423426247427983539E+0,
  -1.120700261622299382716049E+1,
  8.789123535156250000000000E+0,
  -2.364086914062500000000000E+0,
  1.121520996093750000000000E-1
};
static double P5[] = {
  -2.8212072558200244877E1,
  8.4636217674600734632E1,
  -9.1818241543240017361E1,
  4.2534998745388454861E1,
  -7.3687943594796316964E0,
  2.27108001708984375E-1
};
static double P6[] = {
  2.1257013003921712286E2,
  -7.6525246814118164230E2,
  1.0599904525279998779E3,
  -6.9957962737613254123E2,
  2.1819051174421159048E2,
  -2.6491430486951555525E1,
  5.7250142097473144531E-1
};
static double P7[] = {
  -1.9194576623184069963E3,
  8.0617221817373093845E3,
  -1.3586550006434137439E4,
  1.1655393336864533248E4,
  -5.3056469786134031084E3,
  1.2009029132163524628E3,
  -1.0809091978839465550E2,
  1.7277275025844573975E0
};


static double jnx(double n, double x)
{
  double zeta, sqz, zz, zp, np;
  double cbn, n23, t, z, sz;
  double pp, qq, z32i, zzi;
  double ak, bk, akl, bkl;
  int sign, doa, dob, nflg, k, s, tk, tkp1, m;
  static double u[8];
  static double ai, aip, bi, bip;

  /* Test for x very close to n.
   * Use expansion for transition region if so.
   */
  cbn = cbrt(n);
  z = (x - n) / cbn;
  if (fabs(z) <= 0.7)
    return(jnt(n, x));

  z = x / n;
  zz = 1.0 - z * z;
  if (zz == 0.0)
    return(0.0);

  if (zz > 0.0) {
    sz = sqrt(zz);
    t = 1.5 * (log((1.0 + sz) / z) - sz); /* zeta ** 3/2  */
    zeta = cbrt(t * t);
    nflg = 1;
  }
  else {
    sz = sqrt(-zz);
    t = 1.5 * (sz - acos(1.0 / z));
    zeta = -cbrt(t * t);
    nflg = -1;
  }
  z32i = fabs(1.0 / t);
  sqz = cbrt(t);

  /* Airy function */
  n23 = cbrt(n * n);
  t = n23 * zeta;

  airy(t, &ai, &aip, &bi, &bip);

  /* polynomials in expansion */
  u[0] = 1.0;
  zzi = 1.0 / zz;
  u[1] = polevl(zzi, P1, 1) / sz;
  u[2] = polevl(zzi, P2, 2) / zz;
  u[3] = polevl(zzi, P3, 3) / (sz * zz);
  pp = zz * zz;
  u[4] = polevl(zzi, P4, 4) / pp;
  u[5] = polevl(zzi, P5, 5) / (pp * sz);
  pp *= zz;
  u[6] = polevl(zzi, P6, 6) / pp;
  u[7] = polevl(zzi, P7, 7) / (pp * sz);

  pp = 0.0;
  qq = 0.0;
  np = 1.0;
  /* flags to stop when terms get larger */
  doa = 1;
  dob = 1;
  akl = MAXNUM;
  bkl = MAXNUM;

  for (k = 0; k <= 3; k++) {
    tk = 2 * k;
    tkp1 = tk + 1;
    zp = 1.0;
    ak = 0.0;
    bk = 0.0;
    for (s = 0; s <= tk; s++) {
      if (doa) {
        if ((s & 3) > 1)
          sign = nflg;
        else
          sign = 1;
        ak += sign * mu[s] * zp * u[tk-s];
      }

      if (dob) {
        m = tkp1 - s;
        if (((m + 1) & 3) > 1)
          sign = nflg;
        else
          sign = 1;
        bk += sign * lambda[s] * zp * u[m];
      }
      zp *= z32i;
    }

    if (doa) {
      ak *= np;
      t = fabs(ak);
      if (t < akl) {
        akl = t;
        pp += ak;
      }
      else
        doa = 0;
    }

    if (dob) {
      bk += lambda[tkp1] * zp * u[0];
      bk *= -np / sqz;
      t = fabs(bk);
      if (t < bkl) {
        bkl = t;
        qq += bk;
      }
      else
        dob = 0;
    }

    if (np < MACHEP)
      break;
    np /= n * n;
  }

  /* normalizing factor ( 4*zeta/(1 - z**2) )**1/4 */
  t = 4.0 * zeta / zz;
  t = sqrt(sqrt(t));

  t *= ai * pp / cbrt(n)  +  aip * qq / (n23 * n);
  return(t);
}

/* Asymptotic expansion for transition region,
 * n large and x close to n.
 * AMS55 #9.3.23.
 */

static double PF2[] = {
  -9.0000000000000000000e-2,
  8.5714285714285714286e-2
};
static double PF3[] = {
  1.3671428571428571429e-1,
  -5.4920634920634920635e-2,
  -4.4444444444444444444e-3
};
static double PF4[] = {
  1.3500000000000000000e-3,
  -1.6036054421768707483e-1,
  4.2590187590187590188e-2,
  2.7330447330447330447e-3
};
static double PG1[] = {
  -2.4285714285714285714e-1,
  1.4285714285714285714e-2
};
static double PG2[] = {
  -9.0000000000000000000e-3,
  1.9396825396825396825e-1,
  -1.1746031746031746032e-2
};
static double PG3[] = {
  1.9607142857142857143e-2,
  -1.5983694083694083694e-1,
  6.3838383838383838384e-3
};


static double jnt(double n, double x)
{
  double z, zz, z3;
  double cbn, n23, cbtwo;
  double ai, aip, bi, bip; /* Airy functions */
  double nk, fk, gk, pp, qq;
  double F[5], G[4];
  int k;

  cbn = cbrt(n);
  z = (x - n) / cbn;
  cbtwo = cbrt(2.0);

  /* Airy function */
  zz = -cbtwo * z;
  airy(zz, &ai, &aip, &bi, &bip);

  /* polynomials in expansion */
  zz = z * z;
  z3 = zz * z;
  F[0] = 1.0;
  F[1] = -z / 5.0;
  F[2] = polevl(z3, PF2, 1) * zz;
  F[3] = polevl(z3, PF3, 2);
  F[4] = polevl(z3, PF4, 3) * z;
  G[0] = 0.3 * zz;
  G[1] = polevl(z3, PG1, 1);
  G[2] = polevl(z3, PG2, 2) * z;
  G[3] = polevl(z3, PG3, 2) * zz;

  pp = 0.0;
  qq = 0.0;
  nk = 1.0;
  n23 = cbrt(n * n);

  for (k = 0; k <= 4; k++) {
    fk = F[k] * nk;
    pp += fk;
    if (k != 4) {
      gk = G[k] * nk;
      qq += gk;
    }

    nk /= n23;
  }

  fk = cbtwo * ai * pp / cbn  +  cbrt(4.0) * aip * qq / n;
  return(fk);
}
