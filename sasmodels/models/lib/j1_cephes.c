/*							j1.c
 *
 *	Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j1();
 *
 * y = j1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order one of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 24 term Chebyshev
 * expansion is used. In the second, the asymptotic
 * trigonometric representation is employed using two
 * rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       0, 30       10000       4.0e-17     1.1e-17
 *    IEEE      0, 30       30000       2.6e-16     1.1e-16
 *
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

double J1(double x) {

//Cephes double pression function
#if FLOAT_SIZE>4

    double w, z, p, q, xn;

    const double Z1 = 1.46819706421238932572E1;
    const double Z2 = 4.92184563216946036703E1;
    const double THPIO4 =  2.35619449019234492885;
    const double SQ2OPI = 0.79788456080286535588;

    w = x;
    if( x < 0 )
	    w = -x;

    if( w <= 5.0 )
	{
	    z = x * x;
	    w = polevl( z, RPJ1, 3 ) / p1evl( z, RQJ1, 8 );
	    w = w * x * (z - Z1) * (z - Z2);
	    return( w );
	}

    w = 5.0/x;
    z = w * w;

    p = polevl( z, PPJ1, 6)/polevl( z, PQJ1, 6 );
    q = polevl( z, QPJ1, 7)/p1evl( z, QQJ1, 7 );

    xn = x - THPIO4;

    double sn, cn;
    SINCOS(xn, sn, cn);
    p = p * cn - w * q * sn;

    return( p * SQ2OPI / sqrt(x) );


//Single precission version of cephes
#else
    double xx, w, z, p, q, xn;

    const double Z1 = 1.46819706421238932572E1;
    const double THPIO4F =  2.35619449019234492885;    /* 3*pi/4 */


    xx = x;
    if( xx < 0 )
	    xx = -x;

    if( xx <= 2.0 )
	{
	    z = xx * xx;
	    p = (z-Z1) * xx * polevl( z, JPJ1, 4 );
	    return( p );
	}

    q = 1.0/x;
    w = sqrt(q);

    p = w * polevl( q, MO1J1, 7);
    w = q*q;
    xn = q * polevl( w, PH1J1, 7) - THPIO4F;
    p = p * cos(xn + xx);

    return(p);
#endif
}

//Finally J1c function that equals 2*J1(x)/x
double J1c(double x) {
    return (x != 0.0 ) ? 2.0*J1(x)/x : 1.0;
}