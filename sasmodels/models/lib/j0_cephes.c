/*							j0.c
 *
 *	Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j0();
 *
 * y = j0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order zero of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval the following rational
 * approximation is used:
 *
 *
 *        2         2
 * (w - r  ) (w - r  ) P (w) / Q (w)
 *       1         2    3       8
 *
 *            2
 * where w = x  and the two r's are zeros of the function.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30       10000       4.4e-17     6.3e-18
 *    IEEE      0, 30       60000       4.2e-16     1.1e-16
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

/* Note: all coefficients satisfy the relative error criterion
 * except YP, YQ which are designed for absolute error. */

double J0(double x) {

//Cephes single precission
#if FLOAT_SIZE>4
    double w, z, p, q, xn;

    //const double TWOOPI = 6.36619772367581343075535E-1;
    const double SQ2OPI = 7.9788456080286535587989E-1;
    const double PIO4 = 7.85398163397448309616E-1;

    const double DR1 = 5.78318596294678452118E0;
    const double DR2 = 3.04712623436620863991E1;


    if( x < 0 )
	    x = -x;

    if( x <= 5.0 ) {
	    z = x * x;
	    if( x < 1.0e-5 )
		    return( 1.0 - z/4.0 );

	    p = (z - DR1) * (z - DR2);
	    p = p * polevl( z, RPJ0, 3)/p1evl( z, RQJ0, 8 );
	    return( p );
	}

    w = 5.0/x;
    q = 25.0/(x*x);
    p = polevl( q, PPJ0, 6)/polevl( q, PQJ0, 6 );
    q = polevl( q, QPJ0, 7)/p1evl( q, QQJ0, 7 );
    xn = x - PIO4;

    double sn, cn;
    SINCOS(xn, sn, cn);
    p = p * cn - w * q * sn;

    return( p * SQ2OPI / sqrt(x) );
//Cephes single precission
#else
    double xx, w, z, p, q, xn;

    //const double YZ1 =  0.43221455686510834878;
    //const double YZ2 = 22.401876406482861405;
    //const double YZ3 = 64.130620282338755553;
    const double DR1 =  5.78318596294678452118;
    const double PIO4F = 0.7853981633974483096;

    if( x < 0 )
	    xx = -x;
    else
	    xx = x;

    if( x <= 2.0 ) {
	    z = xx * xx;
	    if( x < 1.0e-3 )
		    return( 1.0 - 0.25*z );

	    p = (z-DR1) * polevl( z, JPJ0, 4);
	    return( p );
	}

    q = 1.0/x;
    w = sqrt(q);

    p = w * polevl( q, MOJ0, 7);
    w = q*q;
    xn = q * polevl( w, PHJ0, 7) - PIO4F;
    p = p * cos(xn + xx);
    return(p);
#endif

}

