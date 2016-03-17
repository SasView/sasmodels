/*							lgamma.c
 *
 *	Log Gamma function
 *
 */
/*							lgamma()
 *
 *	Natural logarithm of gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, lgamma();
 * extern int sgngam;
 *
 * y = lgamma( x );
 *
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
 * Arguments greater than MAXLGM return MAXNUM and an error
 * message.  MAXLGM = 2.035093e36 for DEC
 * arithmetic or 2.556348e305 for IEEE arithmetic.
 *
 *
 *
 * ACCURACY:
 *
 *
 * arithmetic      domain        # trials     peak         rms
 *    DEC     0, 3                  7000     5.2e-17     1.3e-17
 *    DEC     2.718, 2.035e36       5000     3.9e-17     9.9e-18
 *    IEEE    0, 3                 28000     5.4e-16     1.1e-16
 *    IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
 * The error criterion was relative when the function magnitude
 * was greater than one but absolute when it was less than one.
 *
 * The following test used the relative error criterion, though
 * at certain points the relative error could be much higher than
 * indicated.
 *    IEEE    -200, -4             10000     4.8e-16     1.3e-16
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
*/


double lgamma( double );

double lgamma( double x) {

#if FLOAT_SIZE > 4
    double p, q, u, w, z;
    int i;
    int sgngam = 1;

    const double LS2PI  =  0.91893853320467274178;
    const double MAXLGM = 2.556348e305;
    const double MAXNUM =  1.79769313486231570815E308;
    const double LOGPI = 1.14472988584940017414;
    const double PI = M_PI;

    double A[8] = {
        8.11614167470508450300E-4,
        -5.95061904284301438324E-4,
        7.93650340457716943945E-4,
        -2.77777777730099687205E-3,
        8.33333333333331927722E-2,
        0.0,
        0.0,
        0.0
    };
    double B[8] = {
        -1.37825152569120859100E3,
        -3.88016315134637840924E4,
        -3.31612992738871184744E5,
        -1.16237097492762307383E6,
        -1.72173700820839662146E6,
        -8.53555664245765465627E5,
        0.0,
        0.0
    };

    double C[8] = {
        /* 1.00000000000000000000E0, */
        -3.51815701436523470549E2,
        -1.70642106651881159223E4,
        -2.20528590553854454839E5,
        -1.13933444367982507207E6,
        -2.53252307177582951285E6,
        -2.01889141433532773231E6,
        0.0,
        0.0
    };

    sgngam = 1;

    if( x < -34.0 ) {
	    q = -x;
	    w = lanczos_gamma(q); /* note this modifies sgngam! */
	    p = floor(q);
	    if( p == q ) {
            lgsing:
        		goto loverf;
		}
	    i = p;
	    if( (i & 1) == 0 )
		    sgngam = -1;
	    else
		    sgngam = 1;
	    z = q - p;
	    if( z > 0.5 ) {
		    p += 1.0;
		    z = p - q;
		}
	    z = q * sin( PI * z );
	    if( z == 0.0 )
		    goto lgsing;
	    z = LOGPI - log( z ) - w;
	    return( z );
	}

    if( x < 13.0 ) {
	    z = 1.0;
	    p = 0.0;
	    u = x;
	    while( u >= 3.0 ) {
		    p -= 1.0;
		    u = x + p;
		    z *= u;
		}
	    while( u < 2.0 ) {
		    if( u == 0.0 )
			    goto lgsing;
		    z /= u;
		    p += 1.0;
		    u = x + p;
		}
	    if( z < 0.0 ) {
		    sgngam = -1;
		    z = -z;
		}
	    else
		    sgngam = 1;
	    if( u == 2.0 )
		    return( log(z) );
	    p -= 2.0;
	    x = x + p;
	    p = x * polevl( x, B, 5 ) / p1evl( x, C, 6);
	    return( log(z) + p );
	}

    if( x > MAXLGM ) {
        loverf:
	        return( sgngam * MAXNUM );
	}

    q = ( x - 0.5 ) * log(x) - x + LS2PI;
    if( x > 1.0e8 )
	    return( q );

    p = 1.0/(x*x);
    if( x >= 1000.0 )
	    q += ((   7.9365079365079365079365e-4 * p
		    - 2.7777777777777777777778e-3) *p
		    + 0.0833333333333333333333) / x;
    else
	    q += polevl( p, A, 4 ) / x;
    return( q );
#else

    double p, q, w, z, xx;
    double nx, tx;
    int i, direction;
    int sgngamf = 1;

     double B[8] = {
        6.055172732649237E-004,
        -1.311620815545743E-003,
        2.863437556468661E-003,
        -7.366775108654962E-003,
        2.058355474821512E-002,
        -6.735323259371034E-002,
        3.224669577325661E-001,
        4.227843421859038E-001
    };

    /* log gamma(x+1), -.25 < x < .25 */
    double C[8] = {
        1.369488127325832E-001,
        -1.590086327657347E-001,
        1.692415923504637E-001,
        -2.067882815621965E-001,
        2.705806208275915E-001,
        -4.006931650563372E-001,
        8.224670749082976E-001,
        -5.772156501719101E-001
    };

    /* log( sqrt( 2*pi ) ) */
    const double LS2PI  =  0.91893853320467274178;
    const double MAXLGM = 2.035093e36;
    const double PIINV =  0.318309886183790671538;
    const double MAXNUMF = 3.4028234663852885981170418348451692544e38;
    const double PIF = 3.141592653589793238;
    /* Logarithm of gamma function */

    sgngamf = 1;

    xx = x;
    if( xx < 0.0 ) {
	    q = -xx;
	    w = lgamma(q); /* note this modifies sgngam! */
	    p = floor(q);
	    if( p == q )
		    goto loverf;
	    i = p;

	    if( (i & 1) == 0 )
		    sgngamf = -1;
	    else
		    sgngamf = 1;

	    z = q - p;

	    if( z > 0.5 ) {
		    p += 1.0;
		    z = p - q;
	    }
	    z = q * sin( PIF * z );
	    if( z == 0.0 )
		    goto loverf;
	    z = -log( PIINV*z ) - w;
	    return( z );
	}

    if( x < 6.5 ) {
	    direction = 0;
	    z = 1.0;
	    tx = x;
	    nx = 0.0;
	    if( x >= 1.5 ) {
		    while( tx > 2.5 ) {
			    nx -= 1.0;
			    tx = x + nx;
			    z *=tx;
			}
		    x += nx - 2.0;
            iv1r5:
		    p = x * polevl( x, B, 7 );
		    goto cont;
	}
	if( x >= 1.25 ) {
		z *= x;
		x -= 1.0; /* x + 1 - 2 */
		direction = 1;
		goto iv1r5;
	}
	if( x >= 0.75 ) {
		x -= 1.0;
		p = x * polevl( x, C, 7 );
		q = 0.0;
		goto contz;
	}
	while( tx < 1.5 ) {
		if( tx == 0.0 )
			goto loverf;
		z *=tx;
		nx += 1.0;
		tx = x + nx;
	}
	direction = 1;
	x += nx - 2.0;
	p = x * polevl( x, B, 7 );

    cont:
	    if( z < 0.0 ) {
		    sgngamf = -1;
		    z = -z;
        }
	    else {
		    sgngamf = 1;
		}
	    q = log(z);
	    if( direction )
		    q = -q;
    contz:
	    return( p + q );
	}

    if( x > MAXLGM ) {
        loverf:
	    return( sgngamf * MAXNUMF );
	}

    /* Note, though an asymptotic formula could be used for x >= 3,
    * there is cancellation error in the following if x < 6.5.  */
    q = LS2PI - x;
    q += ( x - 0.5 ) * log(x);

    if( x <= 1.0e4 ) {
	    z = 1.0/x;
	    p = z * z;
	    q += ((    6.789774945028216E-004 * p
		 - 2.769887652139868E-003 ) * p
		+  8.333316229807355E-002 ) * z;
	}
    return( q );
#endif
}