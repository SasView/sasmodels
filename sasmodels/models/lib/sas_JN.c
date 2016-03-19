/*							jn.c
 *
 *	Bessel function of integer order
 *
 *
 *
 * SYNOPSIS:
 *
 * int n;
 * double x, y, jn();
 *
 * y = jn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order n, where n is a
 * (possibly negative) integer.
 *
 * The ratio of jn(x) to j0(x) is computed by backward
 * recurrence.  First the ratio jn/jn-1 is found by a
 * continued fraction expansion.  Then the recurrence
 * relating successive orders is applied until j0 or j1 is
 * reached.
 *
 * If n = 0 or 1 the routine for j0 or j1 is called
 * directly.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   range      # trials      peak         rms
 *    DEC       0, 30        5500       6.9e-17     9.3e-18
 *    IEEE      0, 30        5000       4.4e-16     7.9e-17
 *
 *
 * Not suitable for large n or x. Use jv() instead.
 *
 */

/*							jn.c
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

double sas_JN( int n, double x );

double sas_JN( int n, double x ) {

    const double MACHEP = 1.11022302462515654042E-16;
    double pkm2, pkm1, pk, xk, r, ans, xinv;
    int k, sign;

    if( n < 0 ) {
	    n = -n;
	    if( (n & 1) == 0 )	/* -1**n */
		    sign = 1;
	    else
		    sign = -1;
	}
    else
	    sign = 1;

    if( x < 0.0 ) {
	    if( n & 1 )
		    sign = -sign;
	    x = -x;
	}

    if( n == 0 )
	    return( sign * sas_J0(x) );
    if( n == 1 )
	    return( sign * sas_J1(x) );
    if( n == 2 )
	    return( sign * (2.0 * sas_J1(x) / x  -  sas_J0(x)) );

    if( x < MACHEP )
	    return( 0.0 );

    #if FLOAT_SIZE > 4
        k = 53;
    #else
        k = 24;
    #endif

    pk = 2 * (n + k);
    ans = pk;
    xk = x * x;

    do {
	    pk -= 2.0;
	    ans = pk - (xk/ans);
	} while( --k > 0 );

    /* backward recurrence */

    pk = 1.0;

    #if FLOAT_SIZE > 4
        ans = x/ans;
        pkm1 = 1.0/ans;

        k = n-1;
        r = 2 * k;

        do {
            pkm2 = (pkm1 * r  -  pk * x) / x;
	        pk = pkm1;
	        pkm1 = pkm2;
	        r -= 2.0;
	    } while( --k > 0 );

        if( fabs(pk) > fabs(pkm1) )
	        ans = sas_J1(x)/pk;
        else
	        ans = sas_J0(x)/pkm1;

	    return( sign * ans );

    #else
        xinv = 1.0/x;
        pkm1 = ans * xinv;
        k = n-1;
        r = (float )(2 * k);

        do {
	        pkm2 = (pkm1 * r  -  pk * x) * xinv;
	        pk = pkm1;
	        pkm1 = pkm2;
	        r -= 2.0;
	    }
        while( --k > 0 );

        r = pk;
        if( r < 0 )
	        r = -r;
        ans = pkm1;
        if( ans < 0 )
	        ans = -ans;

        if( r > ans )  /* if( fabs(pk) > fabs(pkm1) ) */
	        ans = sign * sas_J1(x)/pk;
        else
	        ans = sign * sas_J0(x)/pkm1;
        return( ans );
    #endif
}

