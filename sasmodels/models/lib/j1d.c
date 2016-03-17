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
/*							y1.c
 *
 *	Bessel function of second kind of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y1();
 *
 * y = y1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind of order one
 * of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 25 term Chebyshev
 * expansion is used, and a call to j1() is required.
 * In the second, the asymptotic trigonometric representation
 * is employed using two rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       0, 30       10000       8.6e-17     1.3e-17
 *    IEEE      0, 30       30000       1.0e-15     1.3e-16
 *
 * (error criterion relative when |y1| > 1).
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/
double j1(double );

double j1(double x) {

//Cephes double pression function
#if FLOAT_SIZE>4

    double w, z, p, q, xn;

    const double DR1 = 5.78318596294678452118E0;
    const double DR2 = 3.04712623436620863991E1;
    const double Z1 = 1.46819706421238932572E1;
    const double Z2 = 4.92184563216946036703E1;
    const double THPIO4 =  2.35619449019234492885;
    const double SQ2OPI = 0.79788456080286535588;

    double RP[8] = {
    -8.99971225705559398224E8,
    4.52228297998194034323E11,
    -7.27494245221818276015E13,
    3.68295732863852883286E15,
    0.0,
    0.0,
    0.0,
    0.0
    };

    double RQ[8] = {
    /* 1.00000000000000000000E0,*/
    6.20836478118054335476E2,
    2.56987256757748830383E5,
    8.35146791431949253037E7,
    2.21511595479792499675E10,
    4.74914122079991414898E12,
    7.84369607876235854894E14,
    8.95222336184627338078E16,
    5.32278620332680085395E18,
    };

    double PP[8] = {
    7.62125616208173112003E-4,
    7.31397056940917570436E-2,
    1.12719608129684925192E0,
    5.11207951146807644818E0,
    8.42404590141772420927E0,
    5.21451598682361504063E0,
    1.00000000000000000254E0,
    0.0,
    };
    double PQ[8] = {
    5.71323128072548699714E-4,
    6.88455908754495404082E-2,
    1.10514232634061696926E0,
    5.07386386128601488557E0,
    8.39985554327604159757E0,
    5.20982848682361821619E0,
    9.99999999999999997461E-1,
    0.0,
    };

    double QP[8] = {
    5.10862594750176621635E-2,
    4.98213872951233449420E0,
    7.58238284132545283818E1,
    3.66779609360150777800E2,
    7.10856304998926107277E2,
    5.97489612400613639965E2,
    2.11688757100572135698E2,
    2.52070205858023719784E1,
    };

    double QQ[8] = {
    /* 1.00000000000000000000E0,*/
    7.42373277035675149943E1,
    1.05644886038262816351E3,
    4.98641058337653607651E3,
    9.56231892404756170795E3,
    7.99704160447350683650E3,
    2.82619278517639096600E3,
    3.36093607810698293419E2,
    0.0,
    };

    w = x;
    if( x < 0 )
	    w = -x;

    if( w <= 5.0 )
	{
	    z = x * x;
	    w = polevl( z, RP, 3 ) / p1evl( z, RQ, 8 );
	    w = w * x * (z - Z1) * (z - Z2);
	    return( w );
	}

    w = 5.0/x;
    z = w * w;
    p = polevl( z, PP, 6)/polevl( z, PQ, 6 );
    q = polevl( z, QP, 7)/p1evl( z, QQ, 7 );
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


    double JP[8] = {
        -4.878788132172128E-009,
        6.009061827883699E-007,
        -4.541343896997497E-005,
        1.937383947804541E-003,
        -3.405537384615824E-002,
        0.0,
        0.0,
        0.0
    };

    double MO1[8] = {
        6.913942741265801E-002,
        -2.284801500053359E-001,
        3.138238455499697E-001,
        -2.102302420403875E-001,
        5.435364690523026E-003,
        1.493389585089498E-001,
        4.976029650847191E-006,
        7.978845453073848E-001
    };

    double PH1[8] = {
        -4.497014141919556E+001,
        5.073465654089319E+001,
        -2.485774108720340E+001,
        7.222973196770240E+000,
        -1.544842782180211E+000,
        3.503787691653334E-001,
        -1.637986776941202E-001,
        3.749989509080821E-001
    };

    xx = x;
    if( xx < 0 )
	    xx = -x;

    if( xx <= 2.0 )
	{
	    z = xx * xx;
	    p = (z-Z1) * xx * polevl( z, JP, 4 );
	    return( p );
	}

    q = 1.0/x;
    w = sqrt(q);

    p = w * polevl( q, MO1, 7);
    w = q*q;
    xn = q * polevl( w, PH1, 7) - THPIO4F;
    p = p * cos(xn + xx);

    return(p);
#endif
}

