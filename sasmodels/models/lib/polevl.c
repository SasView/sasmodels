/*							polevl.c
 *							p1evl.c
 *
 *	Evaluate polynomial
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N+1], polevl[];
 *
 * y = polevl( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates polynomial of degree N:
 *
 *                     2          N
 * y  =  C  + C x + C x  +...+ C x
 *        0    1     2          N
 *
 * Coefficients are stored in reverse order:
 *
 * coef[0] = C  , ..., coef[N] = C  .
 *            N                   0
 *
 *  The function p1evl() assumes that coef[N] = 1.0 and is
 * omitted from the array.  Its calling arguments are
 * otherwise the same as polevl().
 *
 *
 * SPEED:
 *
 * In the interest of speed, there are no checks for out
 * of bounds arithmetic.  This routine is used by most of
 * the functions in the library.  Depending on available
 * equipment features, the user may wish to rewrite the
 * program in microcode or assembly language.
 *
 */


/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

double polevlRP(double x, int N ) {

    double coef[8] = {
    -8.99971225705559398224E8,
    4.52228297998194034323E11,
    -7.27494245221818276015E13,
    3.68295732863852883286E15,
    0.0,
    0.0,
    0.0,
    0.0 };

    int i = 0;
    double ans = coef[i];

    while (i < N) {
        i++;
        ans = ans * x + coef[i];
    }
    return ans ;
}

double polevlRQ(double x, int N ) {

    double coef[8] = {
        6.20836478118054335476E2,
        2.56987256757748830383E5,
        8.35146791431949253037E7,
        2.21511595479792499675E10,
        4.74914122079991414898E12,
        7.84369607876235854894E14,
        8.95222336184627338078E16,
        5.32278620332680085395E18
    };

    int i = 0;
    double ans = coef[i];

    while (i < N) {
        i++;
        ans = ans * x + coef[i];
    }
    return ans ;
}

double polevlPP(double x, int N ) {

    double coef[8] = {
    7.62125616208173112003E-4,
    7.31397056940917570436E-2,
    1.12719608129684925192E0,
    5.11207951146807644818E0,
    8.42404590141772420927E0,
    5.21451598682361504063E0,
    1.00000000000000000254E0,
    0.0} ;

    int i = 0;
    double ans = coef[i];

    while (i < N) {
        i++;
        ans = ans * x + coef[i];
    }
    return ans ;
}

double polevlPQ(double x, int N ) {
//4: PQ
    double coef[8] = {
    5.71323128072548699714E-4,
    6.88455908754495404082E-2,
    1.10514232634061696926E0,
    5.07386386128601488557E0,
    8.39985554327604159757E0,
    5.20982848682361821619E0,
    9.99999999999999997461E-1,
    0.0 };

    int i = 0;
    double ans = coef[i];

    while (i < N) {
        i++;
        ans = ans * x + coef[i];
    }
    return ans ;
}

double polevlQP(double x, int N ) {
    double coef[8] = {
    5.10862594750176621635E-2,
    4.98213872951233449420E0,
    7.58238284132545283818E1,
    3.66779609360150777800E2,
    7.10856304998926107277E2,
    5.97489612400613639965E2,
    2.11688757100572135698E2,
    2.52070205858023719784E1 };

    int i = 0;
    double ans = coef[i];

    while (i < N) {
        i++;
        ans = ans * x + coef[i];
    }
    return ans ;

}

double polevlQQ(double x, int N ) {
    double coef[8] = {
    7.42373277035675149943E1,
    1.05644886038262816351E3,
    4.98641058337653607651E3,
    9.56231892404756170795E3,
    7.99704160447350683650E3,
    2.82619278517639096600E3,
    3.36093607810698293419E2,
    0.0 };

    int i = 0;
    double ans = coef[i];

    while (i < N) {
        i++;
        ans = ans * x + coef[i];
    }
    return ans ;

}

double polevlJP( double x, int N ) {
    double coef[8] = {
        -4.878788132172128E-009,
        6.009061827883699E-007,
        -4.541343896997497E-005,
        1.937383947804541E-003,
        -3.405537384615824E-002,
        0.0,
        0.0,
        0.0
    };

    int i = 0;
    double ans = coef[i];

    while (i < N) {
        i++;
        ans = ans * x + coef[i];
    }
    return ans ;

}

double polevlMO1( double x, int N ) {
    double coef[8] = {
        6.913942741265801E-002,
        -2.284801500053359E-001,
        3.138238455499697E-001,
        -2.102302420403875E-001,
        5.435364690523026E-003,
        1.493389585089498E-001,
        4.976029650847191E-006,
        7.978845453073848E-001
    };

    int i = 0;
    double ans = coef[i];

    while (i < N) {
        i++;
        ans = ans * x + coef[i];
    }
    return ans ;
}

double polevlPH1( double x, int N ) {

    double coef[8] = {
        -4.497014141919556E+001,
        5.073465654089319E+001,
        -2.485774108720340E+001,
        7.222973196770240E+000,
        -1.544842782180211E+000,
        3.503787691653334E-001,
        -1.637986776941202E-001,
        3.749989509080821E-001
    };

    int i = 0;
    double ans = coef[i];

    while (i < N) {
        i++;
        ans = ans * x + coef[i];
    }
    return ans ;
}

/*double polevl( double x, double coef[8], int N ) {

    int i = 0;
    double ans = coef[i];

    while (i < N) {
        i++;
        ans = ans * x + coef[i];
    }

    return ans ;

}*/

/*							p1evl()	*/
/*                                          N
 * Evaluate polynomial when coefficient of x  is 1.0.
 * Otherwise same as polevl.
 */

double p1evlRP( double x, int N ) {
    double coef[8] = {
    -8.99971225705559398224E8,
    4.52228297998194034323E11,
    -7.27494245221818276015E13,
    3.68295732863852883286E15,
    0.0,
    0.0,
    0.0,
    0.0 };

    int i=0;
    double ans = x+coef[i];

    while (i < N-1) {
        i++;
        ans = ans*x + coef[i];
    }

    return( ans );

}

double p1evlRQ( double x, int N ) {
//1: RQ
    double coef[8] = {
    6.20836478118054335476E2,
    2.56987256757748830383E5,
    8.35146791431949253037E7,
    2.21511595479792499675E10,
    4.74914122079991414898E12,
    7.84369607876235854894E14,
    8.95222336184627338078E16,
    5.32278620332680085395E18};

    int i=0;
    double ans = x+coef[i];

    while (i < N-1) {
        i++;
        ans = ans*x + coef[i];
    }

    return( ans );
}

double p1evlPP( double x, int N ) {
//3 : PP
    double coef[8] = {
    7.62125616208173112003E-4,
    7.31397056940917570436E-2,
    1.12719608129684925192E0,
    5.11207951146807644818E0,
    8.42404590141772420927E0,
    5.21451598682361504063E0,
    1.00000000000000000254E0,
    0.0};

    int i=0;
    double ans = x+coef[i];

    while (i < N-1) {
        i++;
        ans = ans*x + coef[i];
    }

    return( ans );
}

double p1evlPQ( double x, int N ) {
//4: PQ
    double coef[8] = {
    5.71323128072548699714E-4,
    6.88455908754495404082E-2,
    1.10514232634061696926E0,
    5.07386386128601488557E0,
    8.39985554327604159757E0,
    5.20982848682361821619E0,
    9.99999999999999997461E-1,
    0.0};

    int i=0;
    double ans = x+coef[i];

    while (i < N-1) {
        i++;
        ans = ans*x + coef[i];
    }

    return( ans );
}

double p1evlQP( double x, int N ) {
//5: QP
    double coef[8] = {
    5.10862594750176621635E-2,
    4.98213872951233449420E0,
    7.58238284132545283818E1,
    3.66779609360150777800E2,
    7.10856304998926107277E2,
    5.97489612400613639965E2,
    2.11688757100572135698E2,
    2.52070205858023719784E1 };

    int i=0;
    double ans = x+coef[i];

    while (i < N-1) {
        i++;
        ans = ans*x + coef[i];
    }

    return( ans );
}

double p1evlQQ( double x, int N ) {
    double coef[8] = {
    7.42373277035675149943E1,
    1.05644886038262816351E3,
    4.98641058337653607651E3,
    9.56231892404756170795E3,
    7.99704160447350683650E3,
    2.82619278517639096600E3,
    3.36093607810698293419E2,
    0.0};

    int i=0;
    double ans = x+coef[i];

    while (i < N-1) {
        i++;
        ans = ans*x + coef[i];
    }

    return( ans );

}

double p1evlJP( double x, int N ) {
    double coef[8] = {
        -4.878788132172128E-009,
        6.009061827883699E-007,
        -4.541343896997497E-005,
        1.937383947804541E-003,
        -3.405537384615824E-002,
        0.0,
        0.0,
        0.0};

    int i=0;
    double ans = x+coef[i];

    while (i < N-1) {
        i++;
        ans = ans*x + coef[i];
    }

    return( ans );
}

double p1evlMO1( double x, int N ) {
    double coef[8] = {
        6.913942741265801E-002,
        -2.284801500053359E-001,
        3.138238455499697E-001,
        -2.102302420403875E-001,
        5.435364690523026E-003,
        1.493389585089498E-001,
        4.976029650847191E-006,
        7.978845453073848E-001};

    int i=0;
    double ans = x+coef[i];

    while (i < N-1) {
        i++;
        ans = ans*x + coef[i];
    }

    return( ans );

}

double p1evlPH1( double x, int N ) {
    double coef[8] = {
        -4.497014141919556E+001,
        5.073465654089319E+001,
        -2.485774108720340E+001,
        7.222973196770240E+000,
        -1.544842782180211E+000,
        3.503787691653334E-001,
        -1.637986776941202E-001,
        3.749989509080821E-001};
    int i=0;
    double ans = x+coef[i];

    while (i < N-1) {
        i++;
        ans = ans*x + coef[i];
    }

    return( ans );
}

/*double p1evl( double x, double coef[8], int N ) {
    int i=0;
    double ans = x+coef[i];

    while (i < N-1) {
        i++;
        ans = ans*x + coef[i];
    }

    return( ans );

}*/
