static double checksum = 0.0;

void translate(
    double *z1, double *z2, double *k1, double *k2, double *volumefraction,
    double *a, double *b, double *c1, double *c2, double *d1, double *d2)
{
    int debug = 0;
// Theoretically speaking, z1 and z2 (and k1 and k2) are symetric. 
// Exchanging z1 (k1) with z2 (k2) does not altern the potential. 
// The results should be identical. However, the orginal model proposed 
// by Y. Liu treats z1 and z2 differently when implementing the numerical solution. // Hence, it is in general a good practice to require the z1 > z2. 
// The following code is added here to swap z1 (k1) with z2 (k2) if z1 < z2.

    if ( *z1 < *z2 )
    {
        double temp = *z1;
        *z1 = *z2;
        *z2 = temp;
        temp = *k1;
        *k1 = *k2;
        *k2 = temp;
    }

    else if ( *z1 == *z2) {
        // The calculator produces NaM when z1 == z2 so no need to precompute
        // from the parameter values.
        *a = NAN;
        return;
    }

    TY_SolveEquations(*z1, *z2, *k1, *k2, *volumefraction, a, b, c1, c2, d1, d2, debug );
    int checkFlag = TY_CheckSolution( *z1, *z2, *k1, *k2, *volumefraction, *a, *b, *c1, *c2, *d1, *d2 );
    if (!checkFlag) {
        *a = NAN;
        return;
    }
    //printf("Solving (z1=%g, z2=%g, k1=%g, k2=%g, Vf=%g)\n", *z1, *z2, *k1, *k2, *volumefraction);
    //printf("=> (a=%g, b=%g, c1=%g, c2=%g, d1=%g, d2=%g)\n", *a, *b, *c1, *c2, *d1, *d2);
}

// Normally constructed in generate.py, but we are doing something special here
// so create them ourselves.
// Using k1, k2 negative from the values given in the the model parameters
#define TRANSLATION_VARS(_v) \
   double z1=_v.z1, z2=_v.z2, k1=-_v.k1, k2=-_v.k2, vf=_v.volfraction; \
   double a, b, c1, c2, d1, d2; \
   translate(&z1, &z2, &k1, &k2, &vf, &a, &b, &c1, &c2, &d1, &d2)
#define OVERRIDE_IQ(_q, _v) \
   Iq(_q, _v.radius_effective, vf, k1, k2, z1, z2, a, b, c1, c2, d1, d2)

double Iq(
    double qexp,
    double radius,
    double volumefraction,
    double k1,
    double k2,
    double z1,
    double z2,
    // Hidden parameters set by translate before the q loop
    double a,
    double b,
    double c1,
    double c2,
    double d1,
    double d2
    )
{
    if (isnan(a)) {
          return NAN;
    } else {
        return SqTwoYukawa( qexp * 2 * radius, z1, z2, k1, k2, volumefraction, a, b, c1, c2, d1, d2 );
    }
}
