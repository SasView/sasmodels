/**
* Spherical Bessel function 3*j1(x)/x
*
* Used for low q to avoid cancellation error.
* Note that the values differ from sasview ~ 5e-12 rather than 5e-14, but
* in this case it is likely cancellation errors in the original expression
* using double precision that are the source.  Single precision only
* requires the first 3 terms.  Double precision requires the 4th term.
* The fifth term is not needed, and is commented out.
* Taylor expansion:
*      1.0 + q2*(-3./30. + q2*(3./840.))+ q2*(-3./45360. + q2*(3./3991680.))))
*/

double sph_j1c(double q);
double sph_j1c(double q)
{
    const double q2 = q*q;
    double sin_q, cos_q;

    SINCOS(q, sin_q, cos_q);

#if FLOAT_SIZE>4
// For double precision, choose a cutoff of 0.18, which is the lower limit
// for the trig expression to 14 digits.  If only 12 digits are needed, then
// only 4 terms of the Taylor expansion are needed.
#define CUTOFF 0.18
#else
// For single precision, choose a cutoff halfway between the single precision
// lower limit for the trig expression of 1.03 and the upper limit of 1.3
// for the Taylor series expansion with 5 terms (or 9 if you count the zeros).
// For 5 digits of precision, you can drop two terms of the taylor series
// and choose a cutoff of 0.3.
#define CUTOFF 1.15
#endif

    const double bessel = (q < CUTOFF)
        ? (1.0 + q2*(-3./30. + q2*(3./840. + q2*(-3./45360. + q2*(3./3991680.)))))
        : 3.0*(sin_q/q - cos_q)/q2;
    return bessel;

/*
    // Code to test various expressions
    // tested using sympy.mpmath.mp with
    //    mp.dps = 50 (which is good for q down to 1e-6)
    //    def j1c(x): return 3*(mp.sin(x)/x - mp.cos(x))/x**2
    double ret;
    if (sizeof(q2) > 8) {

        ret = 3.0*(sinl(q)/q - cosl(q))/q2;
printf("%.15Lf %.15Lg\n", q, ret);
    //} else if (q < 0.384038453352533) {
    //} else if (q < 0.18) {
    } else if (q < 18e-10) {
        // NB: all are good to 5 digits below 0.18f
        // good is defined as 14 digits for double and 7 for single
        //ret = 1.0 + q2*q2*(3./840.) - q2*(3./30.); // good below 0.02d 0.34f
        ret = square((1. + 3./5600.*q2*q2) - q2/20.); // good below 0.03d 0.08f
        //ret = 1.0 - q2*(3./30.); // good below 0.02d 0.07f
        //ret = 1.0 + q2*(-3./30. + q2*(3./840.)); // good below 0.02d 0.34f
        //ret = 1.0 + q2*(-3./30. + q2*(3./840. + q2*(-3./45360.))); // good below 0.1d 0.8f, 12 digits below 0.18d
        //ret = 1.0 + q2*(-3./30. + q2*(3./840. + q2*(-3./45360. + q2*(3./3991680.)))); // good below 0.18d 1.3f
printf("%.15g %.15g\n", q, ret);
    } else {
        // NB: can use a cutoff of 0.1f if the goal is 5 digits rather than 7
        ret = 3.0*(sin_q/q - cos_q)/q2; // good above 0.18d, 1.03f
printf("%.15g %.15g\n", q, ret);
    }
    return ret;
*/
}
