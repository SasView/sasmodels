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
* Expression returned from Herbie (herbie.uwpise.org/demo):
*      const double t = ((1. + 3.*q2*q2/5600.) - q2/20.);
*      return t*t;
*/

double sph_j1c(double q);
double sph_j1c(double q)
{
    const double q2 = q*q;
    double sin_q, cos_q;

    SINCOS(q, sin_q, cos_q);

    const double bessel = (q < 0.384038453352533)
        ? (1.0 + q2*(-3./30. + q2*(3./840.))
        : 3.0*(sin_q/q - cos_q)/q2;

    return bessel;

 /*
    // Code to test various expressions
    if (sizeof(q2) > 4) {
        return 3.0*(sin_q/q - cos_q)/q2;
    } else if (q < 0.384038453352533) {
        //const double t = ((1. + 3.*q2*q2/5600.) - q2/20.); return t*t;
        return 1.0 + q2*q2*(3./840.) - q2*(3./30.);
        //return 1.0 + q2*(-3./30. + q2*(3./840.));
        //return 1.0 + q2*(-3./30. + q2*(3./840. + q2*(-3./45360.)));
        //return 1.0 + q2*(-3./30. + q2*(3./840. + q2*(-3./45360. + q2*(3./3991680.))));
    } else {
        return 3.0*(sin_q/q - cos_q)/q2;
    }
*/
}
