// Use taylor series for low q to avoid cancellation error.
// Note that the values differ from sasview ~ 5e-12 rather than 5e-14, but
// in this case it is likely cancellation errors in the original expression
// using double precision that are the source.  Single precision only
// requires the first 3 terms.  Double precision requires the 4th term.
// The fifth term is not needed, and is commented out.
double J1c(double qr);
double J1c(double qr)
{
    const double qr2 = qr*qr;
    double sn, cn;

    SINCOS(qr, sn, cn);

    const double bes = (qr < 1.e-1) ?
        1.0 + qr2*(-3./30. + qr2*(3./840. + qr2*(-3./45360.)))
        // + qr2*(3./3991680.)))
        : 3.0*(sn/qr - cn)/qr2;

    return bes;
}
