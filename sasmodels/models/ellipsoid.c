static double
form_volume(double radius_polar, double radius_equatorial)
{
    return M_4PI_3*radius_polar*radius_equatorial*radius_equatorial;
}

static  double
Iq(double q,
    double sld,
    double sld_solvent,
    double radius_polar,
    double radius_equatorial)
{
    // Using ratio v = Rp/Re, we can implement the form given in Guinier (1955)
    //     i(h) = int_0^pi/2 Phi^2(h a sqrt(cos^2 + v^2 sin^2) cos dT
    //          = int_0^pi/2 Phi^2(h a sqrt((1-sin^2) + v^2 sin^2) cos dT
    //          = int_0^pi/2 Phi^2(h a sqrt(1 + sin^2(v^2-1)) cos dT
    // u-substitution of
    //     u = sin, du = cos dT
    //     i(h) = int_0^1 Phi^2(h a sqrt(1 + u^2(v^2-1)) du
    const double v_square_minus_one = square(radius_polar/radius_equatorial) - 1.0;

    // translate a point in [-1,1] to a point in [0, 1]
    // const double u = GAUSS_Z[i]*(upper-lower)/2 + (upper+lower)/2;
    const double zm = 0.5;
    const double zb = 0.5;
    double total = 0.0;
    for (int i=0;i<GAUSS_N;i++) {
        const double u = GAUSS_Z[i]*zm + zb;
        const double r = radius_equatorial*sqrt(1.0 + u*u*v_square_minus_one);
        const double f = sas_3j1x_x(q*r);
        total += GAUSS_W[i] * f * f;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    const double form = total*zm;
    const double s = (sld - sld_solvent) * form_volume(radius_polar, radius_equatorial);
    return 1.0e-4 * s * s * form;
}

static double
Iqac(double qab, double qc,
    double sld,
    double sld_solvent,
    double radius_polar,
    double radius_equatorial)
{
    const double qr = sqrt(square(radius_equatorial*qab) + square(radius_polar*qc));
    const double f = sas_3j1x_x(qr);
    const double s = (sld - sld_solvent) * form_volume(radius_polar, radius_equatorial);

    return 1.0e-4 * square(f * s);
}
