static double
form_volume(double radius_polar, double radius_equatorial)
{
    return M_4PI_3*radius_polar*radius_equatorial*radius_equatorial;
}

static double
radius_from_volume(double radius_polar, double radius_equatorial)
{
    return cbrt(radius_polar*radius_equatorial*radius_equatorial);
}

static double
radius_from_curvature(double radius_polar, double radius_equatorial)
{
    // Trivial cases
    if (radius_polar == radius_equatorial) return radius_polar;
    if (radius_polar * radius_equatorial == 0.)  return 0.;

    // see equation (26) in A.Isihara, J.Chem.Phys. 18(1950)1446-1449
    const double ratio = (radius_polar < radius_equatorial
                          ? radius_polar / radius_equatorial
                          : radius_equatorial / radius_polar);
    const double e1 = sqrt(1.0 - ratio*ratio);
    const double b1 = 1.0 + asin(e1) / (e1 * ratio);
    const double bL = (1.0 + e1) / (1.0 - e1);
    const double b2 = 1.0 + 0.5 * ratio * ratio / e1 * log(bL);
    const double delta = 0.75 * b1 * b2;
    const double ddd = 2.0 * (delta + 1.0) * radius_polar * radius_equatorial * radius_equatorial;
    return 0.5 * cbrt(ddd);
}

static double
radius_effective(int mode, double radius_polar, double radius_equatorial)
{
    switch (mode) {
    default:
    case 1: // average curvature
        return radius_from_curvature(radius_polar, radius_equatorial);
    case 2: // equivalent volume sphere
        return radius_from_volume(radius_polar, radius_equatorial);
    case 3: // min radius
        return (radius_polar < radius_equatorial ? radius_polar : radius_equatorial);
    case 4: // max radius
        return (radius_polar > radius_equatorial ? radius_polar : radius_equatorial);
    }
}



static void
integrand_ellipsoid(
    const double u,
    const double q,
    double radius_polar,
    double radius_equatorial,
    double* res1,
    double* res2){

    // Using ratio v = Rp/Re, we can implement the form given in Guinier (1955)
    //     i(h) = int_0^pi/2 Phi^2(h a sqrt(cos^2 + v^2 sin^2) cos dT
    //          = int_0^pi/2 Phi^2(h a sqrt((1-sin^2) + v^2 sin^2) cos dT
    //          = int_0^pi/2 Phi^2(h a sqrt(1 + sin^2(v^2-1)) cos dT
    // u-substitution of
    //     u = sin, du = cos dT
    //     i(h) = int_0^1 Phi^2(h a sqrt(1 + u^2(v^2-1)) du

    const double v_square_minus_one = square(radius_polar/radius_equatorial) - 1.0;
    const double r = radius_equatorial*sqrt(1.0 + u*u*v_square_minus_one);
    const double f = sas_3j1x_x(q*r);

    *res1 = f;
    *res2 = f * f;
}

void integrate_ellipsoid(
    const double a,
    const double b,
    const double q,
    double radius_polar,
    double radius_equatorial,
    double* res1,
    double* res2
    ){


    const double A = radius_equatorial;
    const double B = radius_polar;

    int expo = (int)(eval_poly(log2m(max(limits[0][0],min(limits[0][1], q))), log2m(max(limits[1][0],min(limits[1][1], A))), log2m(max(limits[2][0],min(limits[2][1], B)))) + 1);
    int n = (int)(pow(2, max(1, min(15, expo))));

    double *xg, *wg;
    get_gauss_points(n, &xg, &wg);

    // Perform the integration
    *res1 = 0;
    *res2 = 0;

    for (int i = 0; i < n; i++){
        double t1, t2;
        integrand_ellipsoid(a + (b - a) * 0.5 * (xg[i] + 1), q, radius_polar, radius_equatorial, &t1, &t2);
        *res1 += t1 * wg[i];
        *res2 += t2 * wg[i];
    }

    *res1 *= (b - a) * 0.5;
    *res2 *= (b - a) * 0.5;


}


static void
Fq(double q,
    double *F1,
    double *F2,
    double sld,
    double sld_solvent,
    double radius_polar,
    double radius_equatorial)
{


    double total_F1, total_F2;
    integrate_ellipsoid(0, 1, q, radius_polar, radius_equatorial, &total_F1, &total_F2);

    const double s = (sld - sld_solvent) * form_volume(radius_polar, radius_equatorial);
    *F1 = 1e-2 * s * total_F1;
    *F2 = 1e-4 * s * s * total_F2;
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
