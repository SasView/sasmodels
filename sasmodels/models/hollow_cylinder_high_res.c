static double
shell_volume(double radius, double thickness, double length)
{
    return M_PI*length*(square(radius+thickness) - radius*radius);
}

static double
form_volume(double radius, double thickness, double length)
{
    return M_PI*length*square(radius+thickness);
}

static double
radius_from_excluded_volume(double radius, double thickness, double length)
{
    const double radius_tot = radius + thickness;
    return 0.5*cbrt(0.75*radius_tot*(2.0*radius_tot*length + (radius_tot + length)*(M_PI*radius_tot + length)));
}

static double
radius_from_volume(double radius, double thickness, double length)
{
    const double volume_outer_cyl = M_PI*square(radius + thickness)*length;
    return cbrt(volume_outer_cyl/M_4PI_3);
}

static double
radius_from_diagonal(double radius, double thickness, double length)
{
    return sqrt(square(radius + thickness) + 0.25*square(length));
}

static double
radius_effective(int mode, double radius, double thickness, double length)
{
    switch (mode) {
    default:
    case 1: // excluded volume
        return radius_from_excluded_volume(radius, thickness, length);
    case 2: // equivalent volume sphere
        return radius_from_volume(radius, thickness, length);
    case 3: // outer radius
        return radius + thickness;
    case 4: // half length
        return 0.5*length;
    case 5: // half outer min dimension
        return (radius + thickness < 0.5*length ? radius + thickness : 0.5*length);
    case 6: // half outer max dimension
        return (radius + thickness > 0.5*length ? radius + thickness : 0.5*length);
    case 7: // half outer diagonal
        return radius_from_diagonal(radius,thickness,length);
    }
}

static double
_fq(double qab, double qc,
    double radius, double thickness, double length)
{
    const double lam1 = sas_2J1x_x((radius+thickness)*qab);
    const double lam2 = sas_2J1x_x(radius*qab);
    const double gamma_sq = square(radius/(radius+thickness));
    //Note: lim_{thickness -> 0} psi = sas_J0(radius*qab)
    //Note: lim_{radius -> 0} psi = sas_2J1x_x(thickness*qab)
    const double psi = (lam1 - gamma_sq*lam2)/(1.0 - gamma_sq);    //SRK 10/19/00
    const double t2 = sas_sinx_x(0.5*length*qc);
    return psi*t2;
}

void integrandF1F2(
double x,
double q,
double radius,
double thickness,
double length,
double *F1, double *F2, int n, int i){
    double sin_theta, cos_theta;
    cos_theta = x;
    sin_theta = sqrt(1.0 - cos_theta*cos_theta);
    const double form = _fq(q*sin_theta, q*cos_theta, radius, thickness, length);
    *F1 = form;
    *F2 = form * form;
}

void integrate2(
double a,
double b,
double q,
double radius,
double thickness,
double length,
double* res1, double* res2){

    // Determine the number of points for the Gauss quadrature
    int expo = (int)(eval_poly(log2m(max(limits[0][0],min(limits[0][1], q))), log2m(max(limits[1][0],min(limits[1][1], radius))), log2m(max(limits[2][0],min(limits[2][1], thickness))), log2m(max(limits[3][0],min(limits[3][1], length)))) + 1);
    int n = (int)(pow(2, max(1, min(15, expo))));

    double *xg, *wg;
    get_gauss_points(n, &xg, &wg);

    // Perform the integration
    *res1 = 0;
    *res2 = 0;

    for (int i = 0; i < n; i++){
        double t1, t2;
        integrandF1F2(a + (b - a) * 0.5 * (xg[i] + 1), q, radius, thickness, length,  &t1, &t2, n, i);
        *res1 += t1 * wg[i];
        *res2 += t2 * wg[i];
    }

    *res1 *= (b - a) * 0.5;
    *res2 *= (b - a) * 0.5;
}

static void
Fq(double q, double *F1, double *F2, double radius, double thickness, double length,
    double sld, double solvent_sld)
{
    const double lower = 0.0;
    const double upper = 1.0;        //limits of numerical integral

    double total_F1, total_F2;
    integrate2(lower, upper, q, radius, thickness, length, &total_F1, &total_F2);
    const double s = (sld - solvent_sld) * shell_volume(radius, thickness, length);
    *F1 = 1e-2 * s * total_F1;
    *F2 = 1e-4 * s*s * total_F2;
}


static double
Iqac(double qab, double qc,
    double radius, double thickness, double length,
    double sld, double solvent_sld)
{
    const double form = _fq(qab, qc, radius, thickness, length);
    const double s = (sld - solvent_sld) * shell_volume(radius, thickness, length);
    return 1.0e-4*square(s * form);
}
