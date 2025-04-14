static double
form_volume(double radius, double length)
{
    return M_PI*radius*radius*length;
}

static double
_fq(double qab, double qc, double radius, double length)
{
    return sas_2J1x_x(qab*radius) * sas_sinx_x(qc*0.5*length);
}

static double
radius_from_excluded_volume(double radius, double length)
{
    return 0.5*cbrt(0.75*radius*(2.0*radius*length
           + (radius + length)*(M_PI*radius + length)));
}

static double
radius_from_volume(double radius, double length)
{
    return cbrt(M_PI*radius*radius*length/M_4PI_3);
}

static double
radius_from_diagonal(double radius, double length)
{
    return sqrt(radius*radius + 0.25*length*length);
}

static double
radius_effective(int mode, double radius, double length)
{
    switch (mode) {
    default:
    case 1:
        return radius_from_excluded_volume(radius, length);
    case 2:
        return radius_from_volume(radius, length);
    case 3:
        return radius;
    case 4:
        return 0.5*length;
    case 5:
        return (radius < 0.5*length ? radius : 0.5*length);
    case 6:
        return (radius > 0.5*length ? radius : 0.5*length);
    case 7:
        return radius_from_diagonal(radius,length);
    }
}

void integrand_cylinder(
double x,
double q,
double radius,
double length,
double *F1, double *F2, int n, int i){
    double sin_theta, cos_theta;
    get_sin_x(n, i, &sin_theta);
    get_cos_x(n, i, &cos_theta);
    const double form = _fq(q*sin_theta, q*cos_theta, radius, length);
    *F1 = form * sin_theta;
    *F2 = form * form * sin_theta;
}

void integrate_cylinder(
double a,
double b,
double q,
double radius,
double length,
double* res1, double* res2){

    double A = q*length/2;
    double B = q*radius;
    // Determine the number of points for the Gauss quadrature
    int expo = (int)(eval_poly(log2m(max(limits[0][0],min(limits[0][1], A))), log2m(max(limits[1][0],min(limits[1][1], B)))) + 1);
    int n = (int)(pow(2, max(1, min(15, expo))));

    double *xg, *wg;
    get_gauss_points(n, &xg, &wg);

    // Perform the integration
    *res1 = 0;
    *res2 = 0;

    for (int i = 0; i < n; i++){
        double t1, t2;
        integrand_cylinder(a + (b - a) * 0.5 * (xg[i] + 1), q, radius, length,  &t1, &t2, n, i);
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
    double solvent_sld,
    double radius,
    double length)
{
    double total_F1, total_F2;
    integrate_cylinder(0, M_PI_2, q, radius, length, &total_F1, &total_F2);

    const double s = (sld - solvent_sld) * form_volume(radius, length);
    *F1 = 1e-2 * s * total_F1;
    *F2 = 1e-4 * s * s * total_F2;
}



static double
Iqac(double qab, double qc,
    double sld,
    double solvent_sld,
    double radius,
    double length)
{
    const double form = _fq(qab, qc, radius, length);
    const double s = (sld-solvent_sld) * form_volume(radius, length);
    return 1.0e-4 * square(s * form);
}

