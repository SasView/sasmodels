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

void integrandF1F2(
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
    integrate2(integrandF1F2, 0, M_PI_2, q, radius, length, &total_F1, &total_F2);

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

