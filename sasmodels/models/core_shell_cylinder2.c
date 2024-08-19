// vd = volume * delta_rho
// besarg = q * R * sin(theta)
// siarg = q * L/2 * cos(theta)
static double _cyl(double vd, double besarg, double siarg)
{
    return vd * sas_sinx_x(siarg) * sas_2J1x_x(besarg);
}

static double
form_volume(double radius, double thickness, double length)
{
    return M_PI*square(radius+thickness)*(length+2.0*thickness);
}

static double
radius_from_excluded_volume(double radius, double thickness, double length)
{
    const double radius_tot = radius + thickness;
    const double length_tot = length + 2.0*thickness;
    return 0.5*cbrt(0.75*radius_tot*(2.0*radius_tot*length_tot + (radius_tot + length_tot)*(M_PI*radius_tot + length_tot)));
}

static double
radius_from_volume(double radius, double thickness, double length)
{
    const double volume_outer_cyl = form_volume(radius,thickness,length);
    return cbrt(volume_outer_cyl/M_4PI_3);
}

static double
radius_from_diagonal(double radius, double thickness, double length)
{
    const double radius_outer = radius + thickness;
    const double length_outer = length + 2.0*thickness;
    return sqrt(radius_outer*radius_outer + 0.25*length_outer*length_outer);
}

static double
radius_effective(int mode, double radius, double thickness, double length)
{
    switch (mode) {
    default:
    case 1: //cylinder excluded volume
        return radius_from_excluded_volume(radius, thickness, length);
    case 2: // equivalent volume sphere
        return radius_from_volume(radius, thickness, length);
    case 3: // outer radius
        return radius + thickness;
    case 4: // half outer length
        return 0.5*length + thickness;
    case 5: // half min outer length
        return (radius < 0.5*length ? radius + thickness : 0.5*length + thickness);
    case 6: // half max outer length
        return (radius > 0.5*length ? radius + thickness : 0.5*length + thickness);
    case 7: // half outer diagonal
        return radius_from_diagonal(radius,thickness,length);
    }
}

static void
integrand(const double x, const double q, const double core_r, const double core_h, const double core_vd, const double shell_r, const double shell_h, const double shell_vd, double* res1, double* res2, int n, int i){
    double sin_theta, cos_theta;
    get_sin_x(n, i, &sin_theta);
    get_cos_x(n, i, &cos_theta);

    const double qab = q*sin_theta;
    const double qc = q*cos_theta;
    const double fq = _cyl(core_vd, core_r*qab, core_h*qc) + _cyl(shell_vd, shell_r*qab, shell_h*qc);
    *res1 = fq * sin_theta;
    *res2 = fq * fq * sin_theta;
}


static void
Fq(double q,
    double *F1,
    double *F2,
    double core_sld,
    double shell_sld,
    double solvent_sld,
    double radius,
    double thickness,
    double length)
{

    // precalculate constants
    const double core_r = radius;
    const double core_h = 0.5*length;
    const double core_vd = form_volume(radius,0,length) * (core_sld-shell_sld);
    const double shell_r = (radius + thickness);
    const double shell_h = (0.5*length + thickness);
    const double shell_vd = form_volume(radius,thickness,length) * (shell_sld-solvent_sld);
    double total_F1, total_F2;

    integrate2(integrand, 0, M_PI_2, q, core_r, core_h, core_vd, shell_r, shell_h, shell_vd, &total_F1, &total_F2);

    // translate dx in [-1,1] to dx in [lower,upper]
    //const double form = (upper-lower)/2.0*total;
    *F1 = 1.0e-2 * total_F1 * M_PI_4;
    *F2 = 1.0e-4 * total_F2 * M_PI_4;
}

static double
Iqac(double qab, double qc,
    double core_sld,
    double shell_sld,
    double solvent_sld,
    double radius,
    double thickness,
    double length)
{
    const double core_r = radius;
    const double core_h = 0.5*length;
    const double core_vd = form_volume(radius,0,length) * (core_sld-shell_sld);
    const double shell_r = (radius + thickness);
    const double shell_h = (0.5*length + thickness);
    const double shell_vd = form_volume(radius,thickness,length) * (shell_sld-solvent_sld);

    const double fq = _cyl(core_vd, core_r*qab, core_h*qc)
        + _cyl(shell_vd, shell_r*qab, shell_h*qc);
    return 1.0e-4 * fq * fq;
}
