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
radius_from_volume(double radius, double thickness, double length)
{
    const double volume_outer_cyl = form_volume(radius,thickness,length);
    return cbrt(0.75*volume_outer_cyl/M_PI);
}

static double
radius_from_diagonal(double radius, double thickness, double length)
{
    const double radius_outer = radius + thickness;
    const double length_outer = length + thickness;
    return sqrt(radius_outer*radius_outer + 0.25*length_outer*length_outer);
}

static double
effective_radius(int mode, double radius, double thickness, double length)
{
    if (mode == 1) {
        return radius_from_volume(radius, thickness, length);
    } else if (mode == 2) {
        return radius + thickness;
    } else if (mode == 3) {
        return 0.5*length + thickness;
    } else if (mode == 4) {
        return (radius < 0.5*length ? radius + thickness : 0.5*length + thickness);
    } else if (mode == 5) {
        return (radius > 0.5*length ? radius + thickness : 0.5*length + thickness);
    } else {
        return radius_from_diagonal(radius,thickness,length);
    }
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
    double total_F1 = 0.0;
    double total_F2 = 0.0;
    for (int i=0; i<GAUSS_N ;i++) {
        // translate a point in [-1,1] to a point in [0, pi/2]
        //const double theta = ( GAUSS_Z[i]*(upper-lower) + upper + lower )/2.0;
        double sin_theta, cos_theta;
        const double theta = GAUSS_Z[i]*M_PI_4 + M_PI_4;
        SINCOS(theta, sin_theta,  cos_theta);
        const double qab = q*sin_theta;
        const double qc = q*cos_theta;
        const double fq = _cyl(core_vd, core_r*qab, core_h*qc)
            + _cyl(shell_vd, shell_r*qab, shell_h*qc);
        total_F1 += GAUSS_W[i] * fq * sin_theta;
        total_F2 += GAUSS_W[i] * fq * fq * sin_theta;
    }
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
