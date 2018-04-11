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
Iq(double q,
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
    double total = 0.0;
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
        total += GAUSS_W[i] * fq * fq * sin_theta;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    //const double form = (upper-lower)/2.0*total;
    return 1.0e-4 * total * M_PI_4;
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
