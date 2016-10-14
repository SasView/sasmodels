double form_volume(double radius, double thickness, double length);
double Iq(double q, double core_sld, double shell_sld, double solvent_sld,
    double radius, double thickness, double length);
double Iqxy(double qx, double qy, double core_sld, double shell_sld, double solvent_sld,
    double radius, double thickness, double length, double theta, double phi);

// twovd = 2 * volume * delta_rho
// besarg = q * R * sin(alpha)
// siarg = q * L/2 * cos(alpha)
double _cyl(double twovd, double besarg, double siarg);
double _cyl(double twovd, double besarg, double siarg)
{
    return twovd * sinc(siarg) * sas_J1c(besarg);
}

double form_volume(double radius, double thickness, double length)
{
    return M_PI*(radius+thickness)*(radius+thickness)*(length+2.0*thickness);
}

double Iq(double q,
    double core_sld,
    double shell_sld,
    double solvent_sld,
    double radius,
    double thickness,
    double length)
{
    // precalculate constants
    const double core_qr = q*radius;
    const double core_qh = q*0.5*length;
    const double core_twovd = 2.0 * form_volume(radius,0,length)
                            * (core_sld-shell_sld);
    const double shell_qr = q*(radius + thickness);
    const double shell_qh = q*(0.5*length + thickness);
    const double shell_twovd = 2.0 * form_volume(radius,thickness,length)
                             * (shell_sld-solvent_sld);
    double total = 0.0;
    // double lower=0, upper=M_PI_2;
    for (int i=0; i<76 ;i++) {
        // translate a point in [-1,1] to a point in [lower,upper]
        //const double alpha = ( Gauss76Z[i]*(upper-lower) + upper + lower )/2.0;
        double sn, cn;
        const double alpha = 0.5*(Gauss76Z[i]*M_PI_2 + M_PI_2);
        SINCOS(alpha, sn, cn);
        const double fq = _cyl(core_twovd, core_qr*sn, core_qh*cn)
            + _cyl(shell_twovd, shell_qr*sn, shell_qh*cn);
        total += Gauss76Wt[i] * fq * fq * sn;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    //const double form = (upper-lower)/2.0*total;
    return 1.0e-4 * total * M_PI_4;
}


double Iqxy(double qx, double qy,
    double core_sld,
    double shell_sld,
    double solvent_sld,
    double radius,
    double thickness,
    double length,
    double theta,
    double phi)
{
    double q, sin_alpha, cos_alpha;
    ORIENT_SYMMETRIC(qx, qy, theta, phi, q, sin_alpha, cos_alpha);

    const double core_qr = q*radius;
    const double core_qh = q*0.5*length;
    const double core_twovd = 2.0 * form_volume(radius,0,length)
                            * (core_sld-shell_sld);
    const double shell_qr = q*(radius + thickness);
    const double shell_qh = q*(0.5*length + thickness);
    const double shell_twovd = 2.0 * form_volume(radius,thickness,length)
                             * (shell_sld-solvent_sld);

    const double fq = _cyl(core_twovd, core_qr*sin_alpha, core_qh*cos_alpha)
        + _cyl(shell_twovd, shell_qr*sin_alpha, shell_qh*cos_alpha);
    return 1.0e-4 * fq * fq;
}
