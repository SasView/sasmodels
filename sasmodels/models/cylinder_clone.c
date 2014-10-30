double form_volume(double radius, double length);
double Iq(double q, double sld, double solvent_sld, double radius, double length);
double Iqxy(double qx, double qy, double sld, double solvent_sld, double radius, double length, double theta, double phi);


// twovd = 2 * volume * delta_rho
// besarg = q * R * sin(alpha)
// siarg = q * L/2 * cos(alpha)
double _cyl(double twovd, double besarg, double siarg, double alpha);
double _cyl(double twovd, double besarg, double siarg, double alpha)
{
    const double bj = (besarg == 0.0 ? 0.5 : J1(besarg)/besarg);
    const double si = (siarg == 0.0 ? 1.0 : sin(siarg)/siarg);
    return twovd*si*bj;
}

double form_volume(double radius, double length)
{
    return M_PI*radius*radius*length;
}
double Iq(double q,
    double sldCyl,
    double sldSolv,
    double radius,
    double length)
{
    const double qr = q*radius;
    const double qh = q*0.5*length;
    const double twovd = 2.0*(sldCyl-sldSolv)*form_volume(radius, length);
    double total = 0.0;
    // double lower=0, upper=M_PI_2;
    for (int i=0; i<76 ;i++) {
        // translate a point in [-1,1] to a point in [lower,upper]
        //const double alpha = ( Gauss76Z[i]*(upper-lower) + upper + lower )/2.0;
        const double alpha = 0.5*(Gauss76Z[i]*M_PI_2 + M_PI_2);
        double sn, cn;
        SINCOS(alpha, sn, cn);
        const double fq = _cyl(twovd, qr*sn, qh*cn, alpha);
        total += Gauss76Wt[i] * fq * fq * sn;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    //const double form = (upper-lower)/2.0*total;
    return 1.0e8 * total * M_PI_4;
}

double Iqxy(double qx, double qy,
    double sldCyl,
    double sldSolv,
    double radius,
    double length,
    double cyl_theta,
    double cyl_phi)
{
    double sn, cn; // slots to hold sincos function output

    // Compute angle alpha between q and the cylinder axis
    SINCOS(cyl_theta*M_PI_180, sn, cn);
    // # The following correction factor exists in sasview, but it can't be
    // # right, so we are leaving it out for now.
    const double spherical_integration = fabs(cn)*M_PI_2;
    const double q = sqrt(qx*qx+qy*qy);
    const double cos_val = cn*cos(cyl_phi*M_PI_180)*(qx/q) + sn*(qy/q);
    const double alpha = acos(cos_val);

    const double qr = q*radius;
    const double qh = q*0.5*length;
    const double twovd = 2.0*(sldCyl-sldSolv)*form_volume(radius, length);
    SINCOS(alpha, sn, cn);
    const double fq = _cyl(twovd, qr*sn, qh*cn, alpha);
    return 1.0e8 * fq * fq * spherical_integration;
}
