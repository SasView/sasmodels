double form_volume(double radius_polar, double radius_equatorial);
double Iq(double q, double sld, double sld_solvent, double radius_polar, double radius_equatorial);
double Iqxy(double qx, double qy, double sld, double sld_solvent,
    double radius_polar, double radius_equatorial, double theta, double phi);

double _ellipsoid_kernel(double q, double radius_polar, double radius_equatorial, double sin_alpha);
double _ellipsoid_kernel(double q, double radius_polar, double radius_equatorial, double sin_alpha)
{
    double ratio = radius_polar/radius_equatorial;
    const double u = q*radius_equatorial*sqrt(1.0
                   + sin_alpha*sin_alpha*(ratio*ratio - 1.0));
    const double f = sph_j1c(u);

    return f*f;
}

double form_volume(double radius_polar, double radius_equatorial)
{
    return M_4PI_3*radius_polar*radius_equatorial*radius_equatorial;
}

double Iq(double q,
    double sld,
    double sld_solvent,
    double radius_polar,
    double radius_equatorial)
{
    // translate a point in [-1,1] to a point in [0, 1]
    const double zm = 0.5;
    const double zb = 0.5;
    double total = 0.0;
    for (int i=0;i<76;i++) {
        //const double sin_alpha = (Gauss76Z[i]*(upper-lower) + upper + lower)/2;
        const double sin_alpha = Gauss76Z[i]*zm + zb;
        total += Gauss76Wt[i] * _ellipsoid_kernel(q, radius_polar, radius_equatorial, sin_alpha);
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    const double form = total*zm;
    const double s = (sld - sld_solvent) * form_volume(radius_polar, radius_equatorial);
    return 1.0e-4 * s * s * form;
}

double Iqxy(double qx, double qy,
    double sld,
    double sld_solvent,
    double radius_polar,
    double radius_equatorial,
    double theta,
    double phi)
{
    double q, sin_alpha, cos_alpha;
    ORIENT_SYMMETRIC(qx, qy, theta, phi, q, sin_alpha, cos_alpha);
    const double form = _ellipsoid_kernel(q, radius_polar, radius_equatorial, sin_alpha);
    const double s = (sld - sld_solvent) * form_volume(radius_polar, radius_equatorial);

    return 1.0e-4 * form * s * s;
}

