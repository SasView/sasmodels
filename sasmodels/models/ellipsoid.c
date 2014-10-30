double form_volume(double rpolar, double requatorial);
double Iq(double q, double sld, double solvent_sld, double rpolar, double requatorial);
double Iqxy(double qx, double qy, double sld, double solvent_sld,
    double rpolar, double requatorial, double theta, double phi);

double _ellipsoid_kernel(double q, double rpolar, double requatorial, double cos_alpha);
double _ellipsoid_kernel(double q, double rpolar, double requatorial, double cos_alpha)
{
    double sn, cn;
    double ratio = rpolar/requatorial;
    const double u = q*requatorial*sqrt(1.0
                   + cos_alpha*cos_alpha*(ratio*ratio - 1.0));
    SINCOS(u, sn, cn);
    const double f = ( u==0.0 ? 1.0 : 3.0*(sn-u*cn)/(u*u*u) );
    return f*f;
}

double form_volume(double rpolar, double requatorial)
{
    return 1.333333333333333*M_PI*rpolar*requatorial*requatorial;
}

double Iq(double q,
    double sld,
    double solvent_sld,
    double rpolar,
    double requatorial)
{
    //const double lower = 0.0;
    //const double upper = 1.0;
    double total = 0.0;
    for (int i=0;i<76;i++) {
        //const double cos_alpha = (Gauss76Z[i]*(upper-lower) + upper + lower)/2;
        const double cos_alpha = 0.5*(Gauss76Z[i] + 1.0);
        total += Gauss76Wt[i] * _ellipsoid_kernel(q, rpolar, requatorial, cos_alpha);
    }
    //const double form = (upper-lower)/2*total;
    const double form = 0.5*total;
    const double s = (sld - solvent_sld) * form_volume(rpolar, requatorial);
    return 1.0e-4 * form * s * s;
}

double Iqxy(double qx, double qy,
    double sld,
    double solvent_sld,
    double rpolar,
    double requatorial,
    double theta,
    double phi)
{
    double sn, cn;

    const double q = sqrt(qx*qx + qy*qy);
    SINCOS(theta*M_PI_180, sn, cn);
    const double cos_alpha = cn*cos(phi*M_PI_180)*(qx/q) + sn*(qy/q);
    const double form = _ellipsoid_kernel(q, rpolar, requatorial, cos_alpha);
    const double s = (sld - solvent_sld) * form_volume(rpolar, requatorial);

    return 1.0e-4 * form * s * s;
}

