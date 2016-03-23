double form_volume(double radius, double length);
double Iq(double q, double sld, double solvent_sld, double radius, double length);
double Iqxy(double qx, double qy, double sld, double solvent_sld,
    double radius, double length, double theta, double phi);

#define INVALID(v) (v.radius<0 || v.length<0)

double form_volume(double radius, double length)
{
    return M_PI*radius*radius*length;
}

double Iq(double q,
    double sld,
    double solvent_sld,
    double radius,
    double length)
{
    // precompute qr and qh to save time in the loop
    const double qr = q*radius;
    const double qh = q*0.5*length;

    // translate a point in [-1,1] to a point in [0, pi/2]
    const double zm = M_PI_4;
    const double zb = M_PI_4;

    double total = 0.0;
    for (int i=0; i<76 ;i++) {
        const double alpha = Gauss76Z[i]*zm + zb;
        double sn, cn;
        SINCOS(alpha, sn, cn);
        const double fq = sinc(qh*cn) * sas_J1c(qr*sn);
        total += Gauss76Wt[i] * fq*fq * sn;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    const double form = total*zm;
    const double s = (sld - solvent_sld) * form_volume(radius, length);
    return 1.0e-4 * s * s * form;
}


double Iqxy(double qx, double qy,
    double sld,
    double solvent_sld,
    double radius,
    double length,
    double theta,
    double phi)
{
    double sn, cn; // slots to hold sincos function output

    // Compute angle alpha between q and the cylinder axis
    SINCOS(theta*M_PI_180, sn, cn);
    const double q = sqrt(qx*qx + qy*qy);
    const double cos_val = (q==0. ? 1.0 : (cn*cos(phi*M_PI_180)*qx + sn*qy)/q);
    const double alpha = acos(cos_val);

    SINCOS(alpha, sn, cn);
    const double fq = sinc(q*0.5*length*cn) * sas_J1c(q*radius*sn);
    const double s = (sld-solvent_sld) * form_volume(radius, length);
    return 1.0e-4 * square(s * fq);
}
