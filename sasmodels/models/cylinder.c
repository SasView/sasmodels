double form_volume(double radius, double length);
double fq(double q, double sn, double cn,double radius, double length);
double orient_avg_1D(double q, double radius, double length);
double Iq(double q, double sld, double solvent_sld, double radius, double length);
double Iqxy(double qx, double qy, double sld, double solvent_sld,
    double radius, double length, double theta, double phi);

#define INVALID(v) (v.radius<0 || v.length<0)

double form_volume(double radius, double length)
{
    return M_PI*radius*radius*length;
}

double fq(double q, double sn, double cn, double radius, double length)
{
    // precompute qr and qh to save time in the loop
    const double qr = q*radius;
    const double qh = q*0.5*length; 
    return  sas_J1c(qr*sn) * sinc(qh*cn) ;
}

double orient_avg_1D(double q, double radius, double length)
{
    // translate a point in [-1,1] to a point in [0, pi/2]
    const double zm = M_PI_4;
    const double zb = M_PI_4; 

    double total = 0.0;
    for (int i=0; i<76 ;i++) {
        const double alpha = Gauss76Z[i]*zm + zb;
        double sn, cn; // slots to hold sincos function output
        // alpha(theta,phi) the projection of the cylinder on the detector plane
        SINCOS(alpha, sn, cn);
        total += Gauss76Wt[i] * fq(q, sn, cn, radius, length)* fq(q, sn, cn, radius, length) * sn;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    return total*zm;
}

double Iq(double q,
    double sld,
    double solvent_sld,
    double radius,
    double length)
{
    const double s = (sld - solvent_sld) * form_volume(radius, length);
    return 1.0e-4 * s * s * orient_avg_1D(q, radius, length);
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
    SINCOS(phi*M_PI_180, sn, cn);
    const double q = sqrt(qx*qx + qy*qy);
    const double cos_val = (q==0. ? 1.0 : (cn*qx + sn*qy)*sin(theta*M_PI_180)/q);

    const double alpha = acos(cos_val);

    SINCOS(alpha, sn, cn);
    const double s = (sld-solvent_sld) * form_volume(radius, length);
    return 1.0e-4 * square(s * fq(q, sn, cn, radius, length));
}
