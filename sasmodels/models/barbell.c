double form_volume(double radius_bell, double radius, double length);
double Iq(double q, double sld, double solvent_sld,
        double radius_bell, double radius, double length);
double Iqxy(double qx, double qy, double sld, double solvent_sld,
        double radius_bell, double radius, double length,
        double theta, double phi);

#define INVALID(v) (v.radius_bell < v.radius)

//barbell kernel - same as dumbell
static double
_bell_kernel(double q, double h, double radius_bell,
             double half_length, double sin_alpha, double cos_alpha)
{
    // translate a point in [-1,1] to a point in [lower,upper]
    const double upper = 1.0;
    const double lower = h/radius_bell;
    const double zm = 0.5*(upper-lower);
    const double zb = 0.5*(upper+lower);

    // cos term in integral is:
    //    cos (q (R t - h + L/2) cos(alpha))
    // so turn it into:
    //    cos (m t + b)
    // where:
    //    m = q R cos(alpha)
    //    b = q(L/2-h) cos(alpha)
    const double m = q*radius_bell*cos_alpha; // cos argument slope
    const double b = q*(half_length-h)*cos_alpha; // cos argument intercept
    const double qrst = q*radius_bell*sin_alpha; // Q*R*sin(theta)
    double total = 0.0;
    for (int i = 0; i < 76; i++){
        const double t = Gauss76Z[i]*zm + zb;
        const double radical = 1.0 - t*t;
        const double bj = sas_J1c(qrst*sqrt(radical));
        const double Fq = cos(m*t + b) * radical * bj;
        total += Gauss76Wt[i] * Fq;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    const double integral = total*zm;
    const double bell_fq = 2*M_PI*cube(radius_bell)*integral;
    return bell_fq;
}

static double
_fq(double q, double h,
    double radius_bell, double radius, double half_length,
    double sin_alpha, double cos_alpha)
{
    const double bell_fq = _bell_kernel(q, h, radius_bell, half_length, sin_alpha, cos_alpha);
    const double bj = sas_J1c(q*radius*sin_alpha);
    const double si = sinc(q*half_length*cos_alpha);
    const double cyl_fq = 2.0*M_PI*radius*radius*half_length*bj*si;
    const double Aq = bell_fq + cyl_fq;
    return Aq;
}


double form_volume(double radius_bell,
        double radius,
        double length)
{
    // bell radius should never be less than radius when this is called
    const double hdist = sqrt(square(radius_bell) - square(radius));
    const double p1 = 2.0/3.0*cube(radius_bell);
    const double p2 = square(radius_bell)*hdist;
    const double p3 = cube(hdist)/3.0;

    return M_PI*square(radius)*length + 2.0*M_PI*(p1+p2-p3);
}

double Iq(double q, double sld, double solvent_sld,
          double radius_bell, double radius, double length)
{
    const double h = -sqrt(radius_bell*radius_bell - radius*radius);
    const double half_length = 0.5*length;

    // translate a point in [-1,1] to a point in [0, pi/2]
    const double zm = M_PI_4;
    const double zb = M_PI_4;
    double total = 0.0;
    for (int i = 0; i < 76; i++){
        const double alpha= Gauss76Z[i]*zm + zb;
        double sin_alpha, cos_alpha; // slots to hold sincos function output
        SINCOS(alpha, sin_alpha, cos_alpha);
        const double Aq = _fq(q, h, radius_bell, radius, half_length, sin_alpha, cos_alpha);
        total += Gauss76Wt[i] * Aq * Aq * sin_alpha;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    const double form = total*zm;

    //Contrast
    const double s = (sld - solvent_sld);
    return 1.0e-4 * s * s * form;
}


double Iqxy(double qx, double qy,
        double sld, double solvent_sld,
        double radius_bell, double radius, double length,
        double theta, double phi)
{
    double q, sin_alpha, cos_alpha;
    ORIENT_SYMMETRIC(qx, qy, theta, phi, q, sin_alpha, cos_alpha);

    const double h = -sqrt(square(radius_bell) - square(radius));
    const double Aq = _fq(q, h, radius_bell, radius, 0.5*length, sin_alpha, cos_alpha);

    // Multiply by contrast^2 and convert to cm-1
    const double s = (sld - solvent_sld);
    return 1.0e-4 * square(s * Aq);
}
