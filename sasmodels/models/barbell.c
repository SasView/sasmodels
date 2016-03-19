double form_volume(double bell_radius, double radius, double length);
double Iq(double q, double sld, double solvent_sld,
        double bell_radius, double radius, double length);
double Iqxy(double qx, double qy, double sld, double solvent_sld,
        double bell_radius, double radius, double length,
        double theta, double phi);

//barbell kernel - same as dumbell
static double
_bell_kernel(double q, double h, double bell_radius,
             double half_length, double sin_alpha, double cos_alpha)
{
    // translate a point in [-1,1] to a point in [lower,upper]
    const double upper = 1.0;
    const double lower = h/bell_radius;
    const double zm = 0.5*(upper-lower);
    const double zb = 0.5*(upper+lower);

    // cos term in integral is:
    //    cos (q (R t - h + L/2) cos(alpha))
    // so turn it into:
    //    cos (m t + b)
    // where:
    //    m = q R cos(alpha)
    //    b = q(L/2-h) cos(alpha)
    const double m = q*bell_radius*cos_alpha; // cos argument slope
    const double b = q*(half_length-h)*cos_alpha; // cos argument intercept
    const double qrst = q*bell_radius*sin_alpha; // Q*R*sin(theta)
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
    const double bell_Fq = 2*M_PI*cube(bell_radius)*integral;
    return bell_Fq;
}

double form_volume(double bell_radius,
        double radius,
        double length)
{

    // bell radius should never be less than radius when this is called
    const double hdist = sqrt(square(bell_radius) - square(radius));
    const double p1 = 2.0/3.0*cube(bell_radius);
    const double p2 = square(bell_radius)*hdist;
    const double p3 = cube(hdist)/3.0;

    return M_PI*square(radius)*length + 2.0*M_PI*(p1+p2-p3);
}

double Iq(double q, double sld, double solvent_sld,
          double bell_radius, double radius, double length)
{
    // Exclude invalid inputs.
    if (bell_radius < radius) return NAN;
    const double h = -sqrt(bell_radius*bell_radius - radius*radius);
    const double half_length = 0.5*length;

    // translate a point in [-1,1] to a point in [0, pi/2]
    const double zm = M_PI_4;
    const double zb = M_PI_4;
    double total = 0.0;
    for (int i = 0; i < 76; i++){
        const double alpha= Gauss76Z[i]*zm + zb;
        double sin_alpha, cos_alpha; // slots to hold sincos function output
        SINCOS(alpha, sin_alpha, cos_alpha);

        const double bell_Fq = _bell_kernel(q, h, bell_radius, half_length, sin_alpha, cos_alpha);
        const double bj = sas_J1c(q*radius*sin_alpha);
        const double si = sinc(q*half_length*cos_alpha);
        const double cyl_Fq = M_PI*radius*radius*length*bj*si;
        const double Aq = bell_Fq + cyl_Fq;
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
        double bell_radius, double radius, double length,
        double theta, double phi)
{
    // Compute angle alpha between q and the cylinder axis
    double sn, cn; // slots to hold sincos function output
    SINCOS(theta*M_PI_180, sn, cn);
    const double q = sqrt(qx*qx+qy*qy);
    const double cos_val = cn*cos(phi*M_PI_180)*(qx/q) + sn*(qy/q);
    const double alpha = acos(cos_val); // rod angle relative to q

    // Exclude invalid inputs.
    if (bell_radius < radius) return NAN;
    const double h = -sqrt(square(bell_radius) - square(radius));
    const double half_length = 0.5*length;

    double sin_alpha, cos_alpha; // slots to hold sincos function output
    SINCOS(alpha, sin_alpha, cos_alpha);
    const double bell_Fq = _bell_kernel(q, h, bell_radius, half_length, sin_alpha, cos_alpha);
    const double bj = sas_J1c(q*radius*sin_alpha);
    const double si = sinc(q*half_length*cos_alpha);
    const double cyl_Fq = M_PI*radius*radius*length*bj*si;
    const double Aq = cyl_Fq + bell_Fq;

    // Multiply by contrast^2 and convert to cm-1
    const double s = (sld - solvent_sld);
    return 1.0e-4 * square(s * Aq);
}
