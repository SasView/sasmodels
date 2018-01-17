#define INVALID(v) (v.radius_cap < v.radius)

// Integral over a convex lens kernel for t in [h/R,1].  See the docs for
// the definition of the function being integrated.
//   q is the magnitude of the q vector.
//   h is the length of the lens "inside" the cylinder.  This negative wrt the
//       definition of h in the docs.
//   radius_cap is the radius of the lens
//   length is the cylinder length, or the separation between the lens halves
//   theta is the angle of the cylinder wrt q.
static double
_cap_kernel(double qab, double qc, double h, double radius_cap,
    double half_length)
{
    // translate a point in [-1,1] to a point in [lower,upper]
    const double upper = 1.0;
    const double lower = h/radius_cap; // integral lower bound
    const double zm = 0.5*(upper-lower);
    const double zb = 0.5*(upper+lower);

    // cos term in integral is:
    //    cos (q (R t - h + L/2) cos(theta))
    // so turn it into:
    //    cos (m t + b)
    // where:
    //    m = q R cos(theta)
    //    b = q(L/2-h) cos(theta)
    const double m = radius_cap*qc; // cos argument slope
    const double b = (half_length-h)*qc; // cos argument intercept
    const double qab_r = radius_cap*qab; // Q*R*sin(theta)
    double total = 0.0;
    for (int i=0; i<GAUSS_N; i++) {
        const double t = GAUSS_Z[i]*zm + zb;
        const double radical = 1.0 - t*t;
        const double bj = sas_2J1x_x(qab_r*sqrt(radical));
        const double Fq = cos(m*t + b) * radical * bj;
        total += GAUSS_W[i] * Fq;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    const double integral = total*zm;
    const double cap_Fq = 2.0*M_PI*cube(radius_cap)*integral;
    return cap_Fq;
}

static double
_fq(double qab, double qc, double h, double radius_cap, double radius, double half_length)
{
    const double cap_Fq = _cap_kernel(qab, qc, h, radius_cap, half_length);
    const double bj = sas_2J1x_x(radius*qab);
    const double si = sas_sinx_x(half_length*qc);
    const double cyl_Fq = 2.0*M_PI*radius*radius*half_length*bj*si;
    const double Aq = cap_Fq + cyl_Fq;
    return Aq;
}

static double
form_volume(double radius, double radius_cap, double length)
{
    // cap radius should never be less than radius when this is called

    // Note: volume V = 2*V_cap + V_cyl
    //
    // V_cyl = pi r_cyl^2 L
    // V_cap = 1/6 pi h_c (3 r_cyl^2 + h_c^2) = 1/3 pi h_c^2 (3 r_cap - h_c)
    //
    // The docs for capped cylinder give the volume as:
    //    V = pi r^2 L + 2/3 pi (R-h)^2 (2R + h)
    // where r_cap=R and h = R - h_c.
    //
    // The first part is clearly V_cyl.  The second part requires some work:
    //    (R-h)^2 => h_c^2
    //    (2R+h) => 2R+ h_c-h_c + h => 2R + (R-h)-h_c + h => 3R-h_c
    // And so:
    //    2/3 pi (R-h)^2 (2R + h) => 2/3 pi h_c^2 (3 r_cap - h_c)
    // which is 2 V_cap, using the second form above.
    //
    // In this function we are going to use the first form of V_cap
    //
    //      V = V_cyl + 2 V_cap
    //        = pi r^2 L + pi hc (r^2 + hc^2/3)
    //        = pi (r^2 (L+hc) + hc^3/3)
    const double hc = radius_cap - sqrt(radius_cap*radius_cap - radius*radius);
    return M_PI*(radius*radius*(length+hc) + hc*hc*hc/3.0);
}

static double
Iq(double q, double sld, double solvent_sld,
    double radius, double radius_cap, double length)
{
    const double h = sqrt(radius_cap*radius_cap - radius*radius);
    const double half_length = 0.5*length;

    // translate a point in [-1,1] to a point in [0, pi/2]
    const double zm = M_PI_4;
    const double zb = M_PI_4;
    double total = 0.0;
    for (int i=0; i<GAUSS_N ;i++) {
        const double theta = GAUSS_Z[i]*zm + zb;
        double sin_theta, cos_theta; // slots to hold sincos function output
        SINCOS(theta, sin_theta, cos_theta);
        const double qab = q*sin_theta;
        const double qc = q*cos_theta;
        const double Aq = _fq(qab, qc, h, radius_cap, radius, half_length);
        // scale by sin_theta for spherical coord integration
        total += GAUSS_W[i] * Aq * Aq * sin_theta;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    const double form = total * zm;

    // Contrast
    const double s = (sld - solvent_sld);
    return 1.0e-4 * s * s * form;
}


static double
Iqac(double qab, double qc,
    double sld, double solvent_sld, double radius,
    double radius_cap, double length)
{
    const double h = sqrt(radius_cap*radius_cap - radius*radius);
    const double Aq = _fq(qab, qc, h, radius_cap, radius, 0.5*length);

    // Multiply by contrast^2 and convert to cm-1
    const double s = (sld - solvent_sld);
    return 1.0e-4 * square(s * Aq);
}
