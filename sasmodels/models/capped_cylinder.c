double form_volume(double radius, double cap_radius, double length);
double Iq(double q, double sld, double solvent_sld,
    double radius, double cap_radius, double length);
double Iqxy(double qx, double qy, double sld, double solvent_sld,
    double radius, double cap_radius, double length, double theta, double phi);

// Integral over a convex lens kernel for t in [h/R,1].  See the docs for
// the definition of the function being integrated.
//   q is the magnitude of the q vector.
//   h is the length of the lens "inside" the cylinder.  This negative wrt the
//       definition of h in the docs.
//   cap_radius is the radius of the lens
//   length is the cylinder length, or the separation between the lens halves
//   alpha is the angle of the cylinder wrt q.
double _cap_kernel(double q, double h, double cap_radius, double length,
                 double sin_alpha, double cos_alpha);
double _cap_kernel(double q, double h, double cap_radius, double length,
                 double sin_alpha, double cos_alpha)
{
    // For speed, we are pre-calculating terms which are constant over the
    // kernel.
    const double upper = 1.0;
    const double lower = h/cap_radius; // integral lower bound
    // cos term in integral is:
    //    cos (q (R t - h + L/2) cos(alpha))
    // so turn it into:
    //    cos (m t + b)
    // where:
    //    m = q R cos(alpha)
    //    b = q(L/2-h) cos(alpha)
    const double m = q*cap_radius*cos_alpha; // cos argument slope
    const double b = q*(0.5*length-h)*cos_alpha; // cos argument intercept
    const double qrst = q*cap_radius*sin_alpha; // Q*R*sin(theta)
    double total = 0.0;
    for (int i=0; i<76 ;i++) {
        // translate a point in [-1,1] to a point in [lower,upper]
        //const double t = ( Gauss76Z[i]*(upper-lower) + upper + lower )/2.0;
        const double t = 0.5*(Gauss76Z[i]*(upper-lower)+upper+lower);
        const double radical = 1.0 - t*t;
        const double arg = qrst*sqrt(radical); // cap bessel function arg
        const double be = (arg == 0.0 ? 0.5 : J1(arg)/arg);
        const double Fq = cos(m*t + b) * radical * be;
        total += Gauss76Wt[i] * Fq;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    //const double form = (upper-lower)/2.0*total;
    const double integral = 0.5*(upper-lower)*total;
    return 4.0*M_PI*cap_radius*cap_radius*cap_radius*integral;
}

double form_volume(double radius, double cap_radius, double length)
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
    //    (2R+h) => 2R+ h_c-h_c + h => 2R + (R-h)-hc + h => 3R-h_c
    // And so:
    //    2/3 pi (R-h)^2 (2R + h) => 2/3 pi h_c^2 (3 r_cap - h_c)
    // which is 2 V_cap, using the second form above.
    //
    // In this function we are going to use the first form of V_cap
    //
    //      V = V_cyl + 2 V_cap
    //        = pi r^2 L + pi hc (r^2 + hc^2/3)
    //        = pi (r^2 (L+hc) + hc^3/3)
    const double hc = cap_radius - sqrt(cap_radius*cap_radius - radius*radius);
    return M_PI*(radius*radius*(length+hc) + 0.333333333333333*hc*hc*hc);
}

double Iq(double q,
    double sld,
    double solvent_sld,
    double radius,
    double cap_radius,
    double length)
{
    double sn, cn; // slots to hold sincos function output

    // Exclude invalid inputs.
    if (cap_radius < radius) return -1.0;

    const double lower = 0.0;
    const double upper = M_PI_2;
    const double h = sqrt(cap_radius*cap_radius - radius*radius); // negative h
    double total = 0.0;
    for (int i=0; i<76 ;i++) {
        // translate a point in [-1,1] to a point in [lower,upper]
        const double alpha= 0.5*(Gauss76Z[i]*(upper-lower) + upper + lower);
        SINCOS(alpha, sn, cn);

        const double cap_Fq = _cap_kernel(q, h, cap_radius, length, sn, cn);

        // The following is CylKernel() / sin(alpha), but we are doing it in place
        // to avoid sin(alpha)/sin(alpha) for alpha = 0.  It is also a teensy bit
        // faster since we don't multiply and divide sin(alpha).
        const double besarg = q*radius*sn;
        const double siarg = q*0.5*length*cn;
        // lim_{x->0} J1(x)/x = 1/2,   lim_{x->0} sin(x)/x = 1
        const double bj = (besarg == 0.0 ? 0.5 : J1(besarg)/besarg);
        const double si = (siarg == 0.0 ? 1.0 : sin(siarg)/siarg);
        const double cyl_Fq = M_PI*radius*radius*length*2.0*bj*si;

        // Volume weighted average F(q)
        const double Aq = cyl_Fq + cap_Fq;
        total += Gauss76Wt[i] * Aq * Aq * sn; // sn for spherical coord integration
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    const double form = total * 0.5*(upper-lower);

    // Multiply by contrast^2, normalize by cylinder volume and convert to cm-1
    // NOTE that for this (Fournet) definition of the integral, one must MULTIPLY by Vcyl
    // The additional volume factor is for polydisperse volume normalization.
    const double s = (sld - solvent_sld);
    return 1.0e-4 * form * s * s; // form_volume(radius, cap_radius, length);
}


double Iqxy(double qx, double qy,
    double sld,
    double solvent_sld,
    double radius,
    double cap_radius,
    double length,
    double theta,
    double phi)
{
    double sn, cn; // slots to hold sincos function output

    // Exclude invalid inputs.
    if (cap_radius < radius) return -1.0;

    // Compute angle alpha between q and the cylinder axis
    SINCOS(theta*M_PI_180, sn, cn);
    // # The following correction factor exists in sasview, but it can't be
    // # right, so we are leaving it out for now.
    const double q = sqrt(qx*qx+qy*qy);
    const double cos_val = cn*cos(phi*M_PI_180)*(qx/q) + sn*(qy/q);
    const double alpha = acos(cos_val); // rod angle relative to q
    SINCOS(alpha, sn, cn);

    const double h = sqrt(cap_radius*cap_radius - radius*radius); // negative h
    const double cap_Fq = _cap_kernel(q, h, cap_radius, length, sn, cn);

    const double besarg = q*radius*sn;
    const double siarg = q*0.5*length*cn;
    // lim_{x->0} J1(x)/x = 1/2,   lim_{x->0} sin(x)/x = 1
    const double bj = (besarg == 0.0 ? 0.5 : J1(besarg)/besarg);
    const double si = (siarg == 0.0 ? 1.0 : sin(siarg)/siarg);
    const double cyl_Fq = M_PI*radius*radius*length*2.0*bj*si;

    // Volume weighted average F(q)
    const double Aq = cyl_Fq + cap_Fq;

    // Multiply by contrast^2, normalize by cylinder volume and convert to cm-1
    const double s = (sld - solvent_sld);
    return 1.0e-4 * Aq * Aq * s * s; // form_volume(radius, cap_radius, length);
}
