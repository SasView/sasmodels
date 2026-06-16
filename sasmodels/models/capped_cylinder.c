// Integral over a convex lens kernel for t in [h/R,1].  See the docs for
// the definition of the function being integrated.
//   q is the magnitude of the q vector.
//   h is the length of the lens "inside" the cylinder.  This negative wrt the
//       definition of h in the docs.
//   radius_cap is the radius of the lens
//   length is the cylinder length, or the separation between the lens halves
//   theta is the angle of the cylinder wrt q.
static double
_cap_kernel(double qab, double qc, double h, double radius_cap, double radius,
    double half_length, int outer_n)
{
    // translate a point in [-1,1] to a point in [lower,upper]
    const double upper = 1.0;
    const double lower = -h/radius_cap; // integral lower bound
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
    const double b = (half_length+h)*qc; // cos argument intercept
    const double qab_r = radius_cap*qab; // Q*R*sin(theta)

    // m+b = qc*(half_length + radius_cap + h). With h in [-radius_cap, 0] depending
    // on cylinder radius, that means m+b is in qc*[length/2, length_2 + radius_cap].
    // The qab_r term will be very large for mostly flat caps. Since the bj term will
    // oscillate at this frequency, it seems like we should increase the number of
    // gauss points to accomodate. However, if we use the radius of the cylinder
    // we seem to get good results, so use that to set the number of integration points.
    //const double qr_max = fmax(qab_r, m+b);
    const double qr_max = fmax(qab*radius, m+b);
    constant double *w, *z;
    int n = gauss_weights(qr_max, outer_n, &w, &z);

    double total = 0.0;
    for (int i=0; i<n; i++) {
        const double t = z[i]*zm + zb;
        const double radical = 1.0 - t*t;
        const double bj = sas_2J1x_x(qab_r*sqrt(radical));
        const double Fq = cos(m*t + b) * radical * bj;
        total += w[i] * Fq;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    const double integral = total*zm;
    const double cap_Fq = 2.0*M_PI*cube(radius_cap)*integral;
    return cap_Fq;
}

static double
_fq(double qab, double qc, double h, double radius_cap, double radius, double half_length, int n_outer)
{
    const double cap_Fq = _cap_kernel(qab, qc, h, radius_cap, radius, half_length, n_outer);
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
    const double h = -sqrt(square(radius_cap) - square(radius));
    const double slice = M_PI*(square(radius_cap)*h - cube(h)/3.0);
    const double hemisphere = 2.0*M_PI/3.0*cube(radius_cap);
    const double rod = M_PI*square(radius)*length;
    // h < 0 so slice is subtracted from hemisphere
    return rod + 2.0*(hemisphere + slice);
}

static double
radius_from_excluded_volume(double radius, double radius_cap, double length)
{
    const double h = -sqrt(square(radius_cap) - square(radius));
    const double length_tot = length + 2.0*(radius_cap + h);
    // Use cylinder excluded volume with length' = length + caps and
    // radius' = cylinder radius since the lens is smaller than the cylinder.
    return 0.5*cbrt(0.75*radius*(2.0*radius*length_tot
           + (radius + length_tot)*(M_PI*radius + length_tot)));
}

static double
radius_from_volume(double radius, double radius_cap, double length)
{
    const double vol_cappedcyl = form_volume(radius,radius_cap,length);
    return cbrt(vol_cappedcyl/M_4PI_3);
}

static double
radius_from_totallength(double radius, double radius_cap, double length)
{
    const double h = -sqrt(square(radius_cap) - square(radius));
    const double half_length = 0.5*length;
    return half_length + radius_cap - h;
}

static double
radius_effective(int mode, double radius, double radius_cap, double length)
{
    switch (mode) {
    default:
    case 1: // equivalent cylinder excluded volume
        return radius_from_excluded_volume(radius, radius_cap, length);
    case 2: // equivalent volume sphere
        return radius_from_volume(radius, radius_cap, length);
    case 3: // radius
        return radius;
    case 4: // half length
        return 0.5*length;
    case 5: // half total length
        return radius_from_totallength(radius, radius_cap,length);
    }
}

static void
Fq(double q,double *F1, double *F2, double sld, double solvent_sld,
    double radius, double radius_cap, double length)
{
    const double h = -sqrt(square(radius_cap) - square(radius));
    const double half_length = 0.5*length;

    // The term h comes from solving the right triangle with diagonal
    // equal to the cap radius and horizontal equal to the bar radius.
    // The result is the (negative) height of the equator above the end of the rod.
    // To get the total length of bar+cap use bar length + 2*(cap radius + h).
    // We want the radius, so divide that by two.
    // For a lentil with length=0, the radius will be the dominant term, hence
    // the fmax in the calculation below. This isn't needed for the barbell shape
    // since the bell length is always greater than the bar radius.
    const double qr_max = q*fmax(half_length + radius_cap + h, radius);
    //const double qr_max = q*(half_length + radius_cap + h);
    constant double *w_outer, *z_outer;
    // Keep outer loop to 76 or less
    int n_outer = gauss_weights(qr_max, ADAPTIVE_MAX_76, &w_outer, &z_outer);

    // translate a point in [-1,1] to a point in [0, pi/2]
    const double zm = M_PI_4;
    const double zb = M_PI_4;
    double total_F1 = 0.0;
    double total_F2 = 0.0;
    for (int i=0; i<n_outer ;i++) {
        const double theta = z_outer[i]*zm + zb;
        double sin_theta, cos_theta; // slots to hold sincos function output
        SINCOS(theta, sin_theta, cos_theta);
        const double qab = q*sin_theta;
        const double qc = q*cos_theta;
        // Don't constrain size of inner loop to n_outer
        const double Aq = _fq(qab, qc, h, radius_cap, radius, half_length, 1);
        // scale by sin_theta for spherical coord integration
        total_F1 += w_outer[i] * Aq * sin_theta;
        total_F2 += w_outer[i] * Aq * Aq * sin_theta;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    const double form_avg = total_F1 * zm;
    const double form_squared_avg = total_F2 * zm;

    // Contrast
    const double s = (sld - solvent_sld);
    *F1 = 1.0e-2 * s * form_avg;
    *F2 = 1.0e-4 * s * s * form_squared_avg;
}


static double
Iqac(double qab, double qc,
    double sld, double solvent_sld, double radius,
    double radius_cap, double length)
{
    // TODO: For 2D data we may also want to limit the size of the integral at each (qx, qy)
    const int n_outer = 1; // No limits on inner integral
    const double h = -sqrt(square(radius_cap) - square(radius));
    const double Aq = _fq(qab, qc, h, radius_cap, radius, 0.5*length, n_outer);

    // Multiply by contrast^2 and convert to cm-1
    const double s = (sld - solvent_sld);
    return 1.0e-4 * square(s * Aq);
}
