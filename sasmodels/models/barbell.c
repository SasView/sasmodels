//barbell kernel - same as dumbbell
static double
_bell_kernel(double qab, double qc, double h, double radius_bell,
             double half_length)
{
    // translate a point in [-1,1] to a point in [lower,upper]
    const double upper = 1.0;
    const double lower = -h/radius_bell;
    const double zm = 0.5*(upper-lower);
    const double zb = 0.5*(upper+lower);

    // cos term in integral is:
    //    cos (q (R t - h + L/2) cos(alpha))
    // so turn it into:
    //    cos (m t + b)
    // where:
    //    m = q R cos(alpha)
    //    b = q(L/2-h) cos(alpha)
    const double m = radius_bell*qc; // cos argument slope
    const double b = (half_length+h)*qc; // cos argument intercept
    const double qab_r = radius_bell*qab; // Q*R*sin(theta)

    const double qr_max = fmax(qab_r, m);
    constant double *w, *z;
    int n = gauss_weights(qr_max, &w, &z);

    double total = 0.0;
    for (int i = 0; i < n; i++){
        const double t = z[i]*zm + zb;
        const double radical = 1.0 - t*t;
        const double bj = sas_2J1x_x(qab_r*sqrt(radical));
        const double Fq = cos(m*t + b) * radical * bj;
        total += w[i] * Fq;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    const double integral = total*zm;
    const double bell_fq = 2.0*M_PI*cube(radius_bell)*integral;
    return bell_fq;
}

static double
_fq(double qab, double qc, double h,
    double radius_bell, double radius, double half_length)
{
    const double bell_fq = _bell_kernel(qab, qc, h, radius_bell, half_length);
    const double bj = sas_2J1x_x(radius*qab);
    const double si = sas_sinx_x(half_length*qc);
    const double cyl_fq = 2.0*M_PI*radius*radius*half_length*bj*si;
    const double Aq = bell_fq + cyl_fq;
    return Aq;
}

static double
form_volume(double radius_bell,
    double radius,
    double length)
{
    // bell radius should never be less than radius when this is called
    const double h = sqrt(square(radius_bell) - square(radius));
    const double slice = M_PI*(square(radius_bell)*h - cube(h)/3.0);
    const double hemisphere = 2.0*M_PI/3.0*cube(radius_bell);
    const double rod = M_PI*square(radius)*length;
    // h > 0 so slice is added to hemisphere
    return rod + 2.0*(hemisphere + slice);
}

static double
radius_from_excluded_volume(double radius_bell, double radius, double length)
{
    const double h = sqrt(square(radius_bell) - square(radius));
    const double length_tot = length + 2.0*(radius + h);
    // Use cylinder excluded volume with length' = length + caps and
    // radius' = bell radius since the bell is bigger than the cylinder.
    return 0.5*cbrt(0.75*radius_bell*(2.0*radius_bell*length_tot
           + (radius_bell + length_tot)*(M_PI*radius_bell + length_tot)));
}

static double
radius_from_volume(double radius_bell, double radius, double length)
{
    const double vol_barbell = form_volume(radius_bell,radius,length);
    return cbrt(vol_barbell/M_4PI_3);
}

static double
radius_from_totallength(double radius_bell, double radius, double length)
{
    const double h = sqrt(square(radius_bell) - square(radius));
    const double half_length = 0.5*length;
    return half_length + radius_bell + h;
}

static double
radius_effective(int mode, double radius_bell, double radius, double length)
{
    switch (mode) {
    default:
    case 1: // equivalent cylinder excluded volume
        return radius_from_excluded_volume(radius_bell, radius , length);
    case 2: // equivalent volume sphere
        return radius_from_volume(radius_bell, radius , length);
    case 3: // radius
        return radius;
    case 4: // half length
        return 0.5*length;
    case 5: // half total length
        return radius_from_totallength(radius_bell,radius,length);
    }
}

static void
Fq(double q,double *F1, double *F2, double sld, double solvent_sld,
    double radius_bell, double radius, double length)
{
    const double h = sqrt(square(radius_bell) - square(radius));
    const double half_length = 0.5*length;

    const double qr_max = q*fmax(radius, half_length);
    constant double *w, *z;
    int n = gauss_weights(qr_max, &w, &z);

    // translate a point in [-1,1] to a point in [0, pi/2]
    const double zm = M_PI_4;
    const double zb = M_PI_4;
    double total_F1 = 0.0;
    double total_F2 = 0.0;
    for (int i = 0; i < n; i++){
        const double theta = z[i]*zm + zb;
        double sin_theta, cos_theta; // slots to hold sincos function output
        SINCOS(theta, sin_theta, cos_theta);
        const double qab = q*sin_theta;
        const double qc = q*cos_theta;
        const double Aq = _fq(qab, qc, h, radius_bell, radius, half_length);
        // scale by sin_theta for spherical coord integration
        total_F1 += w[i] * Aq * sin_theta;
        total_F2 += w[i] * Aq * Aq * sin_theta;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    const double form_avg = total_F1 * zm;
    const double form_squared_avg = total_F2 * zm;

    //Contrast
    const double s = (sld - solvent_sld);
    *F1 = 1.0e-2 * s * form_avg;
    *F2 = 1.0e-4 * s * s * form_squared_avg;
}

static double
Iqac(double qab, double qc,
    double sld, double solvent_sld,
    double radius_bell, double radius, double length)
{
    const double h = sqrt(square(radius_bell) - square(radius));
    const double Aq = _fq(qab, qc, h, radius_bell, radius, 0.5*length);

    // Multiply by contrast^2 and convert to cm-1
    const double s = (sld - solvent_sld);
    return 1.0e-4 * square(s * Aq);
}
