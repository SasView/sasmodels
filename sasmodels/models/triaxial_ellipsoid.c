static double
form_volume(double radius_equat_minor, double radius_equat_major, double radius_polar)
{
    return M_4PI_3*radius_equat_minor*radius_equat_major*radius_polar;
}

static double
radius_from_curvature(double radius_equat_minor, double radius_equat_major, double radius_polar)
{
    // Trivial cases
    if (radius_equat_minor == radius_equat_major == radius_polar) return radius_polar;
    if (radius_equat_minor * radius_equat_major * radius_polar == 0.)  return 0.;


    double r_equat_equiv, r_polar_equiv;
    double radii[3] = {radius_equat_minor, radius_equat_major, radius_polar};
    double radmax = fmax(radii[0],fmax(radii[1],radii[2]));

    double radius_1 = radmax;
    double radius_2 = radmax;
    double radius_3 = radmax;

    for(int irad=0; irad<3; irad++) {
        if (radii[irad] < radius_1) {
            radius_3 = radius_2;
            radius_2 = radius_1;
            radius_1 = radii[irad];
            } else {
                if (radii[irad] < radius_2) {
                        radius_2 = radii[irad];
                }
            }
    }
    if(radius_2-radius_1 > radius_3-radius_2) {
        r_equat_equiv = sqrt(radius_2*radius_3);
        r_polar_equiv = radius_1;
    } else  {
        r_equat_equiv = sqrt(radius_1*radius_2);
        r_polar_equiv = radius_3;
    }

    // see equation (26) in A.Isihara, J.Chem.Phys. 18(1950)1446-1449
    const double ratio = (r_polar_equiv < r_equat_equiv
                          ? r_polar_equiv / r_equat_equiv
                          : r_equat_equiv / r_polar_equiv);
    const double e1 = sqrt(1.0 - ratio*ratio);
    const double b1 = 1.0 + asin(e1) / (e1 * ratio);
    const double bL = (1.0 + e1) / (1.0 - e1);
    const double b2 = 1.0 + 0.5 * ratio * ratio / e1 * log(bL);
    const double delta = 0.75 * b1 * b2;
    const double ddd = 2.0 * (delta + 1.0) * r_polar_equiv * r_equat_equiv * r_equat_equiv;
    return 0.5 * cbrt(ddd);
}

static double
radius_from_volume(double radius_equat_minor, double radius_equat_major, double radius_polar)
{
    return cbrt(radius_equat_minor*radius_equat_major*radius_polar);
}

static double
radius_from_min_dimension(double radius_equat_minor, double radius_equat_major, double radius_polar)
{
    const double rad_equat_min = (radius_equat_minor < radius_equat_major ? radius_equat_minor : radius_equat_major);
    return (rad_equat_min < radius_polar ? rad_equat_min : radius_polar);
}

static double
radius_from_max_dimension(double radius_equat_minor, double radius_equat_major, double radius_polar)
{
    const double rad_equat_max = (radius_equat_minor < radius_equat_major ? radius_equat_major : radius_equat_minor);
    return (rad_equat_max > radius_polar ? rad_equat_max : radius_polar);
}

static double
radius_effective(int mode, double radius_equat_minor, double radius_equat_major, double radius_polar)
{
    switch (mode) {
    default:
    case 1: // equivalent biaxial ellipsoid average curvature
        return radius_from_curvature(radius_equat_minor,radius_equat_major, radius_polar);
    case 2: // equivalent volume sphere
        return radius_from_volume(radius_equat_minor,radius_equat_major, radius_polar);
    case 3: // min radius
        return radius_from_min_dimension(radius_equat_minor,radius_equat_major, radius_polar);
    case 4: // max radius
        return radius_from_max_dimension(radius_equat_minor,radius_equat_major, radius_polar);
    }
}

static void
Fq(double q,
    double *F1,
    double *F2,
    double sld,
    double sld_solvent,
    double radius_equat_minor,
    double radius_equat_major,
    double radius_polar)
{
    const double pa = square(radius_equat_minor/radius_equat_major) - 1.0;
    const double pc = square(radius_polar/radius_equat_major) - 1.0;

    // TODO: check that c_length dominates the outer loop, and a_length, b_length inner
    // const double qr_max = q*fmax(fmax(radius_equat_minor, radius_equat_major), radius_polar);
    const double qr_max = q*radius_polar;
    constant double *z, *w;
    int n = gauss_weights(qr_max, &w, &z);

    // translate a point in [-1,1] to a point in [0, pi/2]
    const double zm = M_PI_4;
    const double zb = M_PI_4;
    double outer_sum_F1 = 0.0;
    double outer_sum_F2 = 0.0;
    for (int i=0;i<n;i++) {
        //const double u = z[i]*(upper-lower)/2 + (upper + lower)/2;
        const double phi = z[i]*zm + zb;
        const double pa_sinsq_phi = pa*square(sin(phi));

        const double qr_max_inner = q*fmax(radius_equat_major, radius_equat_minor);
        constant double *z_inner, *w_inner;
        int n_inner = gauss_weights(qr_max_inner, &w_inner, &z_inner);

        double inner_sum_F1 = 0.0;
        double inner_sum_F2 = 0.0;
        const double um = 0.5;
        const double ub = 0.5;
        for (int j=0;j<n_inner;j++) {
            // translate a point in [-1,1] to a point in [0, 1]
            const double usq = square(z_inner[j]*um + ub);
            const double r = radius_equat_major*sqrt(pa_sinsq_phi*(1.0-usq) + 1.0 + pc*usq);
            const double fq = sas_3j1x_x(q*r);
            inner_sum_F1 += w_inner[j] * fq;
            inner_sum_F2 += w_inner[j] * fq * fq;
        }
        outer_sum_F1 += w[i] * inner_sum_F1;  // correcting for dx later
        outer_sum_F2 += w[i] * inner_sum_F2;  // correcting for dx later
    }
    // translate integration ranges from [-1,1] to [lower,upper] and normalize by 4 pi
    outer_sum_F1 *= 0.25;  // = outer*um*zm*8.0/(4.0*M_PI);
    outer_sum_F2 *= 0.25;  // = outer*um*zm*8.0/(4.0*M_PI);

    const double volume = form_volume(radius_equat_minor, radius_equat_major, radius_polar);
    const double contrast = (sld - sld_solvent);
    *F1 = 1.0e-2 * contrast * volume * outer_sum_F1;
    *F2 = 1.0e-4 * square(contrast * volume) * outer_sum_F2;
}


static double
Iqabc(double qa, double qb, double qc,
    double sld,
    double sld_solvent,
    double radius_equat_minor,
    double radius_equat_major,
    double radius_polar)
{
    const double qr = sqrt(square(radius_equat_minor*qa)
                           + square(radius_equat_major*qb)
                           + square(radius_polar*qc));
    const double fq = sas_3j1x_x(qr);
    const double vol = form_volume(radius_equat_minor, radius_equat_major, radius_polar);
    const double drho = (sld - sld_solvent);

    return 1.0e-4 * square(vol * drho * fq);
}
