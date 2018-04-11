//#define INVALID(v) (v.radius_equat_minor > v.radius_equat_major || v.radius_equat_major > v.radius_polar)

static double
form_volume(double radius_equat_minor, double radius_equat_major, double radius_polar)
{
    return M_4PI_3*radius_equat_minor*radius_equat_major*radius_polar;
}

static double
Iq(double q,
    double sld,
    double sld_solvent,
    double radius_equat_minor,
    double radius_equat_major,
    double radius_polar)
{
    const double pa = square(radius_equat_minor/radius_equat_major) - 1.0;
    const double pc = square(radius_polar/radius_equat_major) - 1.0;
    // translate a point in [-1,1] to a point in [0, pi/2]
    const double zm = M_PI_4;
    const double zb = M_PI_4;
    double outer = 0.0;
    for (int i=0;i<GAUSS_N;i++) {
        //const double u = GAUSS_Z[i]*(upper-lower)/2 + (upper + lower)/2;
        const double phi = GAUSS_Z[i]*zm + zb;
        const double pa_sinsq_phi = pa*square(sin(phi));

        double inner = 0.0;
        const double um = 0.5;
        const double ub = 0.5;
        for (int j=0;j<GAUSS_N;j++) {
            // translate a point in [-1,1] to a point in [0, 1]
            const double usq = square(GAUSS_Z[j]*um + ub);
            const double r = radius_equat_major*sqrt(pa_sinsq_phi*(1.0-usq) + 1.0 + pc*usq);
            const double fq = sas_3j1x_x(q*r);
            inner += GAUSS_W[j] * fq * fq;
        }
        outer += GAUSS_W[i] * inner;  // correcting for dx later
    }
    // translate integration ranges from [-1,1] to [lower,upper] and normalize by 4 pi
    const double fqsq = outer/4.0;  // = outer*um*zm*8.0/(4.0*M_PI);
    const double vol = form_volume(radius_equat_minor, radius_equat_major, radius_polar);
    const double drho = (sld - sld_solvent);
    return 1.0e-4 * square(vol*drho) * fqsq;
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
