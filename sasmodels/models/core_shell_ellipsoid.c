
// Converted from Igor function gfn4, using the same pattern as ellipsoid
// for evaluating the parts of the integral.
//     FUNCTION gfn4:    CONTAINS F(Q,A,B,MU)**2  AS GIVEN
//                       BY (53) & (58-59) IN CHEN AND
//                       KOTLARCHYK REFERENCE
// see functions in /lib/cs_ellipsoid_funcs.c

static double
form_volume(double radius_equat_core,
    double x_core,
    double thick_shell,
    double x_polar_shell)
{
    const double equat_shell = radius_equat_core + thick_shell;
    const double polar_shell = radius_equat_core*x_core + thick_shell*x_polar_shell;
    double vol = M_4PI_3*equat_shell*equat_shell*polar_shell;
    return vol;
}

static double
radius_from_volume(double radius_equat_core, double x_core, double thick_shell, double x_polar_shell)
{
    const double volume_ellipsoid = form_volume(radius_equat_core, x_core, thick_shell, x_polar_shell);
    return cbrt(volume_ellipsoid/M_4PI_3);
}

static double
radius_effective(int mode, double radius_equat_core, double x_core, double thick_shell, double x_polar_shell)
{
    const double radius_equat_tot = radius_equat_core + thick_shell;
    const double radius_polar_tot = radius_equat_core*x_core + thick_shell*x_polar_shell;

    switch (mode) {
    default:
    case 1: // average outer curvature
        return radius_from_curvature(radius_equat_core, x_core, thick_shell, x_polar_shell);
    case 2: // equivalent volume sphere
        return radius_from_volume(radius_equat_core, x_core, thick_shell, x_polar_shell);
    case 3: // min outer radius
        return (radius_polar_tot < radius_equat_tot ? radius_polar_tot : radius_equat_tot);
    case 4: // max outer radius
        return (radius_polar_tot > radius_equat_tot ? radius_polar_tot : radius_equat_tot);
    }
}

static void
Fq(double q,
    double *F1,
    double *F2,
    double radius_equat_core,
    double x_core,
    double thick_shell,
    double x_polar_shell,
    double core_sld,
    double shell_sld,
    double solvent_sld)
{
    const double sld_core_shell = core_sld - shell_sld;
    const double sld_shell_solvent = shell_sld - solvent_sld;

    const double polar_core = radius_equat_core*x_core;
    const double equat_shell = radius_equat_core + thick_shell;
    const double polar_shell = radius_equat_core*x_core + thick_shell*x_polar_shell;

    // translate from [-1, 1] => [0, 1]
    const double m = 0.5;
    const double b = 0.5;
    double total_F1 = 0.0;     //initialize intergral
    double total_F2 = 0.0;     //initialize intergral
    for(int i=0;i<GAUSS_N;i++) {
        const double cos_theta = GAUSS_Z[i]*m + b;
        const double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
        double fq = _cs_ellipsoid_kernel(q*sin_theta, q*cos_theta,
            radius_equat_core, polar_core,
            equat_shell, polar_shell,
            sld_core_shell, sld_shell_solvent);
        total_F1 += GAUSS_W[i] * fq;
        total_F2 += GAUSS_W[i] * fq * fq;
    }
    total_F1 *= m;
    total_F2 *= m;

    // convert to [cm-1]
    *F1 = 1.0e-2 * total_F1;
    *F2 = 1.0e-4 * total_F2;
}


static double
Iqac(double qab, double qc,
    double radius_equat_core,
    double x_core,
    double thick_shell,
    double x_polar_shell,
    double core_sld,
    double shell_sld,
    double solvent_sld)
{
    const double sld_core_shell = core_sld - shell_sld;
    const double sld_shell_solvent = shell_sld - solvent_sld;

    const double polar_core = radius_equat_core*x_core;
    const double equat_shell = radius_equat_core + thick_shell;
    const double polar_shell = radius_equat_core*x_core + thick_shell*x_polar_shell;

    double fq = _cs_ellipsoid_kernel(qab, qc,
                  radius_equat_core, polar_core,
                  equat_shell, polar_shell,
                  sld_core_shell, sld_shell_solvent);

    //convert to [cm-1]
    return 1.0e-4 * fq * fq;
}
