// tied core/shell ellipsoid, from core_shell_ellipsoid in 4.2.2, RKH 17/07/2019
// RKH 19/09/2023, rewrite for v5 to give Fq and effective radius etc.
// Some common functions now in ../lib/cs_ellipsoid_funcs.c
// though beware of overlapping form_volume functions!

// Converted from Igor function gfn4, using the same pattern as ellipsoid
// for evaluating the parts of the integral.
//     FUNCTION gfn4:    CONTAINS F(Q,A,B,MU)**2  AS GIVEN
//                       BY (53) & (58-59) IN CHEN AND
//                       KOTLARCHYK REFERENCE
//

static double
// compiler needs ALL params here to by of type "volume" in the .py definition, even though we may not really want to be 
// able to makle them all polydisperse.
form_volume(double radius_equat_core,
    double x_core,
    double vol_dry_shell_over_core,
    double x_polar_shell,
    double f_solvent_in_shell)
{
    // Now we have to set up and solve the one real root of a cubic equation for thick_shell
    const double thick_shell = solve_shell_v3(radius_equat_core, x_core, vol_dry_shell_over_core,
    x_polar_shell, fabs(f_solvent_in_shell));
    const double equat_shell = radius_equat_core + thick_shell;
    const double polar_shell = radius_equat_core*x_core + thick_shell*x_polar_shell;
    double vol = M_4PI_3*equat_shell*equat_shell*polar_shell;
    return vol;
}
static double
radius_from_volume(double radius_equat_core, double x_core, double thick_shell, double x_polar_shell)
{
    const double equat_shell = radius_equat_core + thick_shell;
    const double polar_shell = radius_equat_core*x_core + thick_shell*x_polar_shell;
    const double vol = M_4PI_3*equat_shell*equat_shell*polar_shell;

    return cbrt(vol/M_4PI_3);
}

static double
radius_effective(int mode, double radius_equat_core, double x_core, 
    double vol_dry_shell_over_core,
    double x_polar_shell,
    double f_solvent_in_shell)
{
    // Now we have to set up and solve the one real root of a cubic equation for thick_shell
    const double thick_shell = solve_shell_v3(radius_equat_core, x_core, vol_dry_shell_over_core,
    x_polar_shell, fabs(f_solvent_in_shell));
    const double equat_shell = radius_equat_core + thick_shell;
    const double polar_shell = radius_equat_core*x_core + thick_shell*x_polar_shell;

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

static double
Fq(double q, double *F1, double *F2,
    double radius_equat_core,
    double x_core,
    double vol_dry_shell_over_core,
    double x_polar_shell,
    double core_sld,
    double dry_shell_sld,
    double solvent_sld,
    double f_solvent_in_shell)
{
    // sort out the parameterized sld
    const double f_solvent = fabs(f_solvent_in_shell);
    const double shell_sld = f_solvent*solvent_sld + (1-f_solvent)*dry_shell_sld;
    const double sld_core_shell = core_sld - shell_sld;
    const double sld_shell_solvent = shell_sld - solvent_sld;
    // Now we have to set up and solve the one real root of a cubic equation for thick_shell
    const double thick_shell = solve_shell_v3(radius_equat_core, x_core, vol_dry_shell_over_core,
    x_polar_shell, f_solvent);
    if(f_solvent_in_shell < -1.0e-24){
        printf("cs_ellip Q %g",q);
        printf(" ReqC %g",radius_equat_core);
        printf(" Tshell = %g",thick_shell);
        printf(" SLDsh %g\n",shell_sld);}

    // from here on is the same as the original core_shell_ellipsoid
    const double polar_core = radius_equat_core*x_core;
    const double equat_shell = radius_equat_core + thick_shell;
    const double polar_shell = radius_equat_core*x_core + thick_shell*x_polar_shell;

    // translate from [-1, 1] => [0, 1]
    const double m = 0.5;
    const double b = 0.5;
    double total_1 = 0.0;     //initialize intergral
    double total_2 = 0.0;     //initialize intergral
    for(int i=0;i<GAUSS_N;i++) {
        const double cos_theta = GAUSS_Z[i]*m + b;
        const double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
        double fq = _cs_ellipsoid_kernel(q*sin_theta, q*cos_theta,
            radius_equat_core, polar_core,
            equat_shell, polar_shell,
            sld_core_shell, sld_shell_solvent);
        total_1 += GAUSS_W[i] * fq;
        total_2 += GAUSS_W[i] * fq * fq;
    }

    // convert to [cm-1]   TODO check rescale for F1
    *F1 = 1.0e-2*sqrt(m)*total_1;
    *F2 = 1.0e-4*m*total_2;
}

static double
Iqac(double qab, double qc,
    double radius_equat_core,
    double x_core,
    double vol_dry_shell_over_core,
    double x_polar_shell,
    double core_sld,
    double dry_shell_sld,
    double solvent_sld,
    double f_solvent_in_shell)
{
    // sort out the parameterized sld
    const double f_solvent = fabs(f_solvent_in_shell);
    const double shell_sld = f_solvent*solvent_sld + (1-f_solvent)*dry_shell_sld;
    const double sld_core_shell = core_sld - shell_sld;
    const double sld_shell_solvent = shell_sld - solvent_sld;
    // Now we have to set up and solve the one real root of a cubic equation for thick_shell
    const double thick_shell = solve_shell_v3(radius_equat_core, x_core, vol_dry_shell_over_core,
    x_polar_shell, f_solvent);
    // from here on is the same as the original core_shell_ellipsoid
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
