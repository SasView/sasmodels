// some functions from core_shell_ellipsoid, also used by core_shell_ellipsoid_tied
// RKH 19/09/2023
// Converted from Igor function gfn4, using the same pattern as ellipsoid
// for evaluating the parts of the integral.
//     FUNCTION gfn4:    CONTAINS F(Q,A,B,MU)**2  AS GIVEN
//                       BY (53) & (58-59) IN CHEN AND
//                       KOTLARCHYK REFERENCE
//

static double
_cs_ellipsoid_kernel(double qab, double qc,
    double equat_core, double polar_core,
    double equat_shell, double polar_shell,
    double sld_core_shell, double sld_shell_solvent)
{
    const double qr_core = sqrt(square(equat_core*qab) + square(polar_core*qc));
    const double si_core = sas_3j1x_x(qr_core);
    const double volume_core = M_4PI_3*equat_core*equat_core*polar_core;
    const double fq_core = si_core*volume_core*sld_core_shell;

    const double qr_shell = sqrt(square(equat_shell*qab) + square(polar_shell*qc));
    const double si_shell = sas_3j1x_x(qr_shell);
    const double volume_shell = M_4PI_3*equat_shell*equat_shell*polar_shell;
    const double fq_shell = si_shell*volume_shell*sld_shell_solvent;

    return fq_core + fq_shell;
}


static double
radius_from_curvature(double radius_equat_core, double x_core, double thick_shell, double x_polar_shell)
{
    // Trivial cases
    if (1.0 == x_core && 1.0 == x_polar_shell) return radius_equat_core + thick_shell;
    if ((radius_equat_core + thick_shell)*(radius_equat_core*x_core + thick_shell*x_polar_shell) == 0.)  return 0.;

    // see equation (26) in A.Isihara, J.Chem.Phys. 18(1950)1446-1449
    const double radius_equat_tot = radius_equat_core + thick_shell;
    const double radius_polar_tot = radius_equat_core*x_core + thick_shell*x_polar_shell;
    const double ratio = (radius_polar_tot < radius_equat_tot
                          ? radius_polar_tot / radius_equat_tot
                          : radius_equat_tot / radius_polar_tot);
    const double e1 = sqrt(1.0 - ratio*ratio);
    const double b1 = 1.0 + asin(e1) / (e1 * ratio);
    const double bL = (1.0 + e1) / (1.0 - e1);
    const double b2 = 1.0 + 0.5 * ratio * ratio / e1 * log(bL);
    const double delta = 0.75 * b1 * b2;
    const double ddd = 2.0 * (delta + 1.0) * radius_polar_tot * radius_equat_tot * radius_equat_tot;
    return 0.5 * cbrt(ddd);
}

