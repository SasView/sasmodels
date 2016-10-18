double form_volume(double radius_equat_core,
                   double polar_core,
                   double equat_shell,
                   double polar_shell);
double Iq(double q,
          double radius_equat_core,
          double x_core,
          double thick_shell,
          double x_polar_shell,
          double core_sld,
          double shell_sld,
          double solvent_sld);


double Iqxy(double qx, double qy,
          double radius_equat_core,
          double x_core,
          double thick_shell,
          double x_polar_shell,
          double core_sld,
          double shell_sld,
          double solvent_sld,
          double theta,
          double phi);


double form_volume(double radius_equat_core,
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
core_shell_ellipsoid_xt_kernel(double q,
          double radius_equat_core,
          double x_core,
          double thick_shell,
          double x_polar_shell,
          double core_sld,
          double shell_sld,
          double solvent_sld)
{
    const double lolim = 0.0;
    const double uplim = 1.0;


    const double delpc = core_sld - shell_sld; //core - shell
    const double delps = shell_sld - solvent_sld; //shell - solvent


    const double polar_core = radius_equat_core*x_core;
    const double equat_shell = radius_equat_core + thick_shell;
    const double polar_shell = radius_equat_core*x_core + thick_shell*x_polar_shell;

    double summ = 0.0;	 //initialize intergral
    for(int i=0;i<76;i++) {
        double zi = 0.5*( Gauss76Z[i]*(uplim-lolim) + uplim + lolim );
        double yyy = gfn4(zi, radius_equat_core, polar_core, equat_shell,
                          polar_shell, delpc, delps, q);
        summ += Gauss76Wt[i] * yyy;
    }
    summ *= 0.5*(uplim-lolim);

    // convert to [cm-1]
    return 1.0e-4 * summ;
}

static double
core_shell_ellipsoid_xt_kernel_2d(double qx, double qy,
          double radius_equat_core,
          double x_core,
          double thick_shell,
          double x_polar_shell,
          double core_sld,
          double shell_sld,
          double solvent_sld,
          double theta,
          double phi)
{
    double q, sin_alpha, cos_alpha;
    ORIENT_SYMMETRIC(qx, qy, theta, phi, q, sin_alpha, cos_alpha);

    const double sldcs = core_sld - shell_sld;
    const double sldss = shell_sld- solvent_sld;

    const double polar_core = radius_equat_core*x_core;
    const double equat_shell = radius_equat_core + thick_shell;
    const double polar_shell = radius_equat_core*x_core + thick_shell*x_polar_shell;

    // Call the IGOR library function to get the kernel:
    // MUST use gfn4 not gf2 because of the def of params.
    double answer = gfn4(cos_alpha,
                  radius_equat_core,
                  polar_core,
                  equat_shell,
                  polar_shell,
                  sldcs,
                  sldss,
                  q);

    //convert to [cm-1]
    answer *= 1.0e-4;

    return answer;
}

double Iq(double q,
          double radius_equat_core,
          double x_core,
          double thick_shell,
          double x_polar_shell,
          double core_sld,
          double shell_sld,
          double solvent_sld)
{
    double intensity = core_shell_ellipsoid_xt_kernel(q,
           radius_equat_core,
           x_core,
           thick_shell,
           x_polar_shell,
           core_sld,
           shell_sld,
           solvent_sld);

    return intensity;
}


double Iqxy(double qx, double qy,
          double radius_equat_core,
          double x_core,
          double thick_shell,
          double x_polar_shell,
          double core_sld,
          double shell_sld,
          double solvent_sld,
          double theta,
          double phi)
{
    double intensity = core_shell_ellipsoid_xt_kernel_2d(qx, qy,
                       radius_equat_core,
                       x_core,
                       thick_shell,
                       x_polar_shell,
                       core_sld,
                       shell_sld,
                       solvent_sld,
                       theta,
                       phi);

    return intensity;
}
