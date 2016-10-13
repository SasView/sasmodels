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
    double vol = 4.0*M_PI/3.0*equat_shell*equat_shell*polar_shell;
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

    double summ = 0.0;	 //initialize intergral

    const double delpc = core_sld - shell_sld; //core - shell
    const double delps = shell_sld - solvent_sld; //shell - solvent


    const double polar_core = radius_equat_core*x_core;
    const double equat_shell = radius_equat_core + thick_shell;
    const double polar_shell = radius_equat_core*x_core + thick_shell*x_polar_shell;

    for(int i=0;i<N_POINTS_76;i++) {
        double zi = ( Gauss76Z[i]*(uplim-lolim) + uplim + lolim )/2.0;
        double yyy = Gauss76Wt[i] * gfn4(zi,
                                  radius_equat_core,
                                  polar_core,
                                  equat_shell,
                                  polar_shell,
                                  delpc,
                                  delps,
                                  q);
        summ += yyy;
    }

    double answer = (uplim-lolim)/2.0*summ;
    //convert to [cm-1]
    answer *= 1.0e-4;

    return answer;
}

static double
core_shell_ellipsoid_xt_kernel_2d(double q, double q_x, double q_y,
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
    //convert angle degree to radian
    theta = theta * M_PI_180;
    phi = phi * M_PI_180;

    // ellipsoid orientation, the axis of the rotation is consistent with the ploar axis.
    const double cyl_x = sin(theta) * cos(phi);
    const double cyl_y = sin(theta) * sin(phi);

    const double sldcs = core_sld - shell_sld;
    const double sldss = shell_sld- solvent_sld;

    // Compute the angle btw vector q and the
    // axis of the cylinder
    const double cos_val = cyl_x*q_x + cyl_y*q_y;

    const double polar_core = radius_equat_core*x_core;
    const double equat_shell = radius_equat_core + thick_shell;
    const double polar_shell = radius_equat_core*x_core + thick_shell*x_polar_shell;

    // Call the IGOR library function to get the kernel:
    // MUST use gfn4 not gf2 because of the def of params.
    double answer = gfn4(cos_val,
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
    double q;
    q = sqrt(qx*qx+qy*qy);
    double intensity = core_shell_ellipsoid_xt_kernel_2d(q, qx/q, qy/q,
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
