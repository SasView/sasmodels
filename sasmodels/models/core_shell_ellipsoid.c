double form_volume(double radius_equat_core,
                   double radius_polar_core,
                   double radius_equat_shell,
                   double radius_polar_shell);
double Iq(double q,
          double radius_equat_core,
          double radius_polar_core,
          double radius_equat_shell,
          double radius_polar_shell,
          double sld_core,
          double sld_shell,
          double sld_solvent);


double Iqxy(double qx, double qy,
          double radius_equat_core,
          double radius_polar_core,
          double radius_equat_shell,
          double radius_polar_shell,
          double sld_core,
          double sld_shell,
          double sld_solvent,
          double theta,
          double phi);


double form_volume(double radius_equat_core,
                   double radius_polar_core,
                   double radius_equat_shell,
                   double radius_polar_shell)
{
    double vol = 4.0*M_PI/3.0*radius_equat_shell*radius_equat_shell*radius_polar_shell;
    return vol;
}

static double
core_shell_ellipsoid_kernel(double q,
          double radius_equat_core,
          double radius_polar_core,
          double radius_equat_shell,
          double radius_polar_shell,
          double sld_core,
          double sld_shell,
          double sld_solvent)
{

    //upper and lower integration limits
    const double lolim = 0.0;
    const double uplim = 1.0;

    double summ = 0.0;	 //initialize intergral

    const double delpc = sld_core - sld_shell;    //core - shell
    const double delps = sld_shell - sld_solvent; //shell - solvent

    for(int i=0;i<N_POINTS_76;i++) {
        double zi = ( Gauss76Z[i]*(uplim-lolim) + uplim + lolim )/2.0;
        double yyy = Gauss76Wt[i] * gfn4(zi,
                                  radius_equat_core,
                                  radius_polar_core,
                                  radius_equat_shell,
                                  radius_polar_shell,
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
core_shell_ellipsoid_kernel_2d(double q, double q_x, double q_y,
          double radius_equat_core,
          double radius_polar_core,
          double radius_equat_shell,
          double radius_polar_shell,
          double sld_core,
          double sld_shell,
          double sld_solvent,
          double theta,
          double phi)
{
    //convert angle degree to radian
    theta = theta * M_PI_180;
    phi = phi * M_PI_180;


    // ellipsoid orientation, the axis of the rotation is consistent with the ploar axis.
    const double cyl_x = cos(theta) * cos(phi);
    const double cyl_y = sin(theta);

    const double sldcs = sld_core - sld_shell;
    const double sldss = sld_shell- sld_solvent;

    // Compute the angle btw vector q and the
    // axis of the cylinder
    const double cos_val = cyl_x*q_x + cyl_y*q_y;

    // Call the IGOR library function to get the kernel: MUST use gfn4 not gf2 because of the def of params.
    double answer = gfn4(cos_val,
                  radius_equat_core,
                  radius_polar_core,
                  radius_equat_shell,
                  radius_polar_shell,
                  sldcs,
                  sldss,
                  q);

    //convert to [cm-1]
    answer *= 1.0e-4;

    return answer;
}

double Iq(double q,
          double radius_equat_core,
          double radius_polar_core,
          double radius_equat_shell,
          double radius_polar_shell,
          double sld_core,
          double sld_shell,
          double sld_solvent)
{
    double intensity = core_shell_ellipsoid_kernel(q,
           radius_equat_core,
           radius_polar_core,
           radius_equat_shell,
           radius_polar_shell,
           sld_core,
           sld_shell,
           sld_solvent);

    return intensity;
}


double Iqxy(double qx, double qy,
          double radius_equat_core,
          double radius_polar_core,
          double radius_equat_shell,
          double radius_polar_shell,
          double sld_core,
          double sld_shell,
          double sld_solvent,
          double theta,
          double phi)
{
    double q;
    q = sqrt(qx*qx+qy*qy);
    double intensity = core_shell_ellipsoid_kernel_2d(q, qx/q, qy/q,
                       radius_equat_core,
                       radius_polar_core,
                       radius_equat_shell,
                       radius_polar_shell,
                       sld_core,
                       sld_shell,
                       sld_solvent,
                       theta,
                       phi);

    return intensity;
}
