double form_volume(double equat_core,
                   double polar_core,
                   double equat_shell,
                   double polar_shell);
double Iq(double q,
          double equat_core,
          double x_core,
          double t_shell,
          double x_polar_shell,
          double core_sld,
          double shell_sld,
          double solvent_sld);


double Iqxy(double qx, double qy,
          double equat_core,
          double x_core,
          double t_shell,
          double x_polar_shell,
          double core_sld,
          double shell_sld,
          double solvent_sld,
          double theta,
          double phi);


double form_volume(double equat_core,
                   double x_core,
                   double t_shell,
                   double x_polar_shell)
{
	double equat_shell, polar_shell;
    equat_shell = equat_core + t_shell;
    polar_shell = equat_core*x_core + t_shell*x_polar_shell;
    double vol = 4.0*M_PI/3.0*equat_shell*equat_shell*polar_shell;
    return vol;
}

static double
core_shell_ellipsoid_xt_kernel(double q,
          double equat_core,
          double x_core,
          double t_shell,
          double x_polar_shell,
          double core_sld,
          double shell_sld,
          double solvent_sld)
{
	double delpc,delps;
	double uplim,lolim;		//upper and lower integration limits
	double summ,zi,yyy,answer; //running tally of integration
	double polar_core, equat_shell, polar_shell;

	lolim = 0.0;
	uplim = 1.0;

	summ = 0.0;	 //initialize intergral

	delpc = core_sld - shell_sld; //core - shell
	delps = shell_sld - solvent_sld; //shell - solvent


    polar_core = equat_core*x_core;
    equat_shell = equat_core + t_shell;
    polar_shell = equat_core*x_core + t_shell*x_polar_shell;

	for(int i=0;i<76;i++) {
		zi = ( Gauss76Z[i]*(uplim-lolim) + uplim + lolim )/2.0;
		yyy = Gauss76Wt[i] * gfn4(zi,
		                          equat_core,
                                  polar_core,
                                  equat_shell,
                                  polar_shell,
		                          delpc,
		                          delps,
		                          q);
		summ += yyy;
	}

	answer = (uplim-lolim)/2.0*summ;

	//convert to [cm-1]
	answer *= 1.0e-4;

	return answer;
}

static double
core_shell_ellipsoid_xt_kernel_2d(double q, double q_x, double q_y,
          double equat_core,
          double x_core,
          double t_shell,
          double x_polar_shell,
          double core_sld,
          double shell_sld,
          double solvent_sld,
          double theta,
          double phi)
{
    double cyl_x, cyl_y;
    double cos_val;
    double answer;
    double sldcs,sldss;
	double polar_core, equat_shell, polar_shell;

    //convert angle degree to radian
    theta = theta * M_PI/180.0;
    phi = phi * M_PI/180.0;


    // ellipsoid orientation, the axis of the rotation is consistent with the ploar axis.
    cyl_x = cos(theta) * cos(phi);
    cyl_y = sin(theta);

    sldcs = core_sld - shell_sld;
    sldss = shell_sld- solvent_sld;

    // Compute the angle btw vector q and the
    // axis of the cylinder
    cos_val = cyl_x*q_x + cyl_y*q_y;

    polar_core = equat_core*x_core;
    equat_shell = equat_core + t_shell;
    polar_shell = equat_core*x_core + t_shell*x_polar_shell;

    // Call the IGOR library function to get the kernel:
    // MUST use gfn4 not gf2 because of the def of params.
    answer = gfn4(cos_val,
                  equat_core,
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
          double equat_core,
          double x_core,
          double t_shell,
          double x_polar_shell,
          double core_sld,
          double shell_sld,
          double solvent_sld)
{
    double intensity = core_shell_ellipsoid_xt_kernel(q,
           equat_core,
           x_core,
           t_shell,
           x_polar_shell,
           core_sld,
           shell_sld,
           solvent_sld);

    return intensity;
}


double Iqxy(double qx, double qy,
          double equat_core,
          double x_core,
          double t_shell,
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
                       equat_core,
                       x_core,
                       t_shell,
                       x_polar_shell,
                       core_sld,
                       shell_sld,
                       solvent_sld,
                       theta,
                       phi);

    return intensity;
}
