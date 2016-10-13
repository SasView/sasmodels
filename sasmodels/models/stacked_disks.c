double form_volume(double thick_core,
                   double thick_layer,
                   double radius,
                   double n_stacking);

double Iq(double q,
          double thick_core,
          double thick_layer,
          double radius,
          double n_stacking,
          double sigma_dnn,
          double core_sld,
          double layer_sld,
          double solvent_sld);

static
double _kernel(double qq,
               double radius,
               double core_sld,
               double layer_sld,
               double solvent_sld,
               double halfheight,
               double thick_layer,
               double zi,
               double sigma_dnn,
               double d,
               double n_stacking)

{
	// qq is the q-value for the calculation (1/A)
	// radius is the core radius of the cylinder (A)
	// *_sld are the respective SLD's
	// halfheight is the *Half* CORE-LENGTH of the cylinder = L (A)
	// zi is the dummy variable for the integration (x in Feigin's notation)

	const double besarg1 = qq*radius*sin(zi);
	const double besarg2 = qq*radius*sin(zi);

	const double sinarg1 = qq*halfheight*cos(zi);
	const double sinarg2 = qq*(halfheight+thick_layer)*cos(zi);

	const double be1 = sas_J1c(besarg1);
	const double be2 = sas_J1c(besarg2);
	const double si1 = sin(sinarg1)/sinarg1;
	const double si2 = sin(sinarg2)/sinarg2;

	const double dr1 = (core_sld-solvent_sld);
	const double dr2 = (layer_sld-solvent_sld);
	const double area = M_PI*radius*radius;
	const double totald=2.0*(thick_layer+halfheight);

	const double t1 = area*(2.0*halfheight)*dr1*(si1)*(be1);
	const double t2 = area*dr2*(totald*si2-2.0*halfheight*si1)*(be2);


	double retval =((t1+t2)*(t1+t2))*sin(zi);

	// loop for the structure facture S(q)
	double sqq=0.0;
	for(int kk=1;kk<n_stacking;kk+=1) {
		double dexpt=qq*cos(zi)*qq*cos(zi)*d*d*sigma_dnn*sigma_dnn*kk/2.0;
		sqq=sqq+(n_stacking-kk)*cos(qq*cos(zi)*d*kk)*exp(-1.*dexpt);
	}

	// end of loop for S(q)
	sqq=1.0+2.0*sqq/n_stacking;

	retval *= sqq;

	return(retval);
}


static
double stacked_disks_kernel(double q,
                            double thick_core,
                            double thick_layer,
                            double radius,
                            double n_stacking,
                            double sigma_dnn,
                            double core_sld,
                            double layer_sld,
                            double solvent_sld)
{
/*	StackedDiscsX  :  calculates the form factor of a stacked "tactoid" of core shell disks
like clay platelets that are not exfoliated
*/
	double summ = 0.0;	//initialize integral

	double d=2.0*thick_layer+thick_core;
	double halfheight = thick_core/2.0;

	for(int i=0;i<N_POINTS_76;i++) {
		double zi = (Gauss76Z[i] + 1.0)*M_PI/4.0;
		double yyy = Gauss76Wt[i] *
                    _kernel(q,
		                   radius,
		                   core_sld,
		                   layer_sld,
		                   solvent_sld,
		                   halfheight,
		                   thick_layer,
		                   zi,
		                   sigma_dnn,
		                   d,
		                   n_stacking);
		summ += yyy;
	}

	double answer = M_PI/4.0*summ;

	//Convert to [cm-1]
	answer *= 1.0e-4;

	return answer;
}

static double stacked_disks_kernel_2d(double q, double q_x, double q_y,
                            double thick_core,
                            double thick_layer,
                            double radius,
                            double n_stacking,
                            double sigma_dnn,
                            double core_sld,
                            double layer_sld,
                            double solvent_sld,
                            double theta,
                            double phi)
{

    double ct, st, cp, sp;

    //convert angle degree to radian
    theta = theta * M_PI/180.0;
    phi = phi * M_PI/180.0;

    SINCOS(theta, st, ct);
    SINCOS(phi, sp, cp);

    // silence compiler warnings about unused variable
    (void) sp;

    // parallelepiped orientation
    const double cyl_x = st * cp;
    const double cyl_y = st * sp;

    // Compute the angle btw vector q and the
    // axis of the parallelepiped
    const double cos_val = cyl_x*q_x + cyl_y*q_y;

    // Note: cos(alpha) = 0 and 1 will get an
    // undefined value from Stackdisc_kern
    double alpha = acos( cos_val );

    // Call the IGOR library function to get the kernel
    double d = 2 * thick_layer + thick_core;
    double halfheight = thick_core/2.0;
    double answer = _kernel(q,
                     radius,
                     core_sld,
                     layer_sld,
                     solvent_sld,
                     halfheight,
                     thick_layer,
                     alpha,
                     sigma_dnn,
                     d,
                     n_stacking);

    answer /= sin(alpha);
    //convert to [cm-1]
    answer *= 1.0e-4;

    return answer;
}

double form_volume(double thick_core,
                   double thick_layer,
                   double radius,
                   double n_stacking){
    double d = 2 * thick_layer + thick_core;
    return acos(-1.0) * radius * radius * d * n_stacking;
}

double Iq(double q,
          double thick_core,
          double thick_layer,
          double radius,
          double n_stacking,
          double sigma_dnn,
          double core_sld,
          double layer_sld,
          double solvent_sld)
{
    return stacked_disks_kernel(q,
                    thick_core,
                    thick_layer,
                    radius,
                    n_stacking,
                    sigma_dnn,
                    core_sld,
                    layer_sld,
                    solvent_sld);
}
