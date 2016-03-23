double form_volume(double core_thick,
                   double layer_thick,
                   double radius,
                   double n_stacking);

double Iq(double q,
          double core_thick,
          double layer_thick,
          double radius,
          double n_stacking,
          double sigma_d,
          double core_sld,
          double layer_sld,
          double solvent_sld);

double Iqxy(double qx, double qy,
            double core_thick,
            double layer_thick,
            double radius,
            double n_stacking,
            double sigma_d,
            double core_sld,
            double layer_sld,
            double solvent_sld,
            double theta,
            double phi);

static
double _kernel(double qq,
               double radius,
               double core_sld,
               double layer_sld,
               double solvent_sld,
               double halfheight,
               double layer_thick,
               double zi,
               double sigma_d,
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
	const double sinarg2 = qq*(halfheight+layer_thick)*cos(zi);

    const double be1 = sas_J1c(besarg1);
	const double be2 = sas_J1c(besarg2);
	const double si1 = sin(sinarg1)/sinarg1;
	const double si2 = sin(sinarg2)/sinarg2;

	const double dr1 = (core_sld-solvent_sld);
	const double dr2 = (layer_sld-solvent_sld);
	const double area = M_PI*radius*radius;
	const double totald=2.0*(layer_thick+halfheight);

	const double t1 = area*(2.0*halfheight)*dr1*(si1)*(be1);
	const double t2 = area*dr2*(totald*si2-2.0*halfheight*si1)*(be2);


	double retval =((t1+t2)*(t1+t2))*sin(zi);

	// loop for the structure facture S(q)
	double sqq=0.0;
	for(int kk=1;kk<n_stacking;kk+=1) {
		double dexpt=qq*cos(zi)*qq*cos(zi)*d*d*sigma_d*sigma_d*kk/2.0;
		sqq=sqq+(n_stacking-kk)*cos(qq*cos(zi)*d*kk)*exp(-1.*dexpt);
	}

	// end of loop for S(q)
	sqq=1.0+2.0*sqq/n_stacking;

	retval *= sqq;

	return(retval);
}


static
double stacked_disks_kernel(double q,
                            double core_thick,
                            double layer_thick,
                            double radius,
                            double n_stacking,
                            double sigma_d,
                            double core_sld,
                            double layer_sld,
                            double solvent_sld)
{
/*	StackedDiscsX  :  calculates the form factor of a stacked "tactoid" of core shell disks
like clay platelets that are not exfoliated
*/
	double summ = 0.0;	//initialize integral

	double d=2.0*layer_thick+core_thick;
	double halfheight = core_thick/2.0;

	for(int i=0;i<N_POINTS_76;i++) {
		double zi = (Gauss76Z[i] + 1.0)*M_PI/4.0;
		double yyy = Gauss76Wt[i] *
                    _kernel(q,
		                   radius,
		                   core_sld,
		                   layer_sld,
		                   solvent_sld,
		                   halfheight,
		                   layer_thick,
		                   zi,
		                   sigma_d,
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
                            double core_thick,
                            double layer_thick,
                            double radius,
                            double n_stacking,
                            double sigma_d,
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
    const double cyl_x = ct * cp;
    const double cyl_y = st;

    // Compute the angle btw vector q and the
    // axis of the parallelepiped
    const double cos_val = cyl_x*q_x + cyl_y*q_y;

    // Note: cos(alpha) = 0 and 1 will get an
    // undefined value from Stackdisc_kern
    double alpha = acos( cos_val );

    // Call the IGOR library function to get the kernel
    double d = 2 * layer_thick + core_thick;
    double halfheight = core_thick/2.0;
    double answer = _kernel(q,
                     radius,
                     core_sld,
                     layer_sld,
                     solvent_sld,
                     halfheight,
                     layer_thick,
                     alpha,
                     sigma_d,
                     d,
                     n_stacking);

    answer /= sin(alpha);
    //convert to [cm-1]
    answer *= 1.0e-4;

    return answer;
}

double form_volume(double core_thick,
                   double layer_thick,
                   double radius,
                   double n_stacking){
    double d = 2 * layer_thick + core_thick;
    return acos(-1.0) * radius * radius * d * n_stacking;
}

double Iq(double q,
          double core_thick,
          double layer_thick,
          double radius,
          double n_stacking,
          double sigma_d,
          double core_sld,
          double layer_sld,
          double solvent_sld)
{
    return stacked_disks_kernel(q,
                    core_thick,
                    layer_thick,
                    radius,
                    n_stacking,
                    sigma_d,
                    core_sld,
                    layer_sld,
                    solvent_sld);
}

// Iqxy is never called since no orientation or magnetic parameters.
double Iqxy(double qx, double qy,
            double core_thick,
            double layer_thick,
            double radius,
            double n_stacking,
            double sigma_d,
            double core_sld,
            double layer_sld,
            double solvent_sld,
            double theta,
            double phi)
{
    double q = sqrt(qx*qx + qy*qy);
    return stacked_disks_kernel_2d(q, qx/q, qy/q,
                    core_thick,
                    layer_thick,
                    radius,
                    n_stacking,
                    sigma_d,
                    core_sld,
                    layer_sld,
                    solvent_sld,
                    theta,
                    phi);
}

