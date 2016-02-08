double form_volume(double radius);

double Iq(double q,
          double dnn,
          double d_factor,
          double radius,
          double sphere_sld,
          double solvent_sld);

double Iqxy(double qx, double qy,
            double dnn,
            double d_factor,
            double radius,
            double sphere_sld,
            double solvent_sld,
            double theta,
            double phi,
            double psi);

double form_volume(double radius)
{
    return 1.333333333333333*M_PI*radius*radius*radius;
}

static double
sc_eval(double theta, double phi, double temp3, double temp4, double temp5)
{
    double cnt, snt;
    SINCOS(theta, cnt, snt);

    double cnp, snp;
    SINCOS(phi, cnp, snp);

	double temp6 = snt;
	double temp7 = -1.0*temp3*snt*cnp;
	double temp8 = temp3*snt*snp;
	double temp9 = temp3*cnt;
	double result = temp6/((1.0-temp4*cos((temp7))+temp5)*
	                       (1.0-temp4*cos((temp8))+temp5)*
	                       (1.0-temp4*cos((temp9))+temp5));
	return (result);
}

static double
sc_integrand(double dnn, double d_factor, double qq, double xx, double yy)
{
    //Function to calculate integrand values for simple cubic structure

	double da = d_factor*dnn;
	double temp1 = qq*qq*da*da;
	double temp2 = pow( 1.0-exp(-1.0*temp1) ,3);
	double temp3 = qq*dnn;
	double temp4 = 2.0*exp(-0.5*temp1);
	double temp5 = exp(-1.0*temp1);

	double integrand = temp2*sc_eval(yy,xx,temp3,temp4,temp5);
	integrand *= 2.0/M_PI;

	return(integrand);
}

static
double sc_crystal_kernel(double q,
          double dnn,
          double d_factor,
          double radius,
          double sphere_sld,
          double solvent_sld)
{
	const double va = 0.0;
	const double vb = M_PI/2.0; //orientation average, outer integral

    double summ=0.0;
    double answer=0.0;

	for(int i=0;i<150;i++) {
		//setup inner integral over the ellipsoidal cross-section
		double summj=0.0;
		double zi = ( Gauss150Z[i]*(vb-va) + va + vb )/2.0;
		for(int j=0;j<150;j++) {
			//150 gauss points for the inner integral
			double zij = ( Gauss150Z[j]*(vb-va) + va + vb )/2.0;
			double tmp = Gauss150Wt[j] * sc_integrand(dnn,d_factor,q,zi,zij);
			summj += tmp;
		}
		//now calculate the value of the inner integral
		answer = (vb-va)/2.0*summj;

		//now calculate outer integral
		double tmp = Gauss150Wt[i] * answer;
		summ += tmp;
	}		//final scaling is done at the end of the function, after the NT_FP64 case

	answer = (vb-va)/2.0*summ;

	//Volume fraction calculated from lattice symmetry and sphere radius
	const double latticeScale = (4.0/3.0)*M_PI*(radius*radius*radius)/pow(dnn,3);

	//answer *= sphere_form_paracrystal(q, radius,contrast)*latticeScale;
	answer *= sphere_form(q, radius, sphere_sld, solvent_sld)*latticeScale;

	return answer;
}

static
double sc_crystal_kernel_2d(double q, double q_x, double q_y,
          double dnn,
          double d_factor,
          double radius,
          double sphere_sld,
          double solvent_sld,
          double theta,
          double phi,
          double psi)
{
    //convert angle degree to radian
    theta = theta * M_PI_180;
    phi = phi * M_PI_180;
    psi = psi * M_PI_180;

    const double qda_2 = pow(q*d_factor*dnn,2.0);

    double snt, cnt;
    SINCOS(theta, snt, cnt);

    double snp, cnp;
    SINCOS(phi, snp, cnp);

    double sns, cns;
    SINCOS(psi, sns, cns);

    /// Angles here are respect to detector coordinate instead of against
    //  q coordinate(PRB 36, 3, 1754)
    // a3 axis orientation

    const double a3_x = cnt * cnp;
    const double a3_y = snt;

    // Compute the angle btw vector q and the a3 axis
    double cos_val_a3 = a3_x*q_x + a3_y*q_y;

    // a1 axis orientation
    const double a1_x = -cnp*sns * snt+snp*cns;
    const double a1_y = sns*cnt;

    double cos_val_a1 = a1_x*q_x + a1_y*q_y;

    // a2 axis orientation
    const double a2_x = -snt*cns*cnp-sns*snp;
    const double a2_y = cnt*cns;

    // a2 axis
    const double cos_val_a2 =  a2_x*q_x + a2_y*q_y;

    // The following test should always pass
    if (fabs(cos_val_a3)>1.0) {
        //printf("parallel_ana_2D: Unexpected error: cos(alpha)>1\n");
        cos_val_a3 = 1.0;
    }
    if (fabs(cos_val_a1)>1.0) {
        //printf("parallel_ana_2D: Unexpected error: cos(alpha)>1\n");
        cos_val_a1 = 1.0;
    }
    if (fabs(cos_val_a2)>1.0) {
        //printf("parallel_ana_2D: Unexpected error: cos(alpha)>1\n");
        cos_val_a3 = 1.0;
    }

    const double a3_dot_q = dnn*q*cos_val_a3;
    const double a1_dot_q = dnn*q*cos_val_a1;
    const double a2_dot_q = dnn*q*cos_val_a2;

    // Call Zq=Z1*Z2*Z3
    double Zq = (1.0-exp(-qda_2))/(1.0-2.0*exp(-0.5*qda_2)*cos(a1_dot_q)+exp(-qda_2));
    Zq *= (1.0-exp(-qda_2))/(1.0-2.0*exp(-0.5*qda_2)*cos(a2_dot_q)+exp(-qda_2));
    Zq *= (1.0-exp(-qda_2))/(1.0-2.0*exp(-0.5*qda_2)*cos(a3_dot_q)+exp(-qda_2));

    // Use SphereForm directly from libigor
    double answer = sphere_form(q, radius, sphere_sld, solvent_sld)*Zq;

    //consider scales
    const double latticeScale = (4.0/3.0)*M_PI*(radius*radius*radius)/pow(dnn,3.0);
    answer *= latticeScale;

    return answer;
}

double Iq(double q,
          double dnn,
          double d_factor,
          double radius,
          double sphere_sld,
          double solvent_sld)
{
    return sc_crystal_kernel(q,
              dnn,
              d_factor,
              radius,
              sphere_sld,
              solvent_sld);
}

// Iqxy is never called since no orientation or magnetic parameters.
double Iqxy(double qx, double qy,
            double dnn,
            double d_factor,
            double radius,
            double sphere_sld,
            double solvent_sld,
            double theta,
            double phi,
            double psi)
{
    double q = sqrt(qx*qx + qy*qy);


    return sc_crystal_kernel_2d(q, qx/q, qy/q,
                  dnn,
                  d_factor,
                  radius,
                  sphere_sld,
                  solvent_sld,
                  theta,
                  phi,
                  psi);

}

