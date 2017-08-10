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
    return sphere_volume(radius);
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
	double temp2 = cube(-expm1(-temp1));
	double temp3 = qq*dnn;
	double temp4 = 2.0*exp(-0.5*temp1);
	double temp5 = exp(-1.0*temp1);

	double integrand = temp2*sc_eval(yy,xx,temp3,temp4,temp5)/M_PI_2;

	return(integrand);
}

double Iq(double q,
          double dnn,
          double d_factor,
          double radius,
          double sphere_sld,
          double solvent_sld)
{
	const double va = 0.0;
	const double vb = M_PI_2; //orientation average, outer integral

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
	// NB: 4/3 pi r^3 / dnn^3 = 4/3 pi(r/dnn)^3
	const double latticeScale = sphere_volume(radius/dnn);

	answer *= sphere_form(q, radius, sphere_sld, solvent_sld)*latticeScale;

	return answer;
}

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
    double q, zhat, yhat, xhat;
    ORIENT_ASYMMETRIC(qx, qy, theta, phi, psi, q, xhat, yhat, zhat);

    const double qd = q*dnn;
    const double arg = 0.5*square(qd*d_factor);
    const double tanh_qd = tanh(arg);
    const double cosh_qd = cosh(arg);
    const double Zq = tanh_qd/(1. - cos(qd*zhat)/cosh_qd)
                    * tanh_qd/(1. - cos(qd*yhat)/cosh_qd)
                    * tanh_qd/(1. - cos(qd*xhat)/cosh_qd);

    const double Fq = sphere_form(q, radius, sphere_sld, solvent_sld)*Zq;
    //the occupied volume of the lattice
    const double lattice_scale = sphere_volume(radius/dnn);
    return lattice_scale * Fq;
}
