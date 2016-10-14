double form_volume(double radius);
double Iq(double q,double dnn,double d_factor, double radius,double sld, double solvent_sld);
double Iqxy(double qx, double qy, double dnn,
    double d_factor, double radius,double sld, double solvent_sld,
    double theta, double phi, double psi);

double _FCC_Integrand(double q, double dnn, double d_factor, double theta, double phi);
double _FCCeval(double Theta, double Phi, double temp1, double temp3);


double _FCC_Integrand(double q, double dnn, double d_factor, double theta, double phi) {

	const double Da = d_factor*dnn;
	const double temp1 = q*q*Da*Da;
	const double temp3 = q*dnn;

	double retVal = _FCCeval(theta,phi,temp1,temp3)/(4.0*M_PI);
	return(retVal);
}

double _FCCeval(double Theta, double Phi, double temp1, double temp3) {

	double result;
        double sin_theta,cos_theta,sin_phi,cos_phi;
	SINCOS(Theta, sin_theta, cos_theta);
	SINCOS(Phi, sin_phi, cos_phi);

	const double temp6 =  sin_theta;
	const double temp7 =  sin_theta*sin_phi + cos_theta;
	const double temp8 = -sin_theta*cos_phi + cos_theta;
	const double temp9 = -sin_theta*cos_phi + sin_theta*sin_phi;

	const double temp10 = exp((-1.0/8.0)*temp1*((temp7*temp7)+(temp8*temp8)+(temp9*temp9)));
	result = pow((1.0-(temp10*temp10)),3)*temp6
	    / ( (1.0 - 2.0*temp10*cos(0.5*temp3*temp7) + temp10*temp10)
	      * (1.0 - 2.0*temp10*cos(0.5*temp3*temp8) + temp10*temp10)
	      * (1.0 - 2.0*temp10*cos(0.5*temp3*temp9) + temp10*temp10));

	return (result);
}

double form_volume(double radius){
    return sphere_volume(radius);
}


double Iq(double q, double dnn,
  double d_factor, double radius,
  double sld, double solvent_sld){

	//Volume fraction calculated from lattice symmetry and sphere radius
	const double s1 = dnn*sqrt(2.0);
	const double latticescale = 4.0*sphere_volume(radius/s1);

    const double va = 0.0;
    const double vb = 2.0*M_PI;
    const double vaj = 0.0;
    const double vbj = M_PI;

    double summ = 0.0;
    double answer = 0.0;
	for(int i=0; i<150; i++) {
		//setup inner integral over the ellipsoidal cross-section
		double summj=0.0;
		const double zphi = ( Gauss150Z[i]*(vb-va) + va + vb )/2.0;		//the outer dummy is phi
		for(int j=0;j<150;j++) {
			//20 gauss points for the inner integral
			double ztheta = ( Gauss150Z[j]*(vbj-vaj) + vaj + vbj )/2.0;		//the inner dummy is theta
			double yyy = Gauss150Wt[j] * _FCC_Integrand(q,dnn,d_factor,ztheta,zphi);
			summj += yyy;
		}
		//now calculate the value of the inner integral
		double answer = (vbj-vaj)/2.0*summj;

		//now calculate outer integral
		summ = summ+(Gauss150Wt[i] * answer);
	}		//final scaling is done at the end of the function, after the NT_FP64 case

	answer = (vb-va)/2.0*summ;
	answer = answer*sphere_form(q,radius,sld,solvent_sld)*latticescale;

    return answer;


}

double Iqxy(double qx, double qy,
    double dnn, double d_factor, double radius,
    double sld, double solvent_sld,
    double theta, double phi, double psi)
{
    double q, cos_a1, cos_a2, cos_a3;
    ORIENT_ASYMMETRIC(qx, qy, theta, phi, psi, q, cos_a3, cos_a2, cos_a1);

    const double a1 = cos_a2 + cos_a3;
    const double a2 = cos_a3 + cos_a1;
    const double a3 = cos_a2 + cos_a1;
    const double qd = 0.5*q*dnn;
    const double arg = 0.5*square(qd*d_factor)*(a1*a1 + a2*a2 + a3*a3);
    const double tanh_qd = tanh(arg);
    const double cosh_qd = cosh(arg);
    const double Zq = tanh_qd/(1. - cos(qd*a1)/cosh_qd)
                    * tanh_qd/(1. - cos(qd*a2)/cosh_qd)
                    * tanh_qd/(1. - cos(qd*a3)/cosh_qd);

    //if (isnan(Zq)) printf("q:(%g,%g) qd: %g a1: %g a2: %g a3: %g arg: %g\n", qx, qy, qd, a1, a2, a3, arg);

    const double Fq = sphere_form(q,radius,sld,solvent_sld)*Zq;
    //the occupied volume of the lattice
    const double lattice_scale = 4.0*sphere_volume(M_SQRT1_2*radius/dnn);
    return lattice_scale * Fq;
}
