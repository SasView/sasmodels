double form_volume(double radius);
double Iq(double q,double dnn,double d_factor, double radius,double sld, double solvent_sld);
double Iqxy(double qx, double qy, double dnn,
    double d_factor, double radius,double sld, double solvent_sld,
    double theta, double phi, double psi);

double _FCC_Integrand(double q, double dnn, double d_factor, double theta, double phi);
double _FCCeval(double Theta, double Phi, double temp1, double temp3);
double _sphereform(double q, double radius, double sld, double solvent_sld);


double _FCC_Integrand(double q, double dnn, double d_factor, double theta, double phi) {

	const double Da = d_factor*dnn;
	const double temp1 = q*q*Da*Da;
	const double temp3 = q*dnn;

	double retVal = _FCCeval(theta,phi,temp1,temp3)/(4.0*M_PI);
	return(retVal);
}

double _FCCeval(double Theta, double Phi, double temp1, double temp3) {

	double temp6,temp7,temp8,temp9,temp10;
	double result;

	temp6 = sin(Theta);
	temp7 = sin(Theta)*sin(Phi)+cos(Theta);
	temp8 = -1.0*sin(Theta)*cos(Phi)+cos(Theta);
	temp9 = -1.0*sin(Theta)*cos(Phi)+sin(Theta)*sin(Phi);
	temp10 = exp((-1.0/8.0)*temp1*((temp7*temp7)+(temp8*temp8)+(temp9*temp9)));
	result = pow((1.0-(temp10*temp10)),3)*temp6/((1.0-2.0*temp10*cos(0.5*temp3*(temp7))+(temp10*temp10))*(1.0-2.0*temp10*cos(0.5*temp3*(temp8))+(temp10*temp10))*(1.0-2.0*temp10*cos(0.5*temp3*(temp9))+(temp10*temp10)));

	return (result);
}

double _sphereform(double q, double radius, double sld, double solvent_sld){
    const double qr = q*radius;
    double sn, cn;
    SINCOS(qr, sn, cn);
    const double bes = (qr == 0.0 ? 1.0 : 3.0*(sn-qr*cn)/(qr*qr*qr));
    const double fq = bes * (sld - solvent_sld)*form_volume(radius);
    return 1.0e-4*fq*fq;
}

double form_volume(double radius){
    return 1.333333333333333*M_PI*radius*radius*radius;
}


double Iq(double q, double dnn,
  double d_factor, double radius,
  double sld, double solvent_sld){

	//Volume fraction calculated from lattice symmetry and sphere radius
	const double s1 = dnn*sqrt(2.0);
	const double latticescale = 4.0*(4.0/3.0)*M_PI*(radius*radius*radius)/(s1*s1*s1);

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
	answer = answer*_sphereform(q,radius,sld,solvent_sld)*latticescale;

    return answer;


}


double Iqxy(double qx, double qy, double dnn,
    double d_factor, double radius,double sld, double solvent_sld,
    double theta, double phi, double psi){

  double b3_x, b3_y, b1_x, b1_y, b2_x, b2_y; //b3_z,
  double q_z;
  double cos_val_b3, cos_val_b2, cos_val_b1;
  double a1_dot_q, a2_dot_q,a3_dot_q;
  double answer;
  double Zq, Fkq, Fkq_2;

  //convert to q and make scaled values
  double q = sqrt(qx*qx+qy*qy);
  double q_x = qx/q;
  double q_y = qy/q;

  //convert angle degree to radian
  theta = theta * M_PI_180;
  phi = phi * M_PI_180;
  psi = psi * M_PI_180;

  const double Da = d_factor*dnn;
  const double s1 = dnn/sqrt(0.75);


  //the occupied volume of the lattice
  const double latticescale = 2.0*(4.0/3.0)*M_PI*(radius*radius*radius)/(s1*s1*s1);
  // q vector
  q_z = 0.0; // for SANS; assuming qz is negligible
  /// Angles here are respect to detector coordinate
  ///  instead of against q coordinate(PRB 36(46), 3(6), 1754(3854))
    // b3 axis orientation
    b3_x = cos(theta) * cos(phi);
    b3_y = sin(theta);
    //b3_z = -cos(theta) * sin(phi);
    cos_val_b3 =  b3_x*q_x + b3_y*q_y;// + b3_z*q_z;

    //alpha = acos(cos_val_b3);
    // b1 axis orientation
    b1_x = -cos(phi)*sin(psi) * sin(theta)+sin(phi)*cos(psi);
    b1_y = sin(psi)*cos(theta);
    cos_val_b1 = b1_x*q_x + b1_y*q_y;
    // b2 axis orientation
    b2_x = -sin(theta)*cos(psi)*cos(phi)-sin(psi)*sin(phi);
  	b2_y = cos(theta)*cos(psi);
    cos_val_b2 = b2_x*q_x + b2_y*q_y;

    // The following test should always pass
    if (fabs(cos_val_b3)>1.0) {
      //printf("FCC_ana_2D: Unexpected error: cos()>1\n");
      cos_val_b3 = 1.0;
    }
    if (fabs(cos_val_b2)>1.0) {
      //printf("FCC_ana_2D: Unexpected error: cos()>1\n");
      cos_val_b2 = 1.0;
    }
    if (fabs(cos_val_b1)>1.0) {
      //printf("FCC_ana_2D: Unexpected error: cos()>1\n");
      cos_val_b1 = 1.0;
    }
    // Compute the angle btw vector q and the a3 axis
    a3_dot_q = 0.5*dnn*q*(cos_val_b2+cos_val_b1-cos_val_b3);

    // a1 axis
    a1_dot_q = 0.5*dnn*q*(cos_val_b3+cos_val_b2-cos_val_b1);

    // a2 axis
    a2_dot_q = 0.5*dnn*q*(cos_val_b3+cos_val_b1-cos_val_b2);


    // Get Fkq and Fkq_2
    Fkq = exp(-0.5*pow(Da/dnn,2.0)*(a1_dot_q*a1_dot_q+a2_dot_q*a2_dot_q+a3_dot_q*a3_dot_q));
    Fkq_2 = Fkq*Fkq;
    // Call Zq=Z1*Z2*Z3
    Zq = (1.0-Fkq_2)/(1.0-2.0*Fkq*cos(a1_dot_q)+Fkq_2);
    Zq *= (1.0-Fkq_2)/(1.0-2.0*Fkq*cos(a2_dot_q)+Fkq_2);
    Zq *= (1.0-Fkq_2)/(1.0-2.0*Fkq*cos(a3_dot_q)+Fkq_2);

  // Use SphereForm directly from libigor
  answer = _sphereform(q,radius,sld,solvent_sld)*Zq*latticescale;

  return answer;
 }