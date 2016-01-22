static double _hollow_cylinder_kernel(double q, double core_radius, double radius, 
	double length, double dum);
static double hollow_cylinder_analytical_2D_scaled(double q, double q_x, double q_y, double radius, double core_radius, double length, double sld,
	double solvent_sld, double theta, double phi);
static double hollow_cylinder_scaling(double integrand, double delrho, double volume);
	
double form_volume(double radius, double core_radius, double length);

double Iq(double q, double radius, double core_radius, double length, double sld,
	double solvent_sld);
double Iqxy(double qx, double qy, double radius, double core_radius, double length, double sld,
	double solvent_sld, double theta, double phi);

// From Igor library
static double _hollow_cylinder_kernel(double q, double core_radius, double radius, 
	double length, double dum)
{
    double gamma,arg1,arg2,lam1,lam2,psi,sinarg,t2,retval;		//local variables
    
    gamma = core_radius/radius;
    arg1 = q*radius*sqrt(1.0-dum*dum);		//1=shell (outer radius)
    arg2 = q*core_radius*sqrt(1.0-dum*dum);			//2=core (inner radius)
    if (arg1 == 0.0){
    	lam1 = 1.0;
    }else{
    	lam1 = 2.0*J1(arg1)/arg1;
    }
    if (arg2 == 0.0){
    	lam2 = 1.0;
    }else{
    	lam2 = 2.0*J1(arg2)/arg2;
    }
    //Todo: Need to check psi behavior as gamma goes to 1.
    psi = (lam1 -  gamma*gamma*lam2)/(1.0-gamma*gamma);		//SRK 10/19/00
    sinarg = q*length*dum/2.0;
    if (sinarg == 0.0){
    	t2 = 1.0;
    }else{
    	t2 = sin(sinarg)/sinarg;
    }

    retval = psi*psi*t2*t2;
    
    return(retval);
}
static double hollow_cylinder_analytical_2D_scaled(double q, double q_x, double q_y, double radius, double core_radius, double length, double sld,
	double solvent_sld, double theta, double phi) {
	double cyl_x, cyl_y; //, cyl_z
	//double q_z;
	double vol, cos_val, delrho;
	double answer;
	//convert angle degree to radian
	double pi = 4.0*atan(1.0);
	theta = theta * pi/180.0;
	phi = phi * pi/180.0;
	delrho = solvent_sld - sld;

	// Cylinder orientation
	cyl_x = cos(theta) * cos(phi);
	cyl_y = sin(theta);
	//cyl_z = -cos(theta) * sin(phi);

	// q vector
	//q_z = 0;

	// Compute the angle btw vector q and the
	// axis of the cylinder
	cos_val = cyl_x*q_x + cyl_y*q_y;// + cyl_z*q_z;

	// The following test should always pass
	if (fabs(cos_val)>1.0) {
		//printf("core_shell_cylinder_analytical_2D: Unexpected error: cos(alpha)=%g\n", cos_val);
		return NAN;
	}

	answer = _hollow_cylinder_kernel(q, core_radius, radius, length, cos_val);

	vol = form_volume(radius, core_radius, length);
	answer = hollow_cylinder_scaling(answer, delrho, vol);

	return answer;
}
static double hollow_cylinder_scaling(double integrand, double delrho, double volume)
{
	double answer;
	// Multiply by contrast^2
	answer = integrand*delrho*delrho;

	//normalize by cylinder volume
	answer *= volume*volume;

	//convert to [cm-1]
	answer *= 1.0e-4;
	
	return answer;
}


double form_volume(double radius, double core_radius, double length)
{
	double pi = 4.0*atan(1.0);
	double v_shell = pi*length*(radius*radius-core_radius*core_radius);
	return(v_shell);
}


double Iq(double q, double radius, double core_radius, double length, double sld,
	double solvent_sld)
{
    int i;
	int nord=76;			//order of integration
	double lower,upper,zi, inter;		//upper and lower integration limits
	double summ,answer,delrho;			//running tally of integration
	double norm,volume;	//final calculation variables
	
	if (core_radius >= radius || radius >= length) {
		return NAN;
	}
	
	delrho = solvent_sld - sld;
	lower = 0.0;
	upper = 1.0;		//limits of numerical integral

	summ = 0.0;			//initialize intergral
	for(i=0;i<nord;i++) {
		zi = ( Gauss76Z[i] * (upper-lower) + lower + upper )/2.0;
		inter = Gauss76Wt[i] * _hollow_cylinder_kernel(q, core_radius, radius, length, zi);
		summ += inter;
	}
 	
	norm = summ*(upper-lower)/2.0;
	volume = form_volume(radius, core_radius, length);
	answer = hollow_cylinder_scaling(norm, delrho, volume);
	
	return(answer);
}

double Iqxy(double qx, double qy, double radius, double core_radius, double length, double sld,
	double solvent_sld, double theta, double phi)
{
	double q;
	q = sqrt(qx*qx+qy*qy);
	return hollow_cylinder_analytical_2D_scaled(q, qx/q, qy/q, radius, core_radius, length, sld, solvent_sld, theta, phi);
}