double _hollow_cylinder_kernel(double q, double core_radius, double radius, 
	double length, double dum);

double form_volume(double radius, double core_radius, double length);
double Iq(double q, double radius, double core_radius, double length, double sld,
	double solvent_sld);
double Iqxy(double qx, double qy, double radius, double core_radius, double length, double sld,
	double solvent_sld, double theta, double phi);

// From Igor library
double _hollow_cylinder_kernel(double q, double core_radius, double radius, 
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
	double norm,scale,volume,convert;	//final calculation variables
	
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
	// Multiply by contrast^2
	scale = delrho*delrho;
	//normalize by volume
	volume = form_volume(radius, core_radius, length);
	//convert to [cm-1] given sld*1e6
	convert = 1.0e-4;
	answer = norm*scale*convert*volume*volume;
	
	return(answer);
}

//TODO: Add this in
double Iqxy(double qx, double qy, double radius, double core_radius, double length, double sld,
	double solvent_sld, double theta, double phi)
{
    return(0.0);
}
