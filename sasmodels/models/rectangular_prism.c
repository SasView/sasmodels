double form_volume(double length_a, double b2a_ratio, double c2a_ratio);
double Iq(double q, double sld, double solvent_sld, double length_a, 
          double b2a_ratio, double c2a_ratio);

double form_volume(double length_a, double b2a_ratio, double c2a_ratio)
{
    return length_a * (length_a*b2a_ratio) * (length_a*c2a_ratio);
}

double Iq(double q,
    double sld,
    double solvent_sld,
    double length_a,
    double b2a_ratio,
    double c2a_ratio)
{
    double termA, termB, termC;
    
    double b_side = length_a * b2a_ratio;
    double c_side = length_a * c2a_ratio;
    double volume = length_a * b_side * c_side;
    double a_half = 0.5 * length_a;
    double b_half = 0.5 * b_side;
    double c_half = 0.5 * c_side;

   //Integration limits to use in Gaussian quadrature
    double v1a = 0.0;
    double v1b = M_PI_2;  //theta integration limits
    double v2a = 0.0;
    double v2b = M_PI_2;  //phi integration limits
    
    //Order of integration
    int nordi=76;			        
    int nordj=76;

    double sumi = 0.0;
    
    for(int i=0; i<nordi; i++) {

	    double theta = 0.5 * ( Gauss76Z[i]*(v1b-v1a) + v1a + v1b );	

	    double arg = q * c_half * cos(theta);
	    if (fabs(arg) > 1.e-16) {termC = sin(arg)/arg;} else {termC = 1.0;}  

	    double sumj = 0.0;
        
	    for(int j=0; j<nordj; j++) {

            double phi = 0.5 * ( Gauss76Z[j]*(v2b-v2a) + v2a + v2b ); 

	        // Amplitude AP from eqn. (12), rewritten to avoid round-off effects when arg=0

	        arg = q * a_half * sin(theta) * sin(phi); 
	        if (fabs(arg) > 1.e-16) {termA = sin(arg)/arg;} else {termA = 1.0;}
	       
	        arg = q * b_half * sin(theta) * cos(phi); 
	        if (fabs(arg) > 1.e-16) {termB = sin(arg)/arg;} else {termB = 1.0;}	  
               
	        double AP = termA * termB * termC;  

	        sumj += Gauss76Wt[j] * (AP*AP);

	    }

	    sumj = 0.5 * (v2b-v2a) * sumj;
	    sumi += Gauss76Wt[i] * sumj * sin(theta);

    }

    double answer = 0.5*(v1b-v1a)*sumi;

    // Normalize by Pi (Eqn. 16). 
    // The term (ABC)^2 does not appear because it was introduced before on 
    // the definitions of termA, termB, termC.
    // The factor 2 appears because the theta integral has been defined between 
    // 0 and pi/2, instead of 0 to pi.
    answer *= (2.0/M_PI); //Form factor P(q)

    // Multiply by contrast^2 and volume^2
    answer *= (sld-solvent_sld)*(sld-solvent_sld)*volume*volume;

    // Convert from [1e-12 A-1] to [cm-1] 
    answer *= 1.0e-4;

    return answer;
    
}
