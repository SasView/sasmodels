double form_volume(double a_side, double b2a_ratio, double c2a_ratio, double thickness);
double Iq(double q, double sld, double solvent_sld, double a_side, 
          double b2a_ratio, double c2a_ratio, double thickness);
double Iqxy(double qx, double qy, double sld, double solvent_sld, 
            double a_side, double b2a_ratio, double c2a_ratio, double thickness);

double form_volume(double a_side, double b2a_ratio, double c2a_ratio, double thickness)
{
    double b_side = a_side * b2a_ratio;
    double c_side = a_side * c2a_ratio;
    double a_core = a_side - 2.0*thickness;
    double b_core = b_side - 2.0*thickness;
    double c_core = c_side - 2.0*thickness;
    double vol_core = a_core * b_core * c_core;
    double vol_total = a_side * b_side * c_side;
    double vol_shell = vol_total - vol_core;
    return vol_shell;
}

double Iq(double q,
    double sld,
    double solvent_sld,
    double a_side,
    double b2a_ratio,
    double c2a_ratio,
    double thickness)
{
    double termA1, termA2, termB1, termB2, termC1, termC2;
    
    double b_side = a_side * b2a_ratio;
    double c_side = a_side * c2a_ratio;
    double a_half = 0.5 * a_side;
    double b_half = 0.5 * b_side;
    double c_half = 0.5 * c_side;

   //Integration limits to use in Gaussian quadrature
    double v1a = 0.0;
    double v1b = 0.5 * M_PI;  //theta integration limits
    double v2a = 0.0;
    double v2b = 0.5 * M_PI;  //phi integration limits
    
    //Order of integration
    int nordi=76;			        
    int nordj=76;

    double sumi = 0.0;
    
    for(int i=0; i<nordi; i++) {

	    double theta = 0.5 * ( Gauss76Z[i]*(v1b-v1a) + v1a + v1b );	

	    double arg = q * c_half * cos(theta);
	    if (fabs(arg) > 1.e-16) {termC1 = sin(arg)/arg;} else {termC1 = 1.0;}
	    arg = q * (c_half-thickness)*cos(theta);
	    if (fabs(arg) > 1.e-16) {termC2 = sin(arg)/arg;} else {termC2 = 1.0;}

	    double sumj = 0.0;
        
	    for(int j=0; j<nordj; j++) {

            double phi = 0.5 * ( Gauss76Z[j]*(v2b-v2a) + v2a + v2b ); 

            // Amplitude AP from eqn. (13), rewritten to avoid round-off effects when arg=0

	        arg = q * a_half * sin(theta) * sin(phi);
	        if (fabs(arg) > 1.e-16) {termA1 = sin(arg)/arg;} else {termA1 = 1.0;}
	        arg = q * (a_half-thickness) * sin(theta) * sin(phi);
	        if (fabs(arg) > 1.e-16) {termA2 = sin(arg)/arg;} else {termA2 = 1.0;}

	        arg = q * b_half * sin(theta) * cos(phi);
	        if (fabs(arg) > 1.e-16) {termB1 = sin(arg)/arg;} else {termB1 = 1.0;}
	        arg = q * (b_half-thickness) * sin(theta) * cos(phi);
	        if (fabs(arg) > 1.e-16) {termB2 = sin(arg)/arg;} else {termB2 = 1.0;}

            double AP1 = (a_side*b_side*c_side) * termA1 * termB1 * termC1;
            double AP2 = 8.0 * (a_half-thickness) * (b_half-thickness) * (c_half-thickness) * termA2 * termB2 * termC2;
            double AP = AP1 - AP2;

	        sumj += Gauss76Wt[j] * (AP*AP);

	    }

	    sumj = 0.5 * (v2b-v2a) * sumj;
	    sumi += Gauss76Wt[i] * sumj * sin(theta);

    }

    double answer = 0.5*(v1b-v1a)*sumi;

    // Normalize as in Eqn. (15) without the volume factor (as cancels with (V*DelRho)^2 normalization)
    // The factor 2 is due to the different theta integration limit (pi/2 instead of pi)
    answer *= (2.0/M_PI);

    // Multiply by contrast^2. Factor corresponding to volume^2 cancels with previous normalization.
    answer *= (sld-solvent_sld)*(sld-solvent_sld);

    // Convert from [1e-12 A-1] to [cm-1]
    answer *= 1.0e-4;

    return answer;
    
}

double Iqxy(double qx, double qy,
    double sld,
    double solvent_sld,
    double a_side,
    double b2a_ratio,
    double c2a_ratio,
    double thickness)
{
    double q = sqrt(qx*qx + qy*qy);
    double intensity = Iq(q, sld, solvent_sld, a_side, b2a_ratio, c2a_ratio, thickness); 
    return intensity;    
}