double form_volume(double length_a, double b2a_ratio, double c2a_ratio);
double Iq(double q, double sld, double solvent_sld, double length_a, 
          double b2a_ratio, double c2a_ratio);

double form_volume(double length_a, double b2a_ratio, double c2a_ratio)
{
    double b_side = length_a * b2a_ratio;
    double c_side = length_a * c2a_ratio;
    double vol_shell = 2.0 * (length_a*b_side + length_a*c_side + b_side*c_side);
    return vol_shell;
}

double Iq(double q,
    double sld,
    double solvent_sld,
    double length_a,
    double b2a_ratio,
    double c2a_ratio)
{
    double b_side = length_a * b2a_ratio;
    double c_side = length_a * c2a_ratio;
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
        
        // To check potential problems if denominator goes to zero here !!!
        double termAL_theta = 8.0*cos(q*c_half*cos(theta)) / (q*q*sin(theta)*sin(theta));
        double termAT_theta = 8.0*sin(q*c_half*cos(theta)) / (q*q*sin(theta)*cos(theta));

	    double sumj = 0.0;
        
	    for(int j=0; j<nordj; j++) {

            double phi = 0.5 * ( Gauss76Z[j]*(v2b-v2a) + v2a + v2b ); 
            
            // Amplitude AL from eqn. (7c)
            double AL = termAL_theta * sin(q*a_half*sin(theta)*sin(phi)) * 
                sin(q*b_half*sin(theta)*cos(phi)) / (sin(phi)*cos(phi));

            // Amplitude AT from eqn. (9)
            double AT = termAT_theta * (  (cos(q*a_half*sin(theta)*sin(phi))*sin(q*b_half*sin(theta)*cos(phi))/cos(phi)) 
                + (cos(q*b_half*sin(theta)*cos(phi))*sin(q*a_half*sin(theta)*sin(phi))/sin(phi)) );

            sumj += Gauss76Wt[j] * (AL+AT)*(AL+AT);

	    }

	    sumj = 0.5 * (v2b-v2a) * sumj;
	    sumi += Gauss76Wt[i] * sumj * sin(theta);

    }

    double answer = 0.5*(v1b-v1a)*sumi;

    // Normalize as in Eqn. (15) without the volume factor (as cancels with (V*DelRho)^2 normalization)
    // The factor 2 is due to the different theta integration limit (pi/2 instead of pi)
    answer /= M_PI_2;

    // Multiply by contrast^2. Factor corresponding to volume^2 cancels with previous normalization.
    answer *= (sld-solvent_sld)*(sld-solvent_sld);

    // Convert from [1e-12 A-1] to [cm-1]
    answer *= 1.0e-4;

    return answer;
    
}
