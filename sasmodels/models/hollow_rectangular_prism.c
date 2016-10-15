double form_volume(double length_a, double b2a_ratio, double c2a_ratio, double thickness);
double Iq(double q, double sld, double solvent_sld, double length_a, 
          double b2a_ratio, double c2a_ratio, double thickness);

double form_volume(double length_a, double b2a_ratio, double c2a_ratio, double thickness)
{
    double length_b = length_a * b2a_ratio;
    double length_c = length_a * c2a_ratio;
    double a_core = length_a - 2.0*thickness;
    double b_core = length_b - 2.0*thickness;
    double c_core = length_c - 2.0*thickness;
    double vol_core = a_core * b_core * c_core;
    double vol_total = length_a * length_b * length_c;
    double vol_shell = vol_total - vol_core;
    return vol_shell;
}

double Iq(double q,
    double sld,
    double solvent_sld,
    double length_a,
    double b2a_ratio,
    double c2a_ratio,
    double thickness)
{
    double termA1, termA2, termB1, termB2, termC1, termC2;
    
    double length_b = length_a * b2a_ratio;
    double length_c = length_a * c2a_ratio;
    double a_half = 0.5 * length_a;
    double b_half = 0.5 * length_b;
    double c_half = 0.5 * length_c;
    double vol_total = length_a * length_b * length_c;
    double vol_core = 8.0 * (a_half-thickness) * (b_half-thickness) * (c_half-thickness);

    //Integration limits to use in Gaussian quadrature
    double v1a = 0.0;
    double v1b = M_PI_2;  //theta integration limits
    double v2a = 0.0;
    double v2b = M_PI_2;  //phi integration limits
    
    double outer_sum = 0.0;
    
    for(int i=0; i<76; i++) {

        double theta = 0.5 * ( Gauss76Z[i]*(v1b-v1a) + v1a + v1b );    

        double termC1 = sinc(q * c_half * cos(theta));
        double termC2 = sinc(q * (c_half-thickness)*cos(theta));

        double inner_sum = 0.0;
        
        for(int j=0; j<76; j++) {

            double phi = 0.5 * ( Gauss76Z[j]*(v2b-v2a) + v2a + v2b ); 

            // Amplitude AP from eqn. (13), rewritten to avoid round-off effects when arg=0

            termA1 = sinc(q * a_half * sin(theta) * sin(phi));
            termA2 = sinc(q * (a_half-thickness) * sin(theta) * sin(phi));

            termB1 = sinc(q * b_half * sin(theta) * cos(phi));
            termB2 = sinc(q * (b_half-thickness) * sin(theta) * cos(phi));

            double AP1 = vol_total * termA1 * termB1 * termC1;
            double AP2 = vol_core * termA2 * termB2 * termC2;

            inner_sum += Gauss76Wt[j] * square(AP1-AP2);

        }

        inner_sum = 0.5 * (v2b-v2a) * inner_sum;
        outer_sum += Gauss76Wt[i] * inner_sum * sin(theta);

    }

    double answer = 0.5*(v1b-v1a)*outer_sum;

    // Normalize as in Eqn. (15) without the volume factor (as cancels with (V*DelRho)^2 normalization)
    // The factor 2 is due to the different theta integration limit (pi/2 instead of pi)
    answer /= M_PI_2;

    // Multiply by contrast^2. Factor corresponding to volume^2 cancels with previous normalization.
    answer *= square(sld-solvent_sld);

    // Convert from [1e-12 A-1] to [cm-1]
    answer *= 1.0e-4;

    return answer;
    
}
