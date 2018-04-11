static double
form_volume(double length_a, double b2a_ratio, double c2a_ratio)
{
    double length_b = length_a * b2a_ratio;
    double length_c = length_a * c2a_ratio;
    double vol_shell = 2.0 * (length_a*length_b + length_a*length_c + length_b*length_c);
    return vol_shell;
}

static double
Iq(double q,
    double sld,
    double solvent_sld,
    double length_a,
    double b2a_ratio,
    double c2a_ratio)
{
    const double length_b = length_a * b2a_ratio;
    const double length_c = length_a * c2a_ratio;
    const double a_half = 0.5 * length_a;
    const double b_half = 0.5 * length_b;
    const double c_half = 0.5 * length_c;

   //Integration limits to use in Gaussian quadrature
    const double v1a = 0.0;
    const double v1b = M_PI_2;  //theta integration limits
    const double v2a = 0.0;
    const double v2b = M_PI_2;  //phi integration limits

    double outer_sum = 0.0;
    for(int i=0; i<GAUSS_N; i++) {
        const double theta = 0.5 * ( GAUSS_Z[i]*(v1b-v1a) + v1a + v1b );

        double sin_theta, cos_theta;
        double sin_c, cos_c;
        SINCOS(theta, sin_theta, cos_theta);
        SINCOS(q*c_half*cos_theta, sin_c, cos_c);

        // To check potential problems if denominator goes to zero here !!!
        const double termAL_theta = 8.0 * cos_c / (q*q*sin_theta*sin_theta);
        const double termAT_theta = 8.0 * sin_c / (q*q*sin_theta*cos_theta);

        double inner_sum = 0.0;
        for(int j=0; j<GAUSS_N; j++) {
            const double phi = 0.5 * ( GAUSS_Z[j]*(v2b-v2a) + v2a + v2b );

            double sin_phi, cos_phi;
            double sin_a, cos_a;
            double sin_b, cos_b;
            SINCOS(phi, sin_phi, cos_phi);
            SINCOS(q*a_half*sin_theta*sin_phi, sin_a, cos_a);
            SINCOS(q*b_half*sin_theta*cos_phi, sin_b, cos_b);

            // Amplitude AL from eqn. (7c)
            const double AL = termAL_theta
                * sin_a*sin_b / (sin_phi*cos_phi);

            // Amplitude AT from eqn. (9)
            const double AT = termAT_theta
                * ( cos_a*sin_b/cos_phi + cos_b*sin_a/sin_phi );

            inner_sum += GAUSS_W[j] * square(AL+AT);
        }

        inner_sum *= 0.5 * (v2b-v2a);
        outer_sum += GAUSS_W[i] * inner_sum * sin_theta;
    }

    outer_sum *= 0.5*(v1b-v1a);

    // Normalize as in Eqn. (15) without the volume factor (as cancels with (V*DelRho)^2 normalization)
    // The factor 2 is due to the different theta integration limit (pi/2 instead of pi)
    double answer = outer_sum/M_PI_2;

    // Multiply by contrast^2. Factor corresponding to volume^2 cancels with previous normalization.
    answer *= square(sld-solvent_sld);

    // Convert from [1e-12 A-1] to [cm-1]
    answer *= 1.0e-4;

    return answer;
}
