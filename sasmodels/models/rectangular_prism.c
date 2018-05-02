static double
form_volume(double length_a, double b2a_ratio, double c2a_ratio)
{
    return length_a * (length_a*b2a_ratio) * (length_a*c2a_ratio);
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
        SINCOS(theta, sin_theta, cos_theta);

        const double termC = sas_sinx_x(q * c_half * cos_theta);

        double inner_sum = 0.0;
        for(int j=0; j<GAUSS_N; j++) {
            double phi = 0.5 * ( GAUSS_Z[j]*(v2b-v2a) + v2a + v2b );
            double sin_phi, cos_phi;
            SINCOS(phi, sin_phi, cos_phi);

            // Amplitude AP from eqn. (12), rewritten to avoid round-off effects when arg=0
            const double termA = sas_sinx_x(q * a_half * sin_theta * sin_phi);
            const double termB = sas_sinx_x(q * b_half * sin_theta * cos_phi);
            const double AP = termA * termB * termC;
            inner_sum += GAUSS_W[j] * AP * AP;
        }
        inner_sum = 0.5 * (v2b-v2a) * inner_sum;
        outer_sum += GAUSS_W[i] * inner_sum * sin_theta;
    }

    double answer = 0.5*(v1b-v1a)*outer_sum;

    // Normalize by Pi (Eqn. 16).
    // The term (ABC)^2 does not appear because it was introduced before on
    // the definitions of termA, termB, termC.
    // The factor 2 appears because the theta integral has been defined between
    // 0 and pi/2, instead of 0 to pi.
    answer /= M_PI_2; //Form factor P(q)

    // Multiply by contrast^2 and volume^2
    const double volume = length_a * length_b * length_c;
    answer *= square((sld-solvent_sld)*volume);

    // Convert from [1e-12 A-1] to [cm-1]
    answer *= 1.0e-4;

    return answer;
}


static double
Iqabc(double qa, double qb, double qc,
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
    const double volume = length_a * length_b * length_c;

    // Amplitude AP from eqn. (13)

    const double termA = sas_sinx_x(qa * a_half);
    const double termB = sas_sinx_x(qb * b_half);
    const double termC = sas_sinx_x(qc * c_half);

    const double AP = termA * termB * termC;

    // Multiply by contrast^2. Factor corresponding to volume^2 cancels with previous normalization.
    const double delrho = sld - solvent_sld;

    // Convert from [1e-12 A-1] to [cm-1]
    return 1.0e-4 * square(volume * delrho * AP);
}