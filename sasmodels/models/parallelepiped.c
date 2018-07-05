static double
form_volume(double length_a, double length_b, double length_c)
{
    return length_a * length_b * length_c;
}


static double
Iq(double q,
    double sld,
    double solvent_sld,
    double length_a,
    double length_b,
    double length_c)
{
    const double mu = 0.5 * q * length_b;

    // Scale sides by B
    const double a_scaled = length_a / length_b;
    const double c_scaled = length_c / length_b;

    // outer integral (with gauss points), integration limits = 0, 1
    double outer_total = 0; //initialize integral

    for( int i=0; i<GAUSS_N; i++) {
        const double sigma = 0.5 * ( GAUSS_Z[i] + 1.0 );
        const double mu_proj = mu * sqrt(1.0-sigma*sigma);

        // inner integral (with gauss points), integration limits = 0, 1
        // corresponding to angles from 0 to pi/2.
        double inner_total = 0.0;
        for(int j=0; j<GAUSS_N; j++) {
            const double uu = 0.5 * ( GAUSS_Z[j] + 1.0 );
            double sin_uu, cos_uu;
            SINCOS(M_PI_2*uu, sin_uu, cos_uu);
            const double si1 = sas_sinx_x(mu_proj * sin_uu * a_scaled);
            const double si2 = sas_sinx_x(mu_proj * cos_uu);
            inner_total += GAUSS_W[j] * square(si1 * si2);
        }
        // now complete change of inner integration variable (1-0)/(1-(-1))= 0.5
        inner_total *= 0.5;

        const double si = sas_sinx_x(mu * c_scaled * sigma);
        outer_total += GAUSS_W[i] * inner_total * si * si;
    }
    // now complete change of outer integration variable (1-0)/(1-(-1))= 0.5
    outer_total *= 0.5;

    // Multiply by contrast^2 and convert from [1e-12 A-1] to [cm-1]
    const double V = form_volume(length_a, length_b, length_c);
    const double drho = (sld-solvent_sld);
    return 1.0e-4 * square(drho * V) * outer_total;
}


static double
Iqabc(double qa, double qb, double qc,
    double sld,
    double solvent_sld,
    double length_a,
    double length_b,
    double length_c)
{
    const double siA = sas_sinx_x(0.5*length_a*qa);
    const double siB = sas_sinx_x(0.5*length_b*qb);
    const double siC = sas_sinx_x(0.5*length_c*qc);
    const double V = form_volume(length_a, length_b, length_c);
    const double drho = (sld - solvent_sld);
    const double form = V * drho * siA * siB * siC;
    // Square and convert from [1e-12 A-1] to [cm-1]
    return 1.0e-4 * form * form;
}
