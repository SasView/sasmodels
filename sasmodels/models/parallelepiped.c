static double
form_volume(double length_a, double length_b, double length_c)
{
    return length_a * length_b * length_c;
}

static double
effective_radius(int mode, double length_a, double length_b, double length_c)
{
    if (mode == 1) {
        return cbrt(0.75*length_a*length_b*length_c/M_PI);
    } else if (mode == 2) {
        return 0.5 * length_a;
    } else if (mode == 3) {
        return 0.5 * length_b;
    } else if (mode == 4) {
        return 0.5 * length_c;
    } else if (mode == 5) {
        return sqrt(length_a*length_b/M_PI);
    } else if (mode == 6) {
        return 0.5*sqrt(length_a*length_a + length_b*length_b);
    } else {
        return 0.5*sqrt(length_a*length_a + length_b*length_b + length_c*length_c);
    }
}


static void
Fq(double q,
    double *F1,
    double *F2,
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
    double outer_total_F1 = 0.0; //initialize integral
    double outer_total_F2 = 0.0; //initialize integral
    for( int i=0; i<GAUSS_N; i++) {
        const double sigma = 0.5 * ( GAUSS_Z[i] + 1.0 );
        const double mu_proj = mu * sqrt(1.0-sigma*sigma);

        // inner integral (with gauss points), integration limits = 0, 1
        // corresponding to angles from 0 to pi/2.
        double inner_total_F1 = 0.0;
        double inner_total_F2 = 0.0;
        for(int j=0; j<GAUSS_N; j++) {
            const double uu = 0.5 * ( GAUSS_Z[j] + 1.0 );
            double sin_uu, cos_uu;
            SINCOS(M_PI_2*uu, sin_uu, cos_uu);
            const double si1 = sas_sinx_x(mu_proj * sin_uu * a_scaled);
            const double si2 = sas_sinx_x(mu_proj * cos_uu);
            const double fq = si1 * si2;
            inner_total_F1 += GAUSS_W[j] * fq;
            inner_total_F2 += GAUSS_W[j] * fq * fq;
        }
        // now complete change of inner integration variable (1-0)/(1-(-1))= 0.5
        inner_total_F1 *= 0.5;
        inner_total_F2 *= 0.5;

        const double si = sas_sinx_x(mu * c_scaled * sigma);
        outer_total_F1 += GAUSS_W[i] * inner_total_F1 * si;
        outer_total_F2 += GAUSS_W[i] * inner_total_F2 * si * si;
    }
    // now complete change of outer integration variable (1-0)/(1-(-1))= 0.5
    outer_total_F1 *= 0.5;
    outer_total_F2 *= 0.5;

    // Multiply by contrast^2 and convert from [1e-12 A-1] to [cm-1]
    const double V = form_volume(length_a, length_b, length_c);
    const double contrast = (sld-solvent_sld);
    const double s = contrast * V;
    *F1 = 1.0e-2 * s * outer_total_F1;
    *F2 = 1.0e-4 * s * s * outer_total_F2;
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
