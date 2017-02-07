double form_volume(double length_a, double length_b, double length_c);
double Iq(double q, double sld, double solvent_sld,
    double length_a, double length_b, double length_c);
double Iqxy(double qx, double qy, double sld, double solvent_sld,
    double length_a, double length_b, double length_c,
    double theta, double phi, double psi);

double form_volume(double length_a, double length_b, double length_c)
{
    return length_a * length_b * length_c;
}


double Iq(double q,
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

    for( int i=0; i<76; i++) {
        const double sigma = 0.5 * ( Gauss76Z[i] + 1.0 );
        const double mu_proj = mu * sqrt(1.0-sigma*sigma);

        // inner integral (with gauss points), integration limits = 0, 1
        // corresponding to angles from 0 to pi/2.
        double inner_total = 0.0;
        for(int j=0; j<76; j++) {
            const double uu = 0.5 * ( Gauss76Z[j] + 1.0 );
            double sin_uu, cos_uu;
            SINCOS(M_PI_2*uu, sin_uu, cos_uu);
            const double si1 = sas_sinx_x(mu_proj * sin_uu * a_scaled);
            const double si2 = sas_sinx_x(mu_proj * cos_uu);
            inner_total += Gauss76Wt[j] * square(si1 * si2);
        }
        inner_total *= 0.5;

        const double si = sas_sinx_x(mu * c_scaled * sigma);
        outer_total += Gauss76Wt[i] * inner_total * si * si;
    }
    outer_total *= 0.5;

    // Multiply by contrast^2 and convert from [1e-12 A-1] to [cm-1]
    const double V = form_volume(length_a, length_b, length_c);
    const double drho = (sld-solvent_sld);
    return 1.0e-4 * square(drho * V) * outer_total;
}


double Iqxy(double qx, double qy,
    double sld,
    double solvent_sld,
    double length_a,
    double length_b,
    double length_c,
    double theta,
    double phi,
    double psi)
{
    double q, cos_val_a, cos_val_b, cos_val_c;
    ORIENT_ASYMMETRIC(qx, qy, theta, phi, psi, q, cos_val_c, cos_val_b, cos_val_a);

    const double siA = sas_sinx_x(0.5*q*length_a*cos_val_a);
    const double siB = sas_sinx_x(0.5*q*length_b*cos_val_b);
    const double siC = sas_sinx_x(0.5*q*length_c*cos_val_c);
    const double V = form_volume(length_a, length_b, length_c);
    const double drho = (sld - solvent_sld);
    const double form = V * drho * siA * siB * siC;
    // Square and convert from [1e-12 A-1] to [cm-1]
    return 1.0e-4 * form * form;
}
