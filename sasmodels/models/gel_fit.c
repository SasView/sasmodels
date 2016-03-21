static double Iq(double q,
          double guinier_scale,
          double lorentzian_scale,
          double gyration_radius,
          double fractal_exp,
          double cor_length)
{
    // Lorentzian Term
    ////////////////////////double a(x[i]*x[i]*zeta*zeta);
    double lorentzian_term = square(q*cor_length);
    lorentzian_term = 1.0 + ((fractal_exp + 1.0)/3.0)*lorentzian_term;
    lorentzian_term = pow(lorentzian_term, (fractal_exp/2.0) );

    // Exponential Term
    ////////////////////////double d(x[i]*x[i]*rg*rg);
    double exp_term = square(q*gyration_radius);
    exp_term = exp(-1.0*(exp_term/3.0));

    // Scattering Law
    double result = lorentzian_scale/lorentzian_term + guinier_scale*exp_term;
    return result;
}
