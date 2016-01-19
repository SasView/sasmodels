double form_volume(void);

double Iq(double q,
          double guinier_scale,
          double lorentzian_scale,
          double gyration_radius,
          double fractal_exp,
          double cor_length);

double Iqxy(double qx, double qy,
          double guinier_scale,
          double lorentzian_scale,
          double gyration_radius,
          double fractal_exp,
          double cor_length);

static double _gel_fit_kernel(double q,
          double guinier_scale,
          double lorentzian_scale,
          double gyration_radius,
          double fractal_exp,
          double cor_length)
{
    // Lorentzian Term
    ////////////////////////double a(x[i]*x[i]*zeta*zeta);
    double lorentzian_term = q*q*cor_length*cor_length;
    lorentzian_term = 1.0 + ((fractal_exp + 1.0)/3.0)*lorentzian_term;
    lorentzian_term = pow(lorentzian_term, (fractal_exp/2.0) );

    // Exponential Term
    ////////////////////////double d(x[i]*x[i]*rg*rg);
    double exp_term = q*q*gyration_radius*gyration_radius;
    exp_term = exp(-1.0*(exp_term/3.0));

    // Scattering Law
    double result = lorentzian_scale/lorentzian_term + guinier_scale*exp_term;
    return result;
}
double form_volume(void){
    // Unused, so free to return garbage.
    return NAN;
}

double Iq(double q,
          double guinier_scale,
          double lorentzian_scale,
          double gyration_radius,
          double fractal_exp,
          double cor_length)
{
    return _gel_fit_kernel(q,
                          guinier_scale,
                          lorentzian_scale,
                          gyration_radius,
                          fractal_exp,
                          cor_length);
}

// Iqxy is never called since no orientation or magnetic parameters.
double Iqxy(double qx, double qy,
          double guinier_scale,
          double lorentzian_scale,
          double gyration_radius,
          double fractal_exp,
          double cor_length)
{
    double q = sqrt(qx*qx + qy*qy);
    return _gel_fit_kernel(q,
                          guinier_scale,
                          lorentzian_scale,
                          gyration_radius,
                          fractal_exp,
                          cor_length);
}

