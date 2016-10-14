double form_volume(double radius);

double Iq(double q,
          double radius,
          double fractal_dim_mass,
          double cutoff_length);

static double _mass_fractal_kernel(double q,
          double radius,
          double fractal_dim_mass,
          double cutoff_length)
{
    // Actively check the argument.
    if (fractal_dim_mass <= 1.0){
       return 0.0;
    }

    //calculate P(q)
    double pq = sph_j1c(q*radius);
    pq = pq*pq;

    //calculate S(q)
    double mmo = fractal_dim_mass-1.0;
    double sq = sas_gamma(mmo)*sin((mmo)*atan(q*cutoff_length));
    sq *= pow(cutoff_length, mmo);
    sq /= pow((1.0 + (q*cutoff_length)*(q*cutoff_length)),(mmo/2.0));
    sq /= q;

    //combine and return
    double result = pq * sq;

    return result;
}
double form_volume(double radius){

    return M_4PI_3*radius*radius*radius;
}

double Iq(double q,
          double radius,
          double fractal_dim_mass,
          double cutoff_length)
{
    return _mass_fractal_kernel(q,
           radius,
           fractal_dim_mass,
           cutoff_length);
}
