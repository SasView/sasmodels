double form_volume(double radius);

double Iq(double q,
          double radius,
          double mass_dim,
          double cutoff_length);

double Iqxy(double qx, double qy,
          double radius,
          double mass_dim,
          double cutoff_length);


static double _mass_fractal_kernel(double q,
          double radius,
          double mass_dim,
          double cutoff_length)
{

    //calculate P(q) for the spherical subunits; not normalized
    //double pq = pow((3.0*(sin(q*radius) - q*radius*cos(q*radius))/pow((q*radius),3)),2);
    double pq = J1c(q*radius);

    //pq = 1.0;

    //calculate S(q)
    double mmo = mass_dim-1.0;
    double sq = exp(lanczos_gamma(mmo))*sin((mmo)*atan(q*cutoff_length));
    sq *= pow(cutoff_length, mmo);
    sq /= pow((1.0 + (q*cutoff_length)*(q*cutoff_length)),(mmo/2.0));
    sq /= q;

    //sq = 1.0;

    //combine and return
    double result = pq * sq;

    return result;
}
double form_volume(double radius){

    return 1.333333333333333*M_PI*radius*radius*radius;
}

double Iq(double q,
          double radius,
          double mass_dim,
          double cutoff_length)
{
    return _mass_fractal_kernel(q,
           radius,
           mass_dim,
           cutoff_length);
}

// Iqxy is never called since no orientation or magnetic parameters.
double Iqxy(double qx, double qy,
          double radius,
          double mass_dim,
          double cutoff_length)
{
    double q = sqrt(qx*qx + qy*qy);
    return _mass_fractal_kernel(q,
           radius,
           mass_dim,
           cutoff_length);
}

