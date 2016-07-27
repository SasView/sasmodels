double form_volume(double radius);

double Iq(double q,
          double radius,
          double surface_dim,
          double cutoff_length);

static double _surface_fractal_kernel(double q,
    double radius,
    double surface_dim,
    double cutoff_length)
{
    double pq, sq, mmo, result;

    //Replaced the original formula with Taylor expansion near zero.
    //pq = pow((3.0*(sin(q*radius) - q*radius*cos(q*radius))/pow((q*radius),3)),2);

    pq = sph_j1c(q*radius);
    pq = pq*pq;

    //calculate S(q)
    mmo = 5.0 - surface_dim;
    sq  = sas_gamma(mmo)*sin(-(mmo)*atan(q*cutoff_length));
    sq *= pow(cutoff_length, mmo);
    sq /= pow((1.0 + (q*cutoff_length)*(q*cutoff_length)),(mmo/2.0));
    sq /= q;

    //combine and return
    result = pq * sq;

    return result;
}
double form_volume(double radius){

    return 1.333333333333333*M_PI*radius*radius*radius;
}

double Iq(double q,
    double radius,
    double surface_dim,
    double cutoff_length
    )
{
    return _surface_fractal_kernel(q, radius, surface_dim, cutoff_length);
}
