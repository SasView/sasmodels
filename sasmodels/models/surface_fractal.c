double form_volume(double radius);

double Iq(double q, double radius, double surface_dim, double cutoff_length);
double Iqxy(double qx, double qy, double radius, double surface_dim, double cutoff_length);

static double _gamln(double q)
{
    // Lanczos approximation to the Gamma function.
    // Should be refactored out to lib/, if used elsewhere.
    double x,y,tmp,ser;
    double coeff[6]=
        {76.18009172947146,     -86.50532032941677,
         24.01409824083091,     -1.231739572450155,
          0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=q;
    tmp  = x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser  = 1.000000000190015;
    for (j=0; j<=5; j++) {
        y+=1.0;
        ser += coeff[j]/y;
    }
    return -tmp+log(2.5066282746310005*ser/x);
}

static double surface_fractal_kernel(double q,
    double radius,
    double surface_dim,
    double cutoff_length)
{
    double pq, sq, mmo, result;

    //calculate P(q) for the spherical subunits; not normalized
	pq = pow((3.0*(sin(q*radius) - q*radius*cos(q*radius))/pow((q*radius),3)),2);

    //calculate S(q)
    mmo = 5.0 - surface_dim;
    sq  = exp(_gamln(mmo))*sin(-(mmo)*atan(q*cutoff_length));
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
    return surface_fractal_kernel(q, radius, surface_dim, cutoff_length);
}

// Iqxy is never called since no orientation or magnetic parameters.
double Iqxy(double qx, double qy,
    double radius,
    double surface_dim,
    double cutoff_length)
{
    double q = sqrt(qx*qx + qy*qy);
    return surface_fractal_kernel(q, radius, surface_dim, cutoff_length);
}

