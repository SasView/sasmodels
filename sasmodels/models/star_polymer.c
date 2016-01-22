double form_volume(void);

double Iq(double q, double radius2, double arms);

double Iqxy(double qx, double qy, double radius2, double arms);


static double _mass_fractal_kernel(double q, double radius2, double arms)
{

    double u_2 = radius2 * pow(q,2);
    double v = u_2 * arms / (3.0 * arms - 2.0);

    double term1 = v - 1.0 + exp(-v);
    double term2 = ((arms - 1.0)/2.0)* pow((1.0 - exp(-v)),2.0);

    return (2.0 * (term1 + term2)) / (arms * pow(v,2.0));

}

double form_volume(void)
{
    return 1.0;
}

double Iq(double q, double radius2, double arms)
{
    return _mass_fractal_kernel(q, radius2, arms);
}

// Iqxy is never called since no orientation or magnetic parameters.
double Iqxy(double qx, double qy, double radius2, double arms)
{
    double q = sqrt(qx*qx + qy*qy);
    return _mass_fractal_kernel(q, radius2, arms);
}

