static double
fractal_sq(double q, double radius, double fractal_dim, double cor_length)
{
    // Calculate S(q) using Eq(16) from Teixeira (1988), DOI:10.1107/S0021889888000263
    // with terms rearranged so that limiting cases are easier to calculate.
    // Values checked against 500 bit floating point using the original equations.
    double term;
    if (q == 0.) {
        const double D = fractal_dim;
        term = pow(cor_length/radius, fractal_dim)*sas_gamma(fractal_dim+1.);
    } else if (fractal_dim == 0.) {
        term = 1.;
    } else if (fractal_dim == 1.) {
        term = atan(q*cor_length)/(q*radius);
    } else {
        const double t1 = pow(cor_length/radius, fractal_dim) * sas_gamma(fractal_dim+1.);
        const double t2 = sin((fractal_dim-1.)*atan(q*cor_length))/((fractal_dim-1)*q*cor_length);
        const double t3 = pow(1. + square(q*cor_length), -0.5*(fractal_dim-1.));
        term = t1 * t2 * t3;
    }
    return 1. + term;
}
