static double
fractal_sq(double q, double radius, double fractal_dim, double cor_length)
{
    // Calculate S(q) using Eq(16) from Teixeira (1988), DOI:10.1107/S0021889888000263
    // with terms rearranged so that limiting cases are easier to calculate.
    // In particular, multiply by (qξ)^D/(qξ)^D and convert DΓ(D-1) to Γ(D+1)/(D-1)
    // giving (ξ/r)^D Γ(D+1), sin((D-1)atan(qξ))/((D-1)qξ) and ((qξ)² + 1.)^-((D-1)/2)
    // Values checked against 500 bit floating point using the original equations.
    double term;
    if (q == 0.) {
        // lim q->0 = 1 for terms t2 and t3 below
        term = pow(cor_length/radius, fractal_dim)*sas_gamma(fractal_dim + 1.);
    } else if (fractal_dim == 0.) {
        // sin(atan(x)) = x/√(1+x²) so terms t2 and t3 below cancel
        term = 1.;
    } else if (fractal_dim == 1.) {
        // lim x->0 sin(a x)/x = a, so t1 -> ξ/r and t2 -> atan(qξ)/(qξ)
        term = atan(q*cor_length)/(q*radius);
    } else {
        const double Dm1 = fractal_dim - 1.;
        const double t1 = pow(cor_length/radius, fractal_dim) * sas_gamma(fractal_dim + 1.);
        const double t2 = sin(Dm1*atan(q*cor_length)) / (Dm1*q*cor_length);
        const double t3 = pow(square(q*cor_length) + 1., -0.5*Dm1);
        term = t1 * t2 * t3;
    }
    return 1. + term;
}
