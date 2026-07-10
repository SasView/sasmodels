/*
 * Fractal structure factor S(q)
 * Normalization: S(0) = N
 *
 * Parameters:
 *   q  : scattering vector
 *   radius_effective : effective scatterer radius (half center-to-center distance)
 *   volfraction      : unused (required for P@S parameter order)
 *   fractal_dim          : fractal dimension
 *   n_aggreg            : number of particles in cluster
 */

static double fractal_sq_N(double q, double r, double D, double N)
{
    const double eps = 1e-12;

    /* Special cases */
    if (N <= 1.0) return 1.0;     /* no clustering */
    if (D <= 0.0) return 1.0;

    /* q ? 0 limit */
    if (fabs(q) < eps) {
        return N;
    }

    const double Dm1 = D - 1.0;

    /* Avoid singularity at D = 1 */
    if (fabs(Dm1) < 1e-8) {
        /* Limit D ? 1:
           S(q) = 1 + (N-1)/N * atan(q*r*N)/(q*r)
        */
        double x = q * r;
        double y = x * N;  /* since N^(1/D) ? N when D?1 */
        return 1.0 + (N-1.0)/N * atan(y)/(x);
    }

    /* General case */
    double x = q * r;
    double Nd = pow(N, 1.0/D);
    double y = x * Nd;   /* = q*r*N^(1/D) */

    double term =
        sin(Dm1 * atan(y))
        * pow(x, -D)
        * pow(1.0 + 1.0/(y*y), -0.5*Dm1);

    double prefactor = (N - 1.0) / (N * Dm1);

    return 1.0 + prefactor * term;
}

/* ========== Sasmodels interface ========== */
double Iq(double q, double radius_effective, double volfraction,
          double fractal_dim, double n_aggreg)
{
    (void)volfraction;
    return fractal_sq_N(q, radius_effective, fractal_dim, n_aggreg);
}
