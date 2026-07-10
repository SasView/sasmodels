/*
 * Fractal structure factor S(q)
 *
 * Discrete?chain baseline using d directly.
 *
 * High?q limit:
 *     S(q) ? 1 + 2*sin(q d)/(q d) + 2*(sin(q d)/(q d))^2
 *
 * Low?q crossover constant: C ? 4.2
 *
 * Normalization: S(0) = N
 *
 * Parameters:
 *    q : scattering vector
 *    d : center-to-center distance between scatterers (= 2 * radius_effective)
 *    D : fractal dimension
 *    N : number of particles
 */

static double fractal_sq_N(double q, double d, double D, double N)
{
    const double eps = 1e-12;
    const double C   = 5.2;

    /* Trivial cases */
    if (N <= 1.0) return 1.0;
    if (D <= 0.0) return 1.0;

    /* q ? 0 limit */
    if (fabs(q) < eps) {
        return N;
    }

    /* ---------- discrete baseline term ---------- */

    double xd = q * d;
    double A;

    if (fabs(xd) < 1e-8)
        A = 1.0;
    else
        A = sin(xd)/xd;

    double f = 1.0 - exp(-pow(xd / C, 4.0));

    double S_base =  1.0 + (2.0*A + 2.0*A*A) * f;

    /* ---------- fractal correction ---------- */

    const double Dm1 = D - 1.0;

    /* D ? 1 limit */
    if (fabs(Dm1) < 1e-8) {
        double x = q * (d/2.0);     /* effective radius equivalent */
        double y = x * N;
        return S_base + (N - 1.0)/N * atan(y)/x;
    }

    /* General fractal contribution */
    double x  = q * (d/2.0);         /* consistent with original scaling */
    double Nd = pow(N, 1.0/D);
    double y  = x * Nd;

    double term =
        sin(Dm1 * atan(y)) *
        pow(x, -D) *
        pow(1.0 + 1.0/(y*y), -0.5*Dm1);

    double prefactor = (N - 1.0)/(N * Dm1);

    return S_base + prefactor * term;
}

/* ---------------- sasmodels interface ---------------- */

double Iq(double q, double radius_effective, double volfraction,
          double fractal_dim, double n_aggreg)
{
    (void)volfraction;
    return fractal_sq_N(q, 2.0 * radius_effective, fractal_dim, n_aggreg);
}
