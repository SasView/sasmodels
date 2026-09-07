/*
 * linear_aggregate
 *
 * Structure factor for a linear aggregate (Debye sum).
 *
 * Parameters:
 *   Q      - momentum transfer
 *   radius_effective - effective scatterer radius (half center-to-center distance)
 *   volfraction - unused (required for P@S parameter order)
 *   n_aggreg  - effective number of points (may be non-integer)
 *
 * Returns:
 *   linear_aggregate(Q)
 *
 * Notes:
 *   Exits with warning if n_aggreg > 100
 */
static double sq_linN(double Q, int n_aggreg, double sep)
{
    int k;
    double SUM, SN;

    double d = sep;
    
    if (n_aggreg <= 1)
        return 1.0;

/*
    if (n_aggreg > 100) {
        fprintf(stderr,
                "WARNING (linear_aggregate): n_aggreg = %d exceeds maximum allowed value (100).\n",
                n_aggreg);
        exit(EXIT_FAILURE);
    }
*/

    SUM = 0.0;

    for (k = 1; k <= n_aggreg - 1; k++) {
        SUM += (double)(n_aggreg - k) * sas_sinx_x(Q * d * (double)k);
    }

    SN = 1.0 + 2.0 * SUM / (double)n_aggreg;

    return SN;
}

/*
 * linear_aggregate (non-integer N)
 *
 * Linear-chain structure factor with non-integer effective length
 * handled by linear interpolation.
 *
 * Parameters:
 *   Q   - momentum transfer
 *   radius_effective - effective scatterer radius (half center-to-center distance)
 *   volfraction - unused
 *   n_aggreg - real (non-integer) number of points
 *
 * Returns:
 *   linear_aggregate(Q)
 */
double Iq(double Q, double radius_effective, double volfraction, double n_aggreg)
{
    int N;
    double w;
    double intensity;

    (void)volfraction;

    if (n_aggreg <= 1.0)
        return 1.0;

    N = (int)n_aggreg;
    w = n_aggreg - (double)N;

    double sep = 2.0 * radius_effective;

    intensity =
        (1.0 - w) * sq_linN(Q, N, sep)
      +        w  * sq_linN(Q, N + 1, sep);

    return intensity;
}