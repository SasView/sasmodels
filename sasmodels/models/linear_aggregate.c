#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * sinx(x) = sin(x)/x with safe x -> 0 limit
 */
static inline double sinx(double x)
{
    if (fabs(x) < 1e-8)
        return 1.0;
    return sin(x) / x;
}

/*
 * linear_aggregate
 *
 * Structure factor for a linear aggregate (Debye sum).
 *
 * Parameters:
 *   Q      - momentum transfer
 *   radius_effective - effective scatterer radius (half center-to-center distance)
 *   volfraction - unused (required for P@S parameter order)
 *   N_agg  - effective number of points (may be non-integer)
 *
 * Returns:
 *   linear_aggregate(Q)
 *
 * Notes:
 *   Exits with warning if N_agg > 100
 */
static double sq_linN(double Q, int N_agg, double sep)
{
    int k;
    double SUM, SN;

    double d = sep;
    
    if (N_agg <= 1)
        return 1.0;

    if (N_agg > 100) {
        fprintf(stderr,
                "WARNING (linear_aggregate): N_agg = %d exceeds maximum allowed value (100).\n",
                N_agg);
        exit(EXIT_FAILURE);
    }

    SUM = 0.0;

    for (k = 1; k <= N_agg - 1; k++) {
        SUM += (double)(N_agg - k) * sinx(Q * d * (double)k);
    }

    SN = 1.0 + 2.0 * SUM / (double)N_agg;

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
 *   N_agg - real (non-integer) number of points
 *
 * Returns:
 *   linear_aggregate(Q)
 */
double Iq(double Q, double radius_effective, double volfraction, double N_agg)
{
    int N;
    double w;
    double intensity;

    (void)volfraction;

    if (N_agg <= 1.0)
        return 1.0;

    N = (int)N_agg;
    w = N_agg - (double)N;

    double sep = 2.0 * radius_effective;

    intensity =
        (1.0 - w) * sq_linN(Q, N, sep)
      +        w  * sq_linN(Q, N + 1, sep);

    return intensity;
}