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
 *   N_agg  - number of points (integer)
 *   dist_points      - separation between points (fit parameter)
 *
 * Returns:
 *   linear_aggregate(Q)
 *
 * Notes:
 *   Exits with warning if N_agg > 100
 */
static double sq_linN(double Q, int N_agg, double dist_points)
{
    int k;
    double SUM, SN;

    double d = dist_points;
    
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
 *   an  - real (non-integer) number of points
 *   d   - separation between points (fit parameter)
 *
 * Returns:
 *   linear_aggregate(Q)
 */
double Iq(double Q, double an, double d)
{
    int N;
    double w;
    double intensity;

    if (an <= 1.0)
        return 1.0;

    N = (int)an;
    w = an - (double)N;

    intensity =
        (1.0 - w) * sq_linN(Q, N, d)
      +        w  * sq_linN(Q, N + 1, d);

    return intensity;
}