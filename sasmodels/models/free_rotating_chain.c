#include <math.h>

/*
 * free_rotating_chain
 *
 * Computes SQ free rotating points:
 *
 * Burchard, W., & Kajiwara, K. (1970). The statistics of stiff chain molecules I.
 * The particle scattering factor.
 * Proceedings of the Royal Society of London. A. Mathematical and Physical Sciences,
 * 316(1525), 185-199.
 *
 * Inputs:
 *   Q   - momentum transfer
 *   RN  - real number of rotating points
 *   d   - separation between points (fit parameter)
 *
 * Output:
 *   Returns free_rotating_chain(Q, RN, d)
 *
 * Notes:
 *   This routine is a direct translation of the original Fortran code by Jan Skov Pedersen.
 *   All variables are treated as real (double precision).
 */

double Iq(double Q, double N_agg, double dist_points)
{
    int N_SPH, N_SPH1;
    double W;
    double ARG;
    double SN, SN1;
    double x;
    double intensity;

	double RN = N_agg;
	double d = dist_points;

    /* integer and fractional parts of RN */
    N_SPH = (int)fabs(RN);
    W = fabs(RN) - (double)N_SPH;

    /* scaled argument Q*d */
    x = Q * d;

    /* sinc kernel */
    ARG = sin(x) / x;

    /* FOR N */
    SN =  ((double)N_SPH) / (1.0 - ARG)
        - ((double)N_SPH) / 2.0
        - (1.0 - pow(ARG, N_SPH)) / pow(1.0 - ARG, 2.0) * ARG;

    SN = SN * 2.0 / (double)N_SPH;

    /* FOR N+1 */
    N_SPH1 = N_SPH + 1;

    SN1 = ((double)N_SPH1) / (1.0 - ARG)
        - ((double)N_SPH1) / 2.0
        - (1.0 - pow(ARG, N_SPH1)) / pow(1.0 - ARG, 2.0) * ARG;

    SN1 = SN1 * 2.0 / (double)N_SPH1;

    /* linear interpolation */
    intensity = (1.0 - W) * SN + W * SN1;

    return intensity;
}