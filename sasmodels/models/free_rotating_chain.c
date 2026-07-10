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
 * Sasmodels interface:
 *   Q   - momentum transfer
 *   radius_effective - effective scatterer radius (half center-to-center distance)
 *   volfraction - unused (required for P@S parameter order)
 *   N_agg - real number of rotating points (RN)
 */

double Iq(double Q, double radius_effective, double volfraction, double N_agg)
{
    int N_SPH, N_SPH1;
    double W;
    double ARG;
    double SN, SN1;
    double x;
    double intensity;

    double RN = N_agg;
    double d = 2.0 * radius_effective;

    (void)volfraction;

    /* integer and fractional parts of RN */
    N_SPH = (int)fabs(RN);
    W = fabs(RN) - (double)N_SPH;

    /* scaled argument Q*d */
    x = Q * d;

    /* sinc kernel */
    ARG = sas_sinx_x(x);

    /* FOR N */
    SN =  ((double)N_SPH) / (1.0 - ARG)
        - ((double)N_SPH) / 2.0
        - (1.0 - pow(ARG, N_SPH)) / square(1.0 - ARG) * ARG;

    SN = SN * 2.0 / (double)N_SPH;

    /* FOR N+1 */
    N_SPH1 = N_SPH + 1;

    SN1 = ((double)N_SPH1) / (1.0 - ARG)
        - ((double)N_SPH1) / 2.0
        - (1.0 - pow(ARG, N_SPH1)) / square(1.0 - ARG) * ARG;

    SN1 = SN1 * 2.0 / (double)N_SPH1;

    /* linear interpolation */
    intensity = (1.0 - W) * SN + W * SN1;

    return intensity;
}
