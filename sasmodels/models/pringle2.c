double form_volume(double radius, double thickness, double alpha, double beta);

double Iq(double q,
          double radius,
          double thickness,
          double alpha,
          double beta,
          double sld,
          double sld_solvent);


static void
integrand_bessel(
    double r,
    double alpha,
    double beta,
    double q_sin_psi,
    double q_cos_psi,
    double n, // Not the typical n for number of quadrature points. Here is the order of the Bessel function
    double *sumS,
    double *sumC,
)
{
    const double qrs = r*q_sin_psi;
    const double qrrc = r*r*q_cos_psi;

    double y = r * sas_JN(n, beta*qrrc) * sas_JN(2*n, qrs);
    double S, C;
    SINCOS(alpha*qrrc, S, C);
    SumS = y*S;
    SumC = y*C;
}
typedef void (*Integrand2_bessel)(double r, double alpha, double beta, double q_sin_psi, double q_cos_psi, double n, double* res1, double* res2);

void integrate2_bessel(
Integrand2_bessel f,
double a,
double b,
double q, // Required to determine the nodes for integration
double psi, // Required to determine the nodes for integration
double alpha,
double beta,
double q_sin_psi,
double q_sin_psi,
double n,
double* res1, double* res2){
    // Do what you need here
    // Determine the number of points for the Gauss quadrature
    int expo = (int)(eval_poly(log2m(max(limits[0][0],min(limits[0][1], q))), log2m(max(limits[1][0],min(limits[1][1], radius))), log2m(max(limits[2][0],min(limits[2][1], thickness))), log2m(max(limits[3][0],min(limits[3][1], length)))) + 1);
    int ng = (int)(pow(2, max(1, min(15, expo))));

    double *xg, *wg;
    get_gauss_points(n, &xg, &wg);

    // Perform the integration
    *res1 = 0;
    *res2 = 0;

    for (int i = 0; i < ng; i++){
        double t1, t2;
        f(a + (b - a) * 0.5 * (xg[i] + 1), alpha, beta, q_sin_psi, q_sin_psi, n, &t1, &t2);
        *res1 += t1 * wg[i];
        *res2 += t2 * wg[i];
    }

    *res1 *= (b - a) * 0.5;
    *res2 *= (b - a) * 0.5;
}


static
void _integrate_bessel(
    double radius,
    double q, // Required to determine the nodes for integration
    double psi, // Required to determine the nodes for integration
    double alpha,
    double beta,
    double q_sin_psi,
    double q_cos_psi,
    double n,
    double *Sn,
    double *Cn)
{

    // evaluate at Gauss points
    double sumS, sumC;		// initialize integral
    integrate2_bessel(integrand_bessel, 0, radius, q, psi, alpha, beta, q_sin_psi, q_cos_psi, n, &sumS, &sumC);

    *Sn = sumS / (radius*radius);
    *Cn = sumC / (radius*radius);
}

static
double _sum_bessel_orders(
    double radius,
    double q, // Required to determine the nodes for integration
    double psi, // Required to determine the nodes for integration
    double alpha,
    double beta,
    double q_sin_psi,
    double q_cos_psi)
{
    //calculate sum term from n = -3 to 3
    //Note 1:
    //    S_n(-x) = (-1)^S_n(x)
    //    => S_n^2(-x) = S_n^2(x),
    //    => sum_-k^k Sk = S_0^2 + 2*sum_1^kSk^2
    //Note 2:
    //    better precision to sum terms from smaller to larger
    //    though it doesn't seem to make a difference in this case.
    double Sn, Cn, sum;
    sum = 0.0;
    for (int n=3; n>0; n--) {
      _integrate_bessel(radius, q, psi, alpha, beta, q_sin_psi, q_cos_psi, n, &Sn, &Cn);
      sum += 2.0*(Sn*Sn + Cn*Cn);
    }
    _integrate_bessel(radius, q, psi, alpha, beta, q_sin_psi, q_cos_psi, 0, &Sn, &Cn);
    sum += Sn*Sn+ Cn*Cn;
    return sum;
}

static
double _integrate_psi(
    double q,
    double radius,
    double thickness,
    double alpha,
    double beta)
{
    // translate gauss point z in [-1,1] to a point in [0, pi/2]
    const double zm = M_PI_4;
    const double zb = M_PI_4;

    double sum = 0.0;
    for (int i = 0; i < GAUSS_N; i++) {
        double psi = GAUSS_Z[i]*zm + zb;
        double sin_psi, cos_psi;
        SINCOS(psi, sin_psi, cos_psi);
        double bessel_term = _sum_bessel_orders(radius, alpha, beta, q*sin_psi, q*cos_psi);
        double sinc_term = square(sas_sinx_x(q * thickness * cos_psi / 2.0));
        double pringle_kernel = 4.0 * sin_psi * bessel_term * sinc_term;
        sum += GAUSS_W[i] * pringle_kernel;
    }

    return zm * sum;
}

double form_volume(double radius, double thickness, double alpha, double beta)
{
    return M_PI*radius*radius*thickness;
}

static double
radius_from_excluded_volume(double radius, double thickness)
{
    return 0.5*cbrt(0.75*radius*(2.0*radius*thickness + (radius + thickness)*(M_PI*radius + thickness)));
}

static double
radius_effective(int mode, double radius, double thickness, double alpha, double beta)
{
    switch (mode) {
    default:
    case 1: // equivalent cylinder excluded volume
        return radius_from_excluded_volume(radius, thickness);
    case 2: // equivalent volume sphere
        return cbrt(M_PI*radius*radius*thickness/M_4PI_3);
    case 3: // radius
        return radius;
    }
}

double Iq(
    double q,
    double radius,
    double thickness,
    double alpha,
    double beta,
    double sld,
    double sld_solvent)
{
    double form = _integrate_psi(q, radius, thickness, alpha, beta);
    double contrast = sld - sld_solvent;
    double volume = M_PI*radius*radius*thickness;
    return 1.0e-4*form * square(contrast * volume);
}
