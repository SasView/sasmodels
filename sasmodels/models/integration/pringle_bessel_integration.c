

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs[] = {0.0, 1.3527332332606141, 0.1965974372106335, 0.17643293684877298, 0.021315766164999146, 0.364918271604933, -0.0201932167888514, -0.0011367814267716349, -0.0009768819924034647, -0.000720778915708068, -0.016130514250697368, 0.01214417153082178, -0.025962668817754844, -0.0012186157556195718, -0.004861986637280383, 0.012436631293097088, -0.002520753145835999, -0.0034265133379312524, 0.0020146059577315274, -0.00164337014745737, -0.05254320987654153};
const double intercept = 2.93781055736494;
const double limits[][2] = {{1.0, 20000.0}, {5e-05, 1.0}, {5e-05, 1.0}, {5e-05, 1.0}, {2, 16}};


inline float log2m (float val)
{
   int * const    exp_ptr = (int *) (&val);
   int            x = *exp_ptr;
   const int      log_2 = ((x >> 23) & 255) - 128;
   x &= ~(255 << 23);
   x += 127 << 23;
   *exp_ptr = x;

   val = ((-1.0f/3) * val + 2) * val - 2.0f/3;   // (1)

   return (val + log_2);
}



double eval_poly(double var0, double var1, double var2, double var3, double var4){
    return 1 * coeffs[0] + var0 * coeffs[1] + var1 * coeffs[2] + var2 * coeffs[3] + var3 * coeffs[4] + var4 * coeffs[5] + (var0 * var0) * coeffs[6] + var0 * var1 * coeffs[7] + var0 * var2 * coeffs[8] + var0 * var3 * coeffs[9] + var0 * var4 * coeffs[10] + (var1 * var1) * coeffs[11] + var1 * var2 * coeffs[12] + var1 * var3 * coeffs[13] + var1 * var4 * coeffs[14] + (var2 * var2) * coeffs[15] + var2 * var3 * coeffs[16] + var2 * var4 * coeffs[17] + (var3 * var3) * coeffs[18] + var3 * var4 * coeffs[19] + (var4 * var4) * coeffs[20] + intercept;
}




typedef void (*Integrand2_bessel)(double r, double alpha, double beta, double q_sin_psi, double q_cos_psi, double n, double* res1, double* res2);

void integrate2_bessel(
Integrand2_bessel f,
double a,
double b,
double alpha,
double beta,
double q_sin_psi,
double q_sin_psi,
double n,
double* res1, double* res2){

    double A = q_cos_psi*alpha;
    double B = q_cos_psi*beta;
    double C = q_sin_psi;

    // Determine the number of points for the Gauss quadrature
    // don't log2 nj since it is discrete and linear
    int expo = (int)(eval_poly(log2m(max(limits[0][0],min(limits[0][1], R))), log2m(max(limits[1][0],min(limits[1][1], A))), log2m(max(limits[2][0],min(limits[2][1], B))), log2m(max(limits[3][0],min(limits[3][1], C))), max(limits[4][0],min(limits[4][1], nj))) + 1);
    int ng = (int)(pow(2, max(1, min(15, expo))));

    b = R;
    double *xg, *wg;
    get_gauss_points(ng, &xg, &wg);

    // Perform the integration
    *res1 = 0;
    *res2 = 0;

    for (int i = 0; i < ng; i++){
        double t1, t2;
        f(a + (b - a) * 0.5 * (xg[i] + 1), double alpha, double beta, double q_sin_psi, double q_cos_psi, double n, &t1, &t2, ng, i);
        *res1 += t1 * wg[i];
        *res2 += t2 * wg[i];
    }

    *res1 *= (b - a) * 0.5;
    *res2 *= (b - a) * 0.5;
  
    
}