

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs[] = {-6.416688969553514e-14, 0.34954807644022323, 0.26077635791014114, 0.34954807644030794, 0.26077635791007064, -0.00041044938817771123, -0.06333809081229033, 0.021214783074584298, 0.023749082593212678, -0.0018600333610918258, 0.02374908259320729, 0.025064039050265354, -0.00041044938818576415, -0.06333809081228554, -0.001860033361086631, -0.0006281062767819012, 0.0011239112267517868, 9.310990856820164e-05, 0.0005069329096702667, 0.0014994896819265128, -0.0009018569154464675, -0.0011507235633038966, 9.31099085680646e-05, -0.0009018569154467864, -0.00014872597476390823, -0.0005728041354766125, -0.00014872597476366111, 0.00034753441047965627, 0.0005069329096700671, -0.0011507235633039117, 0.00034753441047996695, -0.0006281062767815696, 0.0011239112267517584, 0.0014994896819263553, -0.0005728041354768759};
const double intercept = 0.7979939663528164;
const double limits[][2] = {{1.0, 1000000.0}, {1.0, 1000000.0}, {1.0, 1000000.0}, {1.0, 1000000.0}};


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



double eval_poly(double var0, double var1, double var2, double var3){
    return 1 * coeffs[0] + var0 * coeffs[1] + var1 * coeffs[2] + var2 * coeffs[3] + var3 * coeffs[4] + (var0 * var0) * coeffs[5] + var0 * var1 * coeffs[6] + var0 * var2 * coeffs[7] + var0 * var3 * coeffs[8] + (var1 * var1) * coeffs[9] + var1 * var2 * coeffs[10] + var1 * var3 * coeffs[11] + (var2 * var2) * coeffs[12] + var2 * var3 * coeffs[13] + (var3 * var3) * coeffs[14] + (var0 * var0 * var0) * coeffs[15] + (var0 * var0) * var1 * coeffs[16] + (var0 * var0) * var2 * coeffs[17] + (var0 * var0) * var3 * coeffs[18] + var0 * (var1 * var1) * coeffs[19] + var0 * var1 * var2 * coeffs[20] + var0 * var1 * var3 * coeffs[21] + var0 * (var2 * var2) * coeffs[22] + var0 * var2 * var3 * coeffs[23] + var0 * (var3 * var3) * coeffs[24] + (var1 * var1 * var1) * coeffs[25] + (var1 * var1) * var2 * coeffs[26] + (var1 * var1) * var3 * coeffs[27] + var1 * (var2 * var2) * coeffs[28] + var1 * var2 * var3 * coeffs[29] + var1 * (var3 * var3) * coeffs[30] + (var2 * var2 * var2) * coeffs[31] + (var2 * var2) * var3 * coeffs[32] + var2 * (var3 * var3) * coeffs[33] + (var3 * var3 * var3) * coeffs[34] + intercept;
}


typedef double (*Integrand2)(const double x, const double q, const double core_r, const double core_h, const double core_vd, const double shell_r, const double shell_h, const double shell_vd, double* res1, double* res2, int n, int i);

void integrate2(
    Integrand2 f,
    const double a,
    const double b,
    const double q,
    const double core_r,
    const double core_h,
    const double core_vd,
    const double shell_r,
    const double shell_h,
    const double shell_vd,
    double* res1,
    double* res2
    ){


    const double A = q*core_h;
    const double B = q*core_r;
    const double C = q*shell_h;
    const double D = q*shell_r;

    int expo = (int)(eval_poly(log2m(max(limits[0][0],min(limits[0][1], A))), log2m(max(limits[1][0],min(limits[1][1], B))), log2m(max(limits[2][0],min(limits[2][1], C))), log2m(max(limits[3][0],min(limits[3][1], D)))) + 1);
    int n = (int)(pow(2, max(1, min(15, expo))));

    double *xg, *wg;
    get_gauss_points(n, &xg, &wg);

    // Perform the integration
    *res1 = 0;
    *res2 = 0;

    for (int i = 0; i < n; i++){
        double t1, t2;
        f(a + (b - a) * 0.5 * (xg[i] + 1), q, core_r, core_h, core_vd, shell_r, shell_h, shell_vd, &t1, &t2, n, i);
        *res1 += t1 * wg[i];
        *res2 += t2 * wg[i];
    }

    *res1 *= (b - a) * 0.5;
    *res2 *= (b - a) * 0.5;


}