

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs[] = {0.0, 0.46761427257906757, 0.40545475410012266, 0.2654372551626233, 0.022758141720083266, 0.039037590982758386, 0.019254937359550358, 0.02856956012352786, -0.030895080799509637, 0.021214182194932513};
const double intercept = 2.091842591936018;
const double limits[][2] = {{5e-05, 0.5}, {1.0, 200000.0}, {1.0, 200000.0}};


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



double eval_poly(double var0, double var1, double var2){
    return 1 * coeffs[0] + var0 * coeffs[1] + var1 * coeffs[2] + var2 * coeffs[3] + (var0 * var0) * coeffs[4] + var0 * var1 * coeffs[5] + var0 * var2 * coeffs[6] + (var1 * var1) * coeffs[7] + var1 * var2 * coeffs[8] + (var2 * var2) * coeffs[9] + intercept;
}
    W

typedef double (*Integrand2)(const double x, const double q, const double radius_polar, const double radius_equatorial, double* res1, double* res2);


void integrate2(
    Integrand2 f,
    const double a,
    const double b,
    const double q,
    double radius_polar,
    double radius_equatorial,
    double* res1,
    double* res2
    ){


    const double A = radius_equatorial;
    const double B = radius_polar;

    int expo = (int)(eval_poly(log2m(max(limits[0][0],min(limits[0][1], q))), log2m(max(limits[1][0],min(limits[1][1], A))), log2m(max(limits[2][0],min(limits[2][1], B)))) + 1);
    int n = (int)(pow(2, max(1, min(15, expo))));

    double *xg, *wg;
    get_gauss_points(n, &xg, &wg);

    // Perform the integration
    *res1 = 0;
    *res2 = 0;

    for (int i = 0; i < n; i++){
        double t1, t2;
        f(a + (b - a) * 0.5 * (xg[i] + 1), q, radius_polar, radius_equatorial, &t1, &t2);
        *res1 += t1 * wg[i];
        *res2 += t2 * wg[i];
    }

    *res1 *= (b - a) * 0.5;
    *res2 *= (b - a) * 0.5;


}
