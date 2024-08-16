

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs[] = {-581968.5121967047, 0.5548933974506951, 0.5968599873023919, -0.04381648544730543, -0.19127954929202098, -0.1497354873884559, 0.027267390093687477, 0.02096722822030451, 0.008220670778053996, 0.059824639961972785, -0.0037343880537537037, -0.0002554625787731472, -0.0029986600437200885, 0.0002884765374963119, -0.007189090257242044, 0.00019726025781617072, -2.2406852916695792e-05, 0.00011680153529615087, 5.857137481234467e-05, 8.987540259321886e-06, 0.000353335422321268, -3.631031312678612e-06, 1.4153674057770793e-07, 6.053141692607955e-07, -4.612287085203892e-06, 1.6905875104833434e-06, -1.1904277164298538e-06, -6.23197834831446e-06};
const double intercept = 581971.7934611484;
const double limits[][2] = {{1.0, 1000000.0}, {1.0, 1000000.0}};


        
double eval_poly(double var0, double var1){
    return 1 * coeffs[0] + var0 * coeffs[1] + var1 * coeffs[2] + (var0 * var0) * coeffs[3] + var0 * var1 * coeffs[4] + (var1 * var1) * coeffs[5] + (var0 * var0 * var0) * coeffs[6] + (var0 * var0) * var1 * coeffs[7] + var0 * (var1 * var1) * coeffs[8] + (var1 * var1 * var1) * coeffs[9] + (var0 * var0 * var0 * var0) * coeffs[10] + (var0 * var0 * var0) * var1 * coeffs[11] + (var0 * var0) * (var1 * var1) * coeffs[12] + var0 * (var1 * var1 * var1) * coeffs[13] + (var1 * var1 * var1 * var1) * coeffs[14] + (var0 * var0 * var0 * var0 * var0) * coeffs[15] + (var0 * var0 * var0 * var0) * var1 * coeffs[16] + (var0 * var0 * var0) * (var1 * var1) * coeffs[17] + (var0 * var0) * (var1 * var1 * var1) * coeffs[18] + var0 * (var1 * var1 * var1 * var1) * coeffs[19] + (var1 * var1 * var1 * var1 * var1) * coeffs[20] + (var0 * var0 * var0 * var0 * var0 * var0) * coeffs[21] + (var0 * var0 * var0 * var0 * var0) * var1 * coeffs[22] + (var0 * var0 * var0 * var0) * (var1 * var1) * coeffs[23] + (var0 * var0 * var0) * (var1 * var1 * var1) * coeffs[24] + (var0 * var0) * (var1 * var1 * var1 * var1) * coeffs[25] + var0 * (var1 * var1 * var1 * var1 * var1) * coeffs[26] + (var1 * var1 * var1 * var1 * var1 * var1) * coeffs[27] + intercept;
}
    

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

    
typedef void (*Integrand2) (double, double, double, double*, double*);

void integrate2(Integrand2 f, double a, double b, double A, double B, double *result1, double *result2){


    // Determine the number of points for the Gauss quadrature
    int expo = (int)(eval_poly(log2m(max(limits[0][0],min(limits[0][1], A))), log2m(max(limits[1][0],min(limits[1][1], B)))) + 1);
    int n = (int)(pow(2, max(1, min(15, expo))));

    double *xg, *wg;
    get_gauss_points(n, &xg, &wg);

    // Perform the integration
    *result1 = 0;
    *result2 = 0;

    for (int i = 0; i < n; i++){
        double t1, t2;
        f(a + (b - a) * 0.5 * (xg[i] + 1), A, B, &t1, &t2);
        *result1 += t1 * wg[i];
        *result2 += t2 * wg[i];
    }
    *result1 *= (b - a) * 0.5;
    *result2 *= (b - a) * 0.5;

}