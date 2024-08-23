

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs[] = {-1011830.2761837941, 0.809365194309996, 0.8093651918196798, 0.07797095642576618, -0.5056102268507974, 0.07797095677032968, -0.020912119904878665, 0.0446210404003155, 0.04462104071135143, -0.020912120849020534, 0.0019935969148695875, 0.0010748702336419752, -0.008324594780580306, 0.00107487829698814, 0.0019935900678297144, -0.00010003061559408488, -9.034208509735795e-05, 0.00019405245227936685, 0.000194092950869083, -9.024348547466421e-05, -0.00010006795243483468, 2.119587133198042e-06, -1.5744239736381616e-06, 8.983869866918681e-06, -1.806957091138317e-05, 8.986269732363938e-06, -1.5809980426495507e-06, 2.1216613482755164e-06};
const double intercept = 1011831.369345377;
const double limits[][2] = {{1.0, 1000000.0}, {1.0, 1000000.0}};


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



double eval_poly(double var0, double var1){
    return 1 * coeffs[0] + var0 * coeffs[1] + var1 * coeffs[2] + (var0 * var0) * coeffs[3] + var0 * var1 * coeffs[4] + (var1 * var1) * coeffs[5] + (var0 * var0 * var0) * coeffs[6] + (var0 * var0) * var1 * coeffs[7] + var0 * (var1 * var1) * coeffs[8] + (var1 * var1 * var1) * coeffs[9] + (var0 * var0 * var0 * var0) * coeffs[10] + (var0 * var0 * var0) * var1 * coeffs[11] + (var0 * var0) * (var1 * var1) * coeffs[12] + var0 * (var1 * var1 * var1) * coeffs[13] + (var1 * var1 * var1 * var1) * coeffs[14] + (var0 * var0 * var0 * var0 * var0) * coeffs[15] + (var0 * var0 * var0 * var0) * var1 * coeffs[16] + (var0 * var0 * var0) * (var1 * var1) * coeffs[17] + (var0 * var0) * (var1 * var1 * var1) * coeffs[18] + var0 * (var1 * var1 * var1 * var1) * coeffs[19] + (var1 * var1 * var1 * var1 * var1) * coeffs[20] + (var0 * var0 * var0 * var0 * var0 * var0) * coeffs[21] + (var0 * var0 * var0 * var0 * var0) * var1 * coeffs[22] + (var0 * var0 * var0 * var0) * (var1 * var1) * coeffs[23] + (var0 * var0 * var0) * (var1 * var1 * var1) * coeffs[24] + (var0 * var0) * (var1 * var1 * var1 * var1) * coeffs[25] + var0 * (var1 * var1 * var1 * var1 * var1) * coeffs[26] + (var1 * var1 * var1 * var1 * var1 * var1) * coeffs[27] + intercept;
}


typedef void (*Integrand)(double x, double q, double a, double b, double* out, int n, int i);

void integrate(Integrand f, double al, double bu, double q, double a, double b, double* res){

    double a_q = a*q;
    double b_q = b*q;

    // Determine the number of points for the Gauss quadrature
    int expo = (int)(eval_poly(log2m(max(limits[0][0],min(limits[0][1], a_q))), log2m(max(limits[1][0],min(limits[1][1], b_q)))) + 1);
    int n = (int)(pow(2, max(1, min(15, expo))));

    double *xg, *wg;
    get_gauss_points(n, &xg, &wg);

    // Perform the integration
    *res = 0;
    for (int i = 0; i < n; i++){
        double temp;
        f(al + (bu - al) * 0.5 * (xg[i] + 1), q, a, b, &temp, n, i);
        *res += temp * wg[i];
    }
    *res *= (bu - al) * 0.5;


}