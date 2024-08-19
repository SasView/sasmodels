
#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))



#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs[] = {-2.1085518419661072e-13, 0.17078751320914443, 0.5630988992176875, 0.15708250267740692, 0.004165189443545988, 0.033386503963037197, 0.03849774028610544, 0.04511178363252381, -0.14043467597938902, 0.05834636243673773, -0.002285416559067707, -0.0006445814099293149, 0.001736028627359233, 0.002809936076959517, -0.004322137238567255, 0.004840204354366744, -0.0006931353225303171, 0.0047166814640176, 0.003403919096964598, 5.326798162702834e-05, -0.0001155283372183806, -0.00013283305963653576, -0.00012635952412694206, -0.0002867659736236193, 0.00017289544480935773, -0.0003317667249669333, -0.0002680998065252492, 6.777998111388389e-05, 9.441616281693074e-05, -0.00031374241719131003, -5.8686781060115016e-05, 0.0001561763262468802, -0.0005009287893979028, 0.00017183558245860756, -7.526875168870936e-05};
const double intercept = 1.3741159446473654;
const double limits[][2] = {{5e-05, 0.5}, {1.0, 1000000.0}, {1.0, 1000000.0}};


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
    return 1 * coeffs[0] + var0 * coeffs[1] + var1 * coeffs[2] + var2 * coeffs[3] + (var0 * var0) * coeffs[4] + var0 * var1 * coeffs[5] + var0 * var2 * coeffs[6] + (var1 * var1) * coeffs[7] + var1 * var2 * coeffs[8] + (var2 * var2) * coeffs[9] + (var0 * var0 * var0) * coeffs[10] + (var0 * var0) * var1 * coeffs[11] + (var0 * var0) * var2 * coeffs[12] + var0 * (var1 * var1) * coeffs[13] + var0 * var1 * var2 * coeffs[14] + var0 * (var2 * var2) * coeffs[15] + (var1 * var1 * var1) * coeffs[16] + (var1 * var1) * var2 * coeffs[17] + var1 * (var2 * var2) * coeffs[18] + (var2 * var2 * var2) * coeffs[19] + (var0 * var0 * var0 * var0) * coeffs[20] + (var0 * var0 * var0) * var1 * coeffs[21] + (var0 * var0 * var0) * var2 * coeffs[22] + (var0 * var0) * (var1 * var1) * coeffs[23] + (var0 * var0) * var1 * var2 * coeffs[24] + (var0 * var0) * (var2 * var2) * coeffs[25] + var0 * (var1 * var1 * var1) * coeffs[26] + var0 * (var1 * var1) * var2 * coeffs[27] + var0 * var1 * (var2 * var2) * coeffs[28] + var0 * (var2 * var2 * var2) * coeffs[29] + (var1 * var1 * var1 * var1) * coeffs[30] + (var1 * var1 * var1) * var2 * coeffs[31] + (var1 * var1) * (var2 * var2) * coeffs[32] + var1 * (var2 * var2 * var2) * coeffs[33] + (var2 * var2 * var2 * var2) * coeffs[34] + intercept;
}


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
