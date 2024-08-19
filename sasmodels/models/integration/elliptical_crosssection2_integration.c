
#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs[] = {-2208656.4220807115, 1.2220226966918923, 1.222022692842924, -0.10483263803625356, -0.6573633769947241, -0.10483263755849768, 0.0014198507849557203, 0.06536147090322794, 0.06536147158300348, 0.001419848744428018, 0.001304493575577002, 0.0001298307681846426, -0.011258321675452099, 0.00012984836919656908, 0.0013044786290321308, -0.00011455810310351766, -6.925519949993017e-05, 0.00028979159934117496, 0.0002898800009856692, -6.903997298997322e-05, -0.00011463960318743721, 2.8120088190197734e-06, -2.025824931584097e-06, 9.488224724463201e-06, -2.2037519149864515e-05, 9.493463229925947e-06, -2.0401750259746443e-06, 2.8165364841514773e-06};
const double intercept = 2208657.625436874;
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