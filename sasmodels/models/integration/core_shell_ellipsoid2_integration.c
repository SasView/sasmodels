##include <stdio.h>


#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs[] = {0.0, 0.3150863379101902, 0.16613470175179204, 0.16613470175177467, 0.16613470175177394, 0.1661347017517684, 0.024856896057016388, 0.013270720644418109, 0.013270720644418034, 0.013270720644418244, 0.013270720644418057, 0.0035682981688967726, -0.022708438650798857, 0.012140265072473113, 0.012140265072473217, 0.0035682981688967158, 0.012140265072473257, 0.012140265072473266, 0.00356829816889676, -0.02270843865079874, 0.003568298168896832};
const double intercept = 0.8658426289390544;
const double limits[][2] = {{5e-05, 1.0}, {1.0, 1000000.0}, {1.0, 1000000.0}, {1.0, 1000000.0}, {1.0, 1000000.0}};


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
typedef double (*Integrand2)(const double x, const double q, const double radius_equat_core, const double polar_core, const double equat_shell, const double polar_shell, const double sld_core_shell, const double sld_shell_solvent, double* res1, double* res2, int n, int i);

void integrate2(
    Integrand2 f,
    const double a,
    const double b,
    const double q,
    const double radius_equat_core,
    const double polar_core,
    const double equat_shell,
    const double polar_shell,
    const double sld_core_shell,
    const double sld_shell_solvent,
    double* res1,
    double* res2
    ){


    const double A = radius_equat_core;
    const double B = polar_core;
    const double C = equat_shell;
    const double D = polar_shell;

    int expo = (int)(eval_poly(log2m(max(limits[0][0],min(limits[0][1], q))), log2m(max(limits[1][0],min(limits[1][1], A))), log2m(max(limits[2][0],min(limits[2][1], B))), log2m(max(limits[3][0],min(limits[3][1], C))), log2m(max(limits[4][0],min(limits[4][1], D)))) + 1);
    int n = (int)(pow(2, max(1, min(15, expo))));

    printf("n: %d\n", n);
    double *xg, *wg;
    get_gauss_points(n, &xg, &wg);

    // Perform the integration
    *res1 = 0;
    *res2 = 0;

    for (int i = 0; i < n; i++){
        double t1, t2;
        f(a + (b - a) * 0.5 * (xg[i] + 1), q, radius_equat_core, polar_core, equat_shell, polar_shell, sld_core_shell, sld_shell_solvent, &t1, &t2, n, i);
        *res1 += t1 * wg[i];
        *res2 += t2 * wg[i];
    }

    *res1 *= (b - a) * 0.5;
    *res2 *= (b - a) * 0.5;


}
