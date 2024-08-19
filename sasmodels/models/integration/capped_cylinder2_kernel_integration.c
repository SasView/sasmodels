#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs_kernel[] = {2.3827141882623314e-12, -5.542108933090031, -0.13865531386017627, -0.007153824294720216, -0.03671653954251309, -3.364265512286929, -0.31156604528672893, 0.0008862763456574653, -0.11356619343440391, 0.01685917444484655, 0.0013624422589718112, 0.00776670802553248, 0.0013278081379822014, -0.0007746778783671558, 0.029924326570780025, -0.7064272739798216, -0.07356259477344758, 0.00013144540690817642, 0.00014832965193216504, -0.0016666900705947685, -5.059422439302949e-05, 0.012025686197544232, 5.808034420289607e-05, -0.00010821197569768886, 0.002969450041612194, 0.004209849064297836, -9.112413062332408e-05, -0.0031515188437080594, -8.49087997085475e-05, 6.478559195982787e-06, -0.003855351358646158, -7.719960183685522e-05, 3.0510921246973757e-05, 4.577653108223112e-05, 0.004858334650140711, -0.0479501815963137, -0.004949205901038506, 6.728775343078919e-06, 0.0020820025648621565, 0.0003112412505642398, -4.103806631211571e-06, 0.000998667824300005, 2.001292871978144e-06, -9.795963510306205e-06, 0.0007785234232365323, 0.00019359585561917747, 4.9099643600891696e-08, -0.00011846136916035532, -1.3747582465334363e-06, -2.19406700965874e-06, -7.751739187103479e-05, -3.366080844446895e-06, 5.616409441994152e-06, -4.124051390282102e-07, 0.00016883394537897667, -0.00019274414964742803, 2.7697875750187784e-06, 0.0002477012637321918, 7.04201697454554e-07, -6.621754142788205e-07, -0.0002569426622086035, 2.096654826422295e-06, 4.945116613974676e-07, -6.975960564759021e-07, 0.00029989986945988534, 1.3747850489131963e-06, -1.1777601893747658e-06, 9.457406415780412e-07, -1.5375887975066505e-06, -0.00024483991008139583};
const double intercept_kernel = -0.16997855559810482;
const double limits_kernel[][2] = {{0.01, 0.9999}, {1.0, 1000000.0}, {1.0, 1000000.0}, {1.0, 1000000.0}};


inline float log2m_kernel (float val)
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



double eval_poly_kernel(double var0, double var1, double var2, double var3){
    return 1 * coeffs_kernel[0] + var0 * coeffs_kernel[1] + var1 * coeffs_kernel[2] + var2 * coeffs_kernel[3] + var3 * coeffs_kernel[4] + (var0 * var0) * coeffs_kernel[5] + var0 * var1 * coeffs_kernel[6] + var0 * var2 * coeffs_kernel[7] + var0 * var3 * coeffs_kernel[8] + (var1 * var1) * coeffs_kernel[9] + var1 * var2 * coeffs_kernel[10] + var1 * var3 * coeffs_kernel[11] + (var2 * var2) * coeffs_kernel[12] + var2 * var3 * coeffs_kernel[13] + (var3 * var3) * coeffs_kernel[14] + (var0 * var0 * var0) * coeffs_kernel[15] + (var0 * var0) * var1 * coeffs_kernel[16] + (var0 * var0) * var2 * coeffs_kernel[17] + (var0 * var0) * var3 * coeffs_kernel[18] + var0 * (var1 * var1) * coeffs_kernel[19] + var0 * var1 * var2 * coeffs_kernel[20] + var0 * var1 * var3 * coeffs_kernel[21] + var0 * (var2 * var2) * coeffs_kernel[22] + var0 * var2 * var3 * coeffs_kernel[23] + var0 * (var3 * var3) * coeffs_kernel[24] + (var1 * var1 * var1) * coeffs_kernel[25] + (var1 * var1) * var2 * coeffs_kernel[26] + (var1 * var1) * var3 * coeffs_kernel[27] + var1 * (var2 * var2) * coeffs_kernel[28] + var1 * var2 * var3 * coeffs_kernel[29] + var1 * (var3 * var3) * coeffs_kernel[30] + (var2 * var2 * var2) * coeffs_kernel[31] + (var2 * var2) * var3 * coeffs_kernel[32] + var2 * (var3 * var3) * coeffs_kernel[33] + (var3 * var3 * var3) * coeffs_kernel[34] + (var0 * var0 * var0 * var0) * coeffs_kernel[35] + (var0 * var0 * var0) * var1 * coeffs_kernel[36] + (var0 * var0 * var0) * var2 * coeffs_kernel[37] + (var0 * var0 * var0) * var3 * coeffs_kernel[38] + (var0 * var0) * (var1 * var1) * coeffs_kernel[39] + (var0 * var0) * var1 * var2 * coeffs_kernel[40] + (var0 * var0) * var1 * var3 * coeffs_kernel[41] + (var0 * var0) * (var2 * var2) * coeffs_kernel[42] + (var0 * var0) * var2 * var3 * coeffs_kernel[43] + (var0 * var0) * (var3 * var3) * coeffs_kernel[44] + var0 * (var1 * var1 * var1) * coeffs_kernel[45] + var0 * (var1 * var1) * var2 * coeffs_kernel[46] + var0 * (var1 * var1) * var3 * coeffs_kernel[47] + var0 * var1 * (var2 * var2) * coeffs_kernel[48] + var0 * var1 * var2 * var3 * coeffs_kernel[49] + var0 * var1 * (var3 * var3) * coeffs_kernel[50] + var0 * (var2 * var2 * var2) * coeffs_kernel[51] + var0 * (var2 * var2) * var3 * coeffs_kernel[52] + var0 * var2 * (var3 * var3) * coeffs_kernel[53] + var0 * (var3 * var3 * var3) * coeffs_kernel[54] + (var1 * var1 * var1 * var1) * coeffs_kernel[55] + (var1 * var1 * var1) * var2 * coeffs_kernel[56] + (var1 * var1 * var1) * var3 * coeffs_kernel[57] + (var1 * var1) * (var2 * var2) * coeffs_kernel[58] + (var1 * var1) * var2 * var3 * coeffs_kernel[59] + (var1 * var1) * (var3 * var3) * coeffs_kernel[60] + var1 * (var2 * var2 * var2) * coeffs_kernel[61] + var1 * (var2 * var2) * var3 * coeffs_kernel[62] + var1 * var2 * (var3 * var3) * coeffs_kernel[63] + var1 * (var3 * var3 * var3) * coeffs_kernel[64] + (var2 * var2 * var2 * var2) * coeffs_kernel[65] + (var2 * var2 * var2) * var3 * coeffs_kernel[66] + (var2 * var2) * (var3 * var3) * coeffs_kernel[67] + var2 * (var3 * var3 * var3) * coeffs_kernel[68] + (var3 * var3 * var3 * var3) * coeffs_kernel[69] + intercept_kernel;
}


typedef void (*Integrand)(double x, double , double , double , double , double* out, int n, int i);

void integrate_kernel(Integrand f, double a, double bu, const double lower, const double m, const double b, const double qab_r, double* res){

    double A =  lower;
    double B =  m;
    double C =  b;
    double D =  qab_r;
    // Determine the number of points for the Gauss quadrature
    int expo = (int)(eval_poly_kernel(log2m_kernel(max(limits_kernel[0][0],min(limits_kernel[0][1], A))), log2m_kernel(max(limits_kernel[1][0],min(limits_kernel[1][1], B))), log2m_kernel(max(limits_kernel[2][0],min(limits_kernel[2][1], C))), log2m_kernel(max(limits_kernel[3][0],min(limits_kernel[3][1], D)))) + 1);
    int n = (int)(pow(2, max(1, min(15, expo))));

    double *xg, *wg;
    get_gauss_points(n, &xg, &wg);
    a = A;
    // Perform the integration
    *res = 0;
    for (int i = 0; i < n; i++){
        double temp;
        f(a + (bu - a) * 0.5 * (xg[i] + 1), lower, m, b, qab_r, &temp, n, i);
        *res += temp * wg[i];
    }
    *res *= (bu - a) * 0.5;


}

