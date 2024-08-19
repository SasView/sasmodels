#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs[] = {-3.893195488343438e-12, 0.2860078408129299, 0.09260412944958968, 0.28600784081270336, 0.09260412944749757, 0.015648329140390122, -0.11303835998252526, 0.032763329409162364, 0.04496422015365223, 0.017366747940869176, 0.04496422015375682, 0.06398110215806642, 0.015648329140228026, -0.11303835998218541, 0.017366747941004214, -0.002851696212716831, 0.0026904251152494564, 0.0006960945224957508, 0.0025130255288071442, 0.005320486783323253, -0.0028151853962717598, -0.0029581648410473144, 0.0006960945224982001, -0.002815185396280127, -0.0019443264348413133, -0.002232441583256897, -0.0019443264348468384, -0.0002972461296777189, 0.0025130255287988843, -0.002958164841047644, -0.0002972461296772144, -0.0028516962127129175, 0.0026904251152431663, 0.005320486783321417, -0.0022324415832670735, 5.7687273517813276e-05, 5.607364128473735e-05, -5.853047399737021e-06, -5.78472650323565e-05, -0.00014389455742928534, 1.8943410556720996e-06, -3.951537307178565e-05, -6.318127120662036e-05, 8.15218583180774e-05, -3.488135381455546e-05, 1.8255544276500613e-05, -4.8984837239141044e-05, -0.00010140426101994116, 8.152185831821618e-05, 0.0001231270304931872, 0.00017003845884545205, -5.853047400052741e-06, 1.8943410558594498e-06, -4.8984837239130635e-05, 4.2961872975706694e-05, 3.5677365334614575e-05, 4.2961872975630366e-05, -3.739332735431394e-05, -3.4881353814578825e-05, 0.0001700384588454798, 5.412260598662219e-05, -5.78472650323552e-05, -3.951537307141095e-05, -0.00010140426101989952, -3.739332735426884e-05, 5.768727351808389e-05, 5.607364128473041e-05, -0.00014389455742937554, 1.8255544276233465e-05, 3.567736533493723e-05};
const double intercept = 1.767810299502674;
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
    return 1 * coeffs[0] + var0 * coeffs[1] + var1 * coeffs[2] + var2 * coeffs[3] + var3 * coeffs[4] + (var0 * var0) * coeffs[5] + var0 * var1 * coeffs[6] + var0 * var2 * coeffs[7] + var0 * var3 * coeffs[8] + (var1 * var1) * coeffs[9] + var1 * var2 * coeffs[10] + var1 * var3 * coeffs[11] + (var2 * var2) * coeffs[12] + var2 * var3 * coeffs[13] + (var3 * var3) * coeffs[14] + (var0 * var0 * var0) * coeffs[15] + (var0 * var0) * var1 * coeffs[16] + (var0 * var0) * var2 * coeffs[17] + (var0 * var0) * var3 * coeffs[18] + var0 * (var1 * var1) * coeffs[19] + var0 * var1 * var2 * coeffs[20] + var0 * var1 * var3 * coeffs[21] + var0 * (var2 * var2) * coeffs[22] + var0 * var2 * var3 * coeffs[23] + var0 * (var3 * var3) * coeffs[24] + (var1 * var1 * var1) * coeffs[25] + (var1 * var1) * var2 * coeffs[26] + (var1 * var1) * var3 * coeffs[27] + var1 * (var2 * var2) * coeffs[28] + var1 * var2 * var3 * coeffs[29] + var1 * (var3 * var3) * coeffs[30] + (var2 * var2 * var2) * coeffs[31] + (var2 * var2) * var3 * coeffs[32] + var2 * (var3 * var3) * coeffs[33] + (var3 * var3 * var3) * coeffs[34] + (var0 * var0 * var0 * var0) * coeffs[35] + (var0 * var0 * var0) * var1 * coeffs[36] + (var0 * var0 * var0) * var2 * coeffs[37] + (var0 * var0 * var0) * var3 * coeffs[38] + (var0 * var0) * (var1 * var1) * coeffs[39] + (var0 * var0) * var1 * var2 * coeffs[40] + (var0 * var0) * var1 * var3 * coeffs[41] + (var0 * var0) * (var2 * var2) * coeffs[42] + (var0 * var0) * var2 * var3 * coeffs[43] + (var0 * var0) * (var3 * var3) * coeffs[44] + var0 * (var1 * var1 * var1) * coeffs[45] + var0 * (var1 * var1) * var2 * coeffs[46] + var0 * (var1 * var1) * var3 * coeffs[47] + var0 * var1 * (var2 * var2) * coeffs[48] + var0 * var1 * var2 * var3 * coeffs[49] + var0 * var1 * (var3 * var3) * coeffs[50] + var0 * (var2 * var2 * var2) * coeffs[51] + var0 * (var2 * var2) * var3 * coeffs[52] + var0 * var2 * (var3 * var3) * coeffs[53] + var0 * (var3 * var3 * var3) * coeffs[54] + (var1 * var1 * var1 * var1) * coeffs[55] + (var1 * var1 * var1) * var2 * coeffs[56] + (var1 * var1 * var1) * var3 * coeffs[57] + (var1 * var1) * (var2 * var2) * coeffs[58] + (var1 * var1) * var2 * var3 * coeffs[59] + (var1 * var1) * (var3 * var3) * coeffs[60] + var1 * (var2 * var2 * var2) * coeffs[61] + var1 * (var2 * var2) * var3 * coeffs[62] + var1 * var2 * (var3 * var3) * coeffs[63] + var1 * (var3 * var3 * var3) * coeffs[64] + (var2 * var2 * var2 * var2) * coeffs[65] + (var2 * var2 * var2) * var3 * coeffs[66] + (var2 * var2) * (var3 * var3) * coeffs[67] + var2 * (var3 * var3 * var3) * coeffs[68] + (var3 * var3 * var3 * var3) * coeffs[69] + intercept;
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