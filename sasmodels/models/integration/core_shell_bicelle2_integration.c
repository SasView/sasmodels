


#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs[] = {-4.796453249045279e-12, 0.4030560772171777, -0.037096218390153333, 0.4022721157683786, 0.16121081667516293, -0.01722429713602069, -0.06471335426176533, 0.0046134107854336436, 0.010559345304565204, 0.03636413474651837, -0.002571483690790081, 0.059631669768668204, 0.023160457282325755, -0.07221438734400962, -0.00022771152437634987, -0.0013756874809705832, 0.003862069346776106, 0.0017512540362280245, 0.001133153908017063, 0.001486626745724036, -0.004772883539478669, -0.00040279384486605063, 0.003415780803055544, 0.00029049939360159874, -0.0014204836343689664, -0.0034944600967188465, 0.0011093110854168199, 0.0010290816808170098, 0.003701678086434768, -0.0048993285955974755, -0.001068672386278162, -0.004105947103035079, 0.0006975742872367821, 0.004624441924683564, -0.0010275456263955287, 4.914013967583647e-05, 1.6821461058861426e-05, -7.588823912142717e-05, -4.202767111699188e-05, -0.00015574104268995578, 3.333662948545868e-05, 6.480005795644894e-06, -3.3533683203768386e-05, 6.556211938307799e-05, -2.1988889521339094e-05, 7.430781757472929e-05, 6.929566012464639e-05, -8.191099215584807e-05, 2.4473636770425133e-05, -2.261607000747512e-05, 9.306724372207253e-05, -7.838400365264581e-05, -1.8406296480616424e-05, -3.634833172332097e-05, 3.831878776174415e-05, 6.120917112181767e-05, 3.995449087473846e-05, -9.325498850894925e-05, -0.00016515760864193457, 8.772945202978333e-05, 0.00012771824391254183, 3.3771939425102093e-06, 6.282244456133623e-05, 3.994697097969846e-07, -8.178143243288893e-05, 8.327519466518746e-05, 6.0849628891246343e-05, -0.0001255204292774023, -6.50282841569913e-06, 2.6488564547363577e-05};
const double intercept = 1.9898114292958748;
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
    
    

typedef double (*Integrand2)(double x, double q, double radius, double thick_radius, double thick_face, double halflength, double sld_core, double sld_face, double sld_rim, double sld_solvent, double* res1, double* res2, int n, int i);



void integrate2(
    Integrand2 f,
    double a,
    double b,
    double q,
    double radius,
    double thick_radius,
    double thick_face,
    double halflength,
    double sld_core,
    double sld_face,
    double sld_rim,
    double sld_solvent,
    double* res1,
    double* res2
    ){


    const double A = q*halflength;
    const double B = q*radius;
    const double C = q*(halflength+thick_face);
    const double D = q*(radius+thick_radius);

    int expo = (int)(eval_poly(log2m(max(limits[0][0],min(limits[0][1], A))), log2m(max(limits[1][0],min(limits[1][1], B))), log2m(max(limits[2][0],min(limits[2][1], C))), log2m(max(limits[3][0],min(limits[3][1], D)))) + 1);
    int n = (int)(pow(2, max(1, min(15, expo))));

    double *xg, *wg;
    get_gauss_points(n, &xg, &wg);

    // Perform the integration
    *res1 = 0;
    *res2 = 0;

    for (int i = 0; i < n; i++){
        double t1, t2;
        f(a + (b - a) * 0.5 * (xg[i] + 1), q, radius, thick_radius, thick_face, halflength, sld_core, sld_face, sld_rim, sld_solvent, &t1, &t2, n, i);
        *res1 += t1 * wg[i];
        *res2 += t2 * wg[i];
    }

    *res1 *= (b - a) * 0.5;
    *res2 *= (b - a) * 0.5;


}
