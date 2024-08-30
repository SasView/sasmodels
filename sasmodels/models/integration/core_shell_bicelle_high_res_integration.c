

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs[] = {-5.147118424772234e-14, 0.24061641424503566, 0.14183796960920794, 0.3784966142102907, 0.2187241624334997, -0.009153119723206813, -0.03420886828128656, 0.03282758067061214, 0.002986364139150438, 0.003913792714856429, 0.004704578711861654, 0.027459955873496967, -0.001172682984519246, -0.04372310400765254, -0.0030423394188799696, -0.000424292826953907, 0.0016576272302672125, -0.00020037873509965074, 0.00015632101463564408, 0.00047834883843919803, -0.0024648499273232307, -0.00027666244702073765, 0.00046439348603447917, 0.0002805227700425809, -0.0001478779330928415, -0.0008451160952173944, 0.0005768727908205622, 0.0008446148360184389, 0.001380770959436667, -0.0021160179626081043, -3.66266886512954e-05, -0.0009274240079766123, 0.00045763982083938797, 0.0015699469218321573, -0.000469573471707411};
const double intercept = 1.8984312053702386;
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
    
    
