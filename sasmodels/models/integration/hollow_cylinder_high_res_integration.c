

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs[] = {-4.947090333125777e-15, 0.5650307788861153, 0.3166711615676249, 0.2436727663119887, 0.16859361052538338, 0.013352225246354545, 0.015881407078550572, 0.01157181752170385, 0.013617508843405879, 0.01672260228770812, -0.008413826135788088, -0.013091202622457651, 0.006695877709010346, -0.006464139231893942, 0.015094052566744256, -0.00023295814235393811, -0.0014140558404368388, -0.0003536003752557154, -0.0005301238650319541, -0.0009008642929335271, -0.0004290062892654165, -0.0007700685021726832, -0.0008940285206048428, -4.103480688405019e-05, -0.0008494607178137389, 0.00043294569883312446, -0.0005498761232470379, -0.0008806238587081575, -0.001457729566320206, 0.0026715221234369624, -0.0010629819248146963, 0.000622075299171881, -0.0005752697865157617, -0.0011685378827954482, 0.0005223208296159535};
const double intercept = 2.950898721346258;
const double limits[][2] = {{0.0005, 0.5}, {1.0, 1000000.0}, {1.0, 1000000.0}, {1.0, 1000000.0}};


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
    
    