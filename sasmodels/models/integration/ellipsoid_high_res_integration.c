

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs[] = {-3.4642539074800607e-14, 0.3684780758596396, 0.35752237537072834, 0.2692683938811623, 0.044966067529780614, 0.07629283209630257, 0.05302119319693169, 0.07659139238021585, -0.0986895674230704, 0.06084377611678292, 0.0008423623756695666, 0.0016806597353984942, 0.0019175048498948116, 0.002345110027970389, -0.002630061858112548, 0.0015588577445786694, -0.0024725461924337442, 0.0036659961517777106, 0.003283121229255571, -0.002481660505434884, -3.654880966457715e-05, -0.00025646737934223327, -6.400505516347556e-05, -0.0005405758397514663, 8.476260898502353e-05, -0.00033919909231283396, -0.00036626624065601776, -8.881020770035688e-05, 7.550570094533425e-05, -0.00027835434588737253, -1.0021005859930985e-05, 0.00010039314607481115, -0.0005181254455630624, 0.00018395995151758494, -1.2261849453824331e-05};
const double intercept = 1.2026407426875778;
const double limits[][2] = {{5e-05, 0.5}, {1.0, 200000.0}, {1.0, 200000.0}};


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
