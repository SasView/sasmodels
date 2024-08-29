

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double coeffs[] = {51484.937074370806, 0.4660998468515879, 0.29585870214499826, 0.1561930621920928, -0.14762617861970184, 0.23117554547634084, -0.027738185367707127, 0.002808724188318949, -0.0057509136696648495, -0.03219955697332339, 0.0024801121719956994, 0.0014628462506479437, -0.0016954314294066605, 0.0022577622331882733, 0.0018880514698128993, -0.0001166007879894213, -7.727642724140362e-05, 3.4507104384076497e-05, 2.164564280124001e-05, -0.00010757930961477653, -5.042408959986602e-05, 2.168721230877413e-06, 5.262656228854912e-07, 2.562857471063995e-06, -4.279565356410386e-06, 2.5001213676822953e-06, 1.144056424484674e-06, 4.92101268345646e-07};
const double intercept = -51483.939312803006;
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
    

typedef double (*Integrand2)(double x , double q, double radius, double length, double* res1, double* res2, int n, int i);


