#ifdef __OPENCL_VERSION__
# define USE_OPENCL
#elif defined(_OPENMP)
# define USE_OPENMP
#endif

// If opencl is not available, then we are compiling a C function
// Note: if using a C++ compiler, then define kernel as extern "C"
#ifndef USE_OPENCL
#  ifdef __cplusplus
      #include <cstdio>
      #include <cmath>
      using namespace std;
      #if defined(_MSC_VER)
         #include <limits>
         #include <float.h>
         #define kernel extern "C" __declspec( dllexport )
         static double trunc(double x) { return x>=0?floor(x):-floor(-x); }
         static double fmin(double x, double y) { return x>y ? y : x; }
         static double fmax(double x, double y) { return x<y ? y : x; }
         static double isnan(double x) { return _isnan(x); }
         #define NAN (std::numeric_limits<double>::quiet_NaN()) // non-signalling NaN
         static double cephes_expm1(double x) {
            // Adapted from the cephes math library.
            // Copyright 1984 - 1992 by Stephen L. Moshier
            if (x != x || x == 0.0) {
               return x; // NaN and +/- 0
            } else if (x < -0.5 || x > 0.5) {
               return exp(x) - 1.0;
            } else {
               const double xsq = x*x;
               const double p = (((
                  +1.2617719307481059087798E-4)*xsq
                  +3.0299440770744196129956E-2)*xsq
                  +9.9999999999999999991025E-1);
               const double q = ((((
                  +3.0019850513866445504159E-6)*xsq
                  +2.5244834034968410419224E-3)*xsq
                  +2.2726554820815502876593E-1)*xsq
                  +2.0000000000000000000897E0);
               double r = x * p;
               r =  r / (q - r);
               return r+r;
             }
         }
         #define expm1 cephes_expm1
         typedef __int32 int32_t
     #else
         #define kernel extern "C"
         #include <cstdint>
     #endif
     static void SINCOS(double angle, double &svar, double &cvar) { svar=sin(angle); cvar=cos(angle); }
#  else
     #include <inttypes.h>  // C99 guarantees that int32_t types is here
     #include <stdio.h>
     #include <tgmath.h> // C99 type-generic math, so sin(float) => sinf
     // MSVC doesn't support C99, so no need for dllexport on C99 branch
     #define kernel
     #define SINCOS(angle,svar,cvar) do {const double _t_=angle; svar=sin(_t_);cvar=cos(_t_);} while (0)
#  endif
#  define global
#  define local
#  define constant const
// OpenCL powr(a,b) = C99 pow(a,b), b >= 0
// OpenCL pown(a,b) = C99 pow(a,b), b integer
#  define powr(a,b) pow(a,b)
#  define pown(a,b) pow(a,b)
#else
   typedef int int32_t;
#  if defined(USE_SINCOS)
#    define SINCOS(angle,svar,cvar) svar=sincos(angle,&cvar)
#  else
#    define SINCOS(angle,svar,cvar) do {const double _t_=angle; svar=sin(_t_);cvar=cos(_t_);} while (0)
#  endif
#endif

// Standard mathematical constants:
//   M_E, M_LOG2E, M_LOG10E, M_LN2, M_LN10, M_PI, M_PI_2=pi/2, M_PI_4=pi/4,
//   M_1_PI=1/pi, M_2_PI=2/pi, M_2_SQRTPI=2/sqrt(pi), SQRT2, SQRT1_2=sqrt(1/2)
// OpenCL defines M_constant_F for float constants, and nothing if double
// is not enabled on the card, which is why these constants may be missing
#ifndef M_PI
#  define M_PI 3.141592653589793
#endif
#ifndef M_PI_2
#  define M_PI_2 1.570796326794897
#endif
#ifndef M_PI_4
#  define M_PI_4 0.7853981633974483
#endif
#ifndef M_E
#  define M_E 2.718281828459045091
#endif

// Non-standard function library
// pi/180, used for converting between degrees and radians
// 4/3 pi for computing sphere volumes
// square and cube for computing squares and cubes
#ifndef M_PI_180
#  define M_PI_180 0.017453292519943295
#endif
#ifndef M_4PI_3
#  define M_4PI_3 4.18879020478639
#endif
static inline double square(double x) { return x*x; }
static inline double cube(double x) { return x*x*x; }
static inline double sinc(double x) { return x==0 ? 1.0 : sin(x)/x; }

