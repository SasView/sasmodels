#ifdef __OPENCL_VERSION__
# define USE_OPENCL
#endif

// If opencl is not available
#ifndef USE_OPENCL
#  include <math.h>
#  define REAL(x) (x)
#  define real double
#  define global
#  define local
#  define kernel extern "C"
#endif

#ifndef M_PI
#  define M_PI REAL(3.141592653589793)
#endif
#define RADIANS REAL(0.017453292519943295)
#define PI_E8 REAL(3.141592653589793e8)

real f(real qx, real qy, real sub, real radius_a, real radius_b, real scale, real axis_theta, real axis_phi);
real f(real qx, real qy, real sub, real radius_a, real radius_b, real scale, real axis_theta, real axis_phi)
{
    real ret = 0;
    const real q = sqrt(qx*qx + qy*qy);
    const real theta = axis_theta*RADIANS;
    const real cos_val = cos(theta)*cos(axis_phi*RADIANS)*(qx/q) + sin(theta)*(qy/q);

    const real arg = q*radius_b*sqrt(1.0+(cos_val*cos_val*(((radius_a*radius_a/(radius_b*radius_b))-1.0))));
    if(arg == 0.0){
     ret = 1.0/3.0;
    }
    else{
     ret = (sin(arg)-arg*cos(arg))/(arg*arg*arg);
    }
    ret*=ret*9.0*sub*sub*4.0/3.0*acos(-1.0)*radius_b*radius_b*radius_a*(1.0e8);

    return radius_a*scale*ret*pow(radius_b, 2);
}

kernel void EllipsoidKernel(global const real *qx, global const real *qy, global real *result,
#ifdef USE_OPENCL
    global real *loops_g,
#else
    const int Nq,
#endif
    local real *loops, const real cutoff,
    const real scale, const real background,
    const real sub,
    const int Nradius_a, const int Nradius_b, const int Ntheta, const int Nphi)
{
#ifdef USE_OPENCL
    // copy loops info to local memory
    event_t e = async_work_group_copy(loops, loops_g, (Nradius_a+Nradius_b+Ntheta+Nphi)*2, 0);
    wait_group_events(1, &e);

    int i = get_global_id(0);
    int count = get_global_size(0);

    if(i < count)
#else
    #pragma omp parallel for
    for (int i=0; i<Nq; i++)
#endif
{
        const real qxi=qx[i];
        const real qyi=qy[i];
        real ret=REAL(0.0), norm=REAL(0.0), norm_vol=REAL(0.0), vol=REAL(0.0);
        for (int ai=0; ai < Nradius_a; ai++) {
            const real rav = loops[2*ai];
            const real raw = loops[2*ai+1];
            for (int bi=0; bi < Nradius_b; bi++) {
                const real rbv = loops[2*(bi+Nradius_a)];
                const real rbw = loops[2*(bi+Nradius_a)+1];
                for (int thi=0; thi < Ntheta; thi++) {
                    const real thv = loops[2*(thi+Nradius_a+Nradius_b)];
                    const real thw = loops[2*(thi+Nradius_a+Nradius_b)+1];
                    // #pragma unroll
                    for (int phi=0; phi < Nphi; phi++) {
                        const real phv = loops[2*(phi+Nradius_a+Nradius_b+Ntheta)];
                        const real phw = loops[2*(phi+Nradius_a+Nradius_b+Ntheta)+1];

                        const real weight = raw*rbw*thw*phw;
                        //ret += qxi + qyi + sub + rav + rbv + weight + thv + phv;
                        if (weight > cutoff) {
                            ret += f(qxi, qyi, sub, rav, rbv, weight, thv, phv);
                            norm += weight;
                            vol += raw*rbw*rav*rbv*rbv;
                            norm_vol += raw*rbw;
                        }
                    }
                }
            }
        }
        //if (Ntheta>1) norm = norm/(M_PI/2);
        if (vol != REAL(0.0) && norm_vol != REAL(0.0)) {
            ret *= norm_vol/vol;
        }
        result[i] = scale*ret/norm+background;
    }
}
