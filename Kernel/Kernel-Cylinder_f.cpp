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
#  include "NR_BessJ1.cpp"
#endif

#ifndef M_PI
#  define M_PI REAL(3.141592653589793)
#endif
#define RADIANS REAL(0.017453292519943295)
#define PI_E8 REAL(3.141592653589793e8)
real f(real qx, real qy, real sub, real rr, real h, real scale, real cyl_theta, real cyl_phi);
real f(real qx, real qy, real sub, real rr, real h, real scale, real cyl_theta, real cyl_phi)
{
    const real qq = sqrt(qx*qx+qy*qy);

    const real theta = cyl_theta*RADIANS;
    const real cos_val = cos(theta)*cos(cyl_phi*RADIANS)*(qx/qq) + sin(theta)*(qy/qq);

    real alpha = acos(cos_val);
    if(alpha == REAL(0.0)) {
        alpha = REAL(1.0e-26);
    }
    const real sin_alpha = sin(alpha);
    const real besarg = qq*rr*sin_alpha;
    const real siarg = qq*h/2*cos(alpha);

    real be, si;
    if (besarg == REAL(0.0)) {
        be = sin_alpha;
    } else {
        const real bj = NR_BessJ1(besarg)/besarg;
        be = REAL(4.0)*bj*bj*sin_alpha;
    }
    if (siarg == REAL(0.0)) {
        si = REAL(1.0);
    } else {
        si = sin(siarg)/siarg;
    }
    const real form = be*si*si/sin_alpha;
    return PI_E8*scale*form*sub*sub*rr*rr*rr*rr*h*h;
}

kernel void CylinderKernel(global const real *qx, global const real *qy, global real *result,
#ifdef USE_OPENCL
    global real *loops_g,
#else
    const int Nq,
#endif
    local real *loops, const real cutoff,
    const real scale, const real background,
    const real sub,
    const int Nradius, const int Nlength, const int Ntheta, const int Nphi)
{
#ifdef USE_OPENCL
    // copy loops info to local memory
    event_t e = async_work_group_copy(loops, loops_g, (Nradius+Nlength+Ntheta+Nphi)*2, 0);
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
        for (int ri=0; ri < Nradius; ri++) {
            const real rv = loops[2*ri];
            const real rw = loops[2*ri+1];
            for (int li=0; li < Nlength; li++) {
                const real lv = loops[2*(li+Nradius)];
                const real lw = loops[2*(li+Nradius)+1];
                for (int thi=0; thi < Ntheta; thi++) {
                    const real thv = loops[2*(thi+Nradius+Nlength)];
                    const real thw = loops[2*(thi+Nradius+Nlength)+1];
                    // #pragma unroll
                    for (int phi=0; phi < Nphi; phi++) {
                        const real phv = loops[2*(phi+Nradius+Nlength+Ntheta)];
                        const real phw = loops[2*(phi+Nradius+Nlength+Ntheta)+1];

                        const real weight = rw*lw*thw*phw;
                        //ret += qxi + qyi + sub + rv + lv + weight + thv + phv;
                        if (weight > cutoff) {
                            ret += f(qxi, qyi, sub, rv, lv, weight, thv, phv);
                            norm += weight;
                            vol += rw*lw*rv*rv*lv;
                            norm_vol += rw*lw;
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
