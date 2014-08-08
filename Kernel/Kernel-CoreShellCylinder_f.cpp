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

real f(real qx, real qy, real subone, real subtwo, real radius, real length, real thickness, real weight, real axis_theta, real axis_phi);
real f(real qx, real qy, real subone, real subtwo, real radius, real length, real thickness, real weight, real axis_theta, real axis_phi)
{
    const real q = sqrt(qx*qx+qy*qy);
    const real theta = axis_theta*RADIANS;
    real alpha = acos(cos(theta)*cos(axis_phi*RADIANS)*qx/q + sin(theta)*qy/q);

    if (alpha == 0.0){
    alpha = 1.0e-26;
    }

    real si1=0; real si2=0; real be1=0; real be2=0;

    const real vol2 = M_PI*(radius+thickness)*(radius+thickness)*(length+2.0*thickness);

    const real besarg1 = q*radius*sin(alpha);
    const real besarg2 = q*(radius+thickness)*sin(alpha);
    const real sinarg1 = q*length/2.0*cos(alpha);
    const real sinarg2 = q*(length/2.0+thickness)*cos(alpha);

    if (besarg1 == 0.0){be1 = 0.5;}
    else{be1 = NR_BessJ1(besarg1)/besarg1;}

    if (besarg2 == 0.0){be2 = 0.5;}
    else{be2 = NR_BessJ1(besarg2)/besarg2;}

    if (sinarg1 == 0.0){si1 = 1.0;}
    else{si1 = sin(sinarg1)/sinarg1;}

    if (sinarg2 == 0.0){si2 = 1.0;}
    else{si2 = sin(sinarg2)/sinarg2;}

    const real tt = 2.0*vol2*(subone)*si2*be2+2.0*(M_PI*radius*radius*(length))*(subtwo)*si1*be1;

    real answer = (tt*tt)*sin(alpha)/fabs(sin(alpha));
    //answer *= answer/(PI_E8*(radius+thickness)*(radius+thickness)*(length+2.0*thickness));

    return weight*answer*answer/PI_E8;

}

kernel void CoreShellCylinderKernel(global const real *qx, global const real *qy, global real *result,
#ifdef USE_OPENCL
    global real *loops_g,
#else
    const int Nq,
#endif
    local real *loops, const real cutoff,
    const real scale, const real background,
    const real subone, const real subtwo,
    const int Nradius, const int Nlength, const int Nthick, const int Ntheta, const int Nphi)
{
#ifdef USE_OPENCL
    // copy loops info to local memory
    event_t e = async_work_group_copy(loops, loops_g, (Nradius+Nlength+Nthick+Ntheta+Nphi)*2, 0);
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
                for (int ti=0; ti < Nthick; ti++) {
                    const real tv = loops[2*(ti+Nradius+Nlength)]; //////////////////////////
                    const real tw = loops[2*(ti+Nradius+Nlength)+1];
                    for (int thi=0; thi < Ntheta; thi++) {
                        const real thv = loops[2*(thi+Nradius+Nlength+Nthick)];
                        const real thw = loops[2*(thi+Nradius+Nlength+Nthick)+1];
                        // #pragma unroll
                        for (int phi=0; phi < Nphi; phi++) {
                            const real phv = loops[2*(phi+Nradius+Nlength+Ntheta+Nthick)];
                            const real phw = loops[2*(phi+Nradius+Nlength+Ntheta+Nthick)+1];

                            const real weight = rw*lw*tw*thw*phw;
                            //ret += qxi + qyi + sub + rv + lv + weight + thv + phv;
                            if (weight > cutoff) {
                                ret += f(qxi, qyi, subone, subtwo, rv, lv, tv, weight, thv, phv);
                                norm += weight;
                                vol += rw*lw*tw*(rv+tv)*(rv+tv)*(lv+2.0*tv);
                                norm_vol += rw*lw*tw;
                            }
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
