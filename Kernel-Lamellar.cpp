#ifndef real
# define real float
#endif

__kernel void LamellarKernel(__global const real *qx, global const real *qy, __global real *ret, const real bi_thick,
 const real scale, const real sub, const real background, const int length)
{
    int i = get_global_id(0);
    if(i < length)
    {
        real q = sqrt(qx[i]*qx[i]+qy[i]*qy[i]);
        real pi = 4.0*atan(1.0);
        real Pq = 2.0*sub*(sub/q)/q*(1.0-cos(q*bi_thick));
        ret[i] = 2.0*pi*scale*Pq/(q*q)/bi_thick*1.0e8;
        ret[i] += background;
    }
}