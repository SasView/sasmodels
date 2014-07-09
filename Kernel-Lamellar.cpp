__kernel void LamellarKernel(__global const float *qx, global const float *qy, __global float *ret, const float bi_thick,
 const float scale, const float sub, const float background, const int length)
{
    int i = get_global_id(0);
    if(i < length)
    {
        float q = sqrt(qx[i]*qx[i]+qy[i]*qy[i]);
        float pi = 4.0*atan(1.0);
        float Pq = 2.0*sub*(sub/q)/q*(1.0-cos(q*bi_thick));
        ret[i] = 2.0*pi*scale*Pq/(q*q)/bi_thick*1.0e8;
        ret[i] += background;
    }
}