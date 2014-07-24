__kernel void OneDCylKernel(__global const real *q, __global real *answer, const real sub, const real length,
const real radius, const real scale, const real size, const real uplim, const real lolim)
{
    int k = get_global_id(0);
    if(k < size)
    {
        real Pi = 4.0*atan(1.0);
        real up = uplim*Pi/180;
        real low = lolim*Pi/180;

        real summ = 0.0; real zi = 0;
        real halfheight = length/2.0;

        for(int i=0; i<76; i++)
        {
            zi = (Gauss76Z(i)*(up-low) + up + low )/2.0;
            summ += Gauss76Wt(i)*CylKernel(q[k], radius, halfheight, zi);
        }

        answer[k] = (up-low)/2.0*summ*sub*sub*1.0e8*scale*Pi*radius*radius*length;

    }
}


