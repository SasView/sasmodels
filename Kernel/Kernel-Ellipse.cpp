__kernel void EllipsoidKernel(const real radius_a_weight, const real radius_b_weight,
const real axis_theta_weight,
const real axis_phi_weight, const real scale, const real radius_a, const real radius_b,
 const real sub, const real axis_theta, const real axis_phi, __global const real *qx,
__global const real *qy, __global real *_ptvalue, const int length, const int size)
{
     int i = get_global_id(0);
     if(i < length){
         real ret = 0;
         real q = sqrt(qx[i]*qx[i] + qy[i]*qy[i]);
         real pi = 4.0*atan(1.0);
         real theta = axis_theta*pi/180.0;
         real cos_val = cos(theta)*cos(axis_phi*pi/180.0)*(qx[i]/q) + sin(theta)*(qy[i]/q);

         real arg = q*radius_b*sqrt(1.0+(cos_val*cos_val*(((radius_a*radius_a/(radius_b*radius_b))-1.0))));
         if(arg == 0.0){
             ret = 1.0/3.0;
         }
         else{
             ret = (sin(arg)-arg*cos(arg))/(arg*arg*arg);
         }
         ret*=ret*9.0*sub*sub*4.0/3.0*acos(-1.0)*radius_b*radius_b*radius_a*scale*(1.0e8);

         _ptvalue[i] += radius_a_weight*radius_b_weight*axis_theta_weight*radius_a*axis_phi_weight*ret*pow(radius_b, 2);
         //if(size > 1){
          //  _ptvalue[i] *= fabs(cos(axis_theta*pi/180.0));
         //}
     }
}
