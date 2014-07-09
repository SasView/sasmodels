__kernel void EllipsoidKernel(const float radius_a_weight, const float radius_b_weight, const float axis_theta_weight,
const float axis_phi_weight, const float scale, const float radius_a, const float radius_b, const float sub, const float background, const float axis_theta, const float axis_phi, __global const float *qx,
__global const float *qy, __global float *_ptvalue, const int length, const int size)
{
     int i = get_global_id(0);
     if(i < length){
         float ret = 0;
         float q = sqrt(qx[i]*qx[i] + qy[i]*qy[i]);
         float pi = 4.0*atan(1.0);
         float theta = axis_theta*pi/180.0;
         float h = axis_phi*pi/180.0;
         float cyl_x = cos(theta)*cos(h);
         float cyl_y = sin(theta);
         float cos_val = cyl_x*(qx[i]/q) + cyl_y*(qy[i]/q);

         float nu = radius_a/radius_b;
         float arg = q*radius_b*sqrt(1.0+(cos_val*cos_val*((nu*nu)-1.0)));
         if(arg == 0.0){
             ret = 1.0/3.0;
         }
         else{
             ret = (sin(arg)-arg*cos(arg))/(arg*arg*arg);
         }
         ret*=ret*9.0*sub*sub;
         ret*=(4.0/3.0*acos(-1.0)*radius_b*radius_b*radius_a)*scale*(1.0e8);
         ret+=background;

         _ptvalue[i] = radius_a_weight*radius_b_weight*axis_theta_weight*radius_a*axis_phi_weight*ret*pow(radius_b, 2);

         if(size > 1){
            _ptvalue[i] *= fabs(cos(axis_theta*pi/180.0));
         }
     }
}
