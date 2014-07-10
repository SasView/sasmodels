__kernel void CoreShellCylinderKernel(__global const float *qx, __global const float *qy, __global float *_ptvalue,
const float axis_theta, const float axis_phi, const float thickness, const float length, const float radius,
const float scale, const float radius_weight, const float length_weight, const float thickness_weight,
const float theta_weight, const float phi_weight, const float core_sld, const float shell_sld, const float solvent_sld,
const int size, const int total)
{
    int i = get_global_id(0);
    if(i < total)
    {
        float q = sqrt(qx[i]*qx[i]+qy[i]*qy[i]);
        float pi = 4.0*atan(1.0);
        float theta = axis_theta*pi/180.0;
        float phi = axis_phi*pi/180.0;
        float cyl_x = cos(theta)*cos(phi);
        float cyl_y = sin(theta);
        float cos_val = cyl_x*qx[i]/q + cyl_y*qx[i]/q;
        float alpha = acos(cos_val);

        if (alpha == 0.0){
        alpha = 1.0e-26;
        }

        float si1=0; float si2=0; float be1=0; float be2=0;

        float dr1 = core_sld-shell_sld;
        float dr2 = shell_sld-solvent_sld;
        float vol1 = pi*radius*radius*(length);
        float vol2 = pi*(radius+thickness)*(radius+thickness)*(length+2.0*thickness);

        float besarg1 = q*radius*sin(alpha);
        float besarg2 = q*(radius+thickness)*sin(alpha);
        float sinarg1 = q*length/2.0*cos(alpha);
        float sinarg2 = q*(length/2.0+thickness)*cos(alpha);


        if (besarg1 == 0.0){be1 = 0.5;}
        else{
            be1 = NR_BessJ1(besarg1)/besarg1;
        }
        if (besarg2 == 0.0){be2 = 0.5;}
        else{
            be2 = NR_BessJ1(besarg2)/besarg2;
        }
        if (sinarg1 == 0.0){
            si1 = 1.0;
        }
        else{
            si1 = sin(sinarg1)/sinarg1;
        }
        if (sinarg2 == 0.0){
            si2 = 1.0;
        }
        else{
            si2 = sin(sinarg2)/sinarg2;
        }

        float t1 = 2.0*vol1*dr1*si1*be1;
        float t2 = 2.0*vol2*dr2*si2*be2;

        float answer = ((t1+t2)*(t1+t2))*sin(alpha)/fabs(sin(alpha));
        float vol=pi*(radius+thickness)*(radius+thickness)*(length+2.0*thickness);
        answer = answer/vol*1.0e8*scale;

        _ptvalue[i] = radius_weight*length_weight*thickness_weight*theta_weight*phi_weight*answer;
        _ptvalue[i] *= pow(radius+thickness,2)*(length+2.0*thickness);

        if (size>1) {
        _ptvalue[i] *= fabs(cos(axis_theta*pi/180.0));
        }

    }
}



