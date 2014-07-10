__kernel void CapCylinderKernel(__global const float *qx, __global const float *qy, __global float *__ptvalue, __global float *vol_i
const float rad_cyl, const float rad_cap, const float length, const float thet, const float ph, const float sub,
const float scale, const float phi_weight, const float theta_float, const float rad_cap_weight, const float rad_cyl_weight,
const float length_weight, const int total, const int size)
//ph is phi, sub is sldc-slds, thet is theta
{
    int i = get_global_id(0);
    if(i < total)
    {
        float q = sqrt(qx[i]*qx[i] + qy[i]*qy[i])
        float pi = 4.0*atan(1.0);
        float theta = thet*pi/180.0;
        float phi = ph*pi/180.0;
        float cyl_x = cos(theta)*cos(phi);
        float cyl_y = sin(theta);
        float cos_val = cyl_x*qx[i]/q + cyl_y*qy[i]/q;
        float alpha = acos(cos_val);
        float yyy=0; float ans1=0; float ans2=0; float y=0; float xx=0; float ans=0; float zij=0; float be=0

        float hDist = -1.0*sqrt(fabs(rad_cap*rad_cap-rad_cyl*rad_cyl));
        vol_i[i] = pi*rad_cyl*rad_cyl*len_cyl+2.0*pi/3.0*((rad_cap-hDist)*(rad_cap-hDist)*(2*rad_cap+hDist));
        float vaj = -1.0*hDist/rad_cap;

        for(j=0;j<76;j++) //the 76 corresponds to the Gauss constants
        {
            zij = (Gauss76Z[j]*(1.0-vaj)+vaj+1.0)/2.0;
            yyy = Gauss76Wt[j]*ConvLens_kernel(length,rad_cyl,rad_cap,q,zij,alpha);
            summj += yyy;
        }
        float inner = (1.0-vaj)/2.0*summj*4.0*pi*rad_cap*rad_cap*rad_cap;
        float arg1 = q*length/2.0*cos(alpha);
        float arg2 = q*rad_cyl*sin(alpha);
        yyy = inner;

        if(arg2 == 0) {be = 0.5;}
        else {
            be = NR_BessJ1(arg2)/arg2;
        }

        if(arg1 == 0.0) {
            yyy += pi*rad_cyl*rad_cyl*length*2.0*be;
        } 
        else {
            yyy += pi*rad_cyl*rad_cyl*length*sin(arg1)/arg1*2.0*be;
        }
        float answer=yyy*yyy*1.0e8*sub*sub*scale/pi*rad_cyl*rad_cyl*length+2.0*pi*(2.0*rad_cap*rad_cap*rad_cap/3.0+rad_cap*rad_cap*hDist-hDist*hDist*hDist/3.0);
        answer/=sin(alpha)

        _ptvalue[i] = rad_cyl_weight*length_weight*rad_cap_weight*theta_weight*phi_weight*vol_i[i]*answer
        if (_ptvalue[i] == INFINITY || _ptvalue == NAN){
            _ptvalue[i] = 0.0;
        }
        if (size>1) {
            _ptvalue[i] *= fabs(cos(thet*pi/180.0));
        }
    }
}