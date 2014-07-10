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

        for(j=0;j<76;j++) //the value 76 was pre-set in the kernel...
        {
            zij = ( Gauss76Z[j]*(1.0-vaj) + vaj + 1.0 )/2.0;    //the "t" dummy
            yyy = Gauss76Wt[j]*ConvLens_kernel(dp,q,zij,alpha);    //uses the same Kernel as the Dumbbell, here L>0
            summj += yyy;
        }
        float inner = (1.0-vaj)/2.0*summj*4.0*pi*rad_cap*rad_cap*rad_cap;
        float arg1 = q*length/2.0*cos(alpha);
        float arg2 = q*rad_cyl*sin(alpha);
        yyy = inner;

        if(arg2 == 0) {be = 0.5;}
        else {
            float ax=fabs(arg2);
            if ((ax < 8.0) {
                y=arg2*arg2;
                ans1=arg2*(72362614232.0+y*(-7895059235.0+y*(242396853.1+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
                ans2=144725228442.0+y*(2300535178.0+y*(18583304.74+y*(99447.43394+y*(376.9991397+y*1.0))));
                ans=ans1/ans2;
            } 
            else {
                y=64.0/(ax*ax);
                xx=ax-2.356194491;
                ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))));
                ans2=0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
                ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-(8.0/ax)*sin(xx)*ans2);
                if (arg2 < 0.0) {ans *= -1;}
	        }
	        be = ans/arg2
        }

        if(arg1 == 0.0) {   //limiting value of sinc(0) is 1; sinc is not defined in math.h
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