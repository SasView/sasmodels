__kernel void CapCylinderKernel(__global const real *qx, __global const real *qy, __global real *_ptvalue, __global real *vol_i,
const real rad_cyl, const real rad_cap, const real length, const real thet, const real ph, const real sub,
const real scale, const real phi_weight, const real theta_float, const real rad_cap_weight, const real rad_cyl_weight,
const real length_weight, const int total, const int size, __const real Gauss76Wt, __const real Gauss76Z)
//ph is phi, sub is sldc-slds, thet is theta
{
    int i = get_global_id(0);
    if(i < total)
    {
        real q = sqrt(qx[i]*qx[i] + qy[i]*qy[i]);
        real pi = 4.0*atan(1.0);
        real theta = thet*pi/180.0;
        real phi = ph*pi/180.0;
        real cyl_x = cos(theta)*cos(phi);
        real cyl_y = sin(theta);
        real cos_val = cyl_x*qx[i]/q + cyl_y*qy[i]/q;
        real alpha = acos(cos_val);
        real yyy=0; real ans1=0; real ans2=0; real y=0; real xx=0; real ans=0; real zij=0; real be=0; real summj=0;

        real hDist = -1.0*sqrt(fabs(rad_cap*rad_cap-rad_cyl*rad_cyl));
        vol_i[i] = pi*rad_cyl*rad_cyl*length+2.0*pi/3.0*((rad_cap-hDist)*(rad_cap-hDist)*(2*rad_cap+hDist));
        real vaj = -1.0*hDist/rad_cap;

        for(int j=0;j<76;j++) //the 76 corresponds to the Gauss constants
        {
            zij = (Gauss76Z[j]*(1.0-vaj)+vaj+1.0)/2.0;
            yyy = Gauss76Wt[j]*ConvLens_kernel(length,rad_cyl,rad_cap,q,zij,alpha);
            summj += yyy;
        }
        real inner = (1.0-vaj)/2.0*summj*4.0*pi*rad_cap*rad_cap*rad_cap;
        real arg1 = q*length/2.0*cos(alpha);
        real arg2 = q*rad_cyl*sin(alpha);
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
        real answer=yyy*yyy*1.0e8*sub*sub*scale/pi*rad_cyl*rad_cyl*length+2.0*pi*(2.0*rad_cap*rad_cap*rad_cap/3.0+rad_cap*rad_cap*hDist-hDist*hDist*hDist/3.0);
        answer/=sin(alpha);

        _ptvalue[i] = rad_cyl_weight*length_weight*rad_cap_weight*theta_weight*phi_weight*vol_i[i]*answer
        if (_ptvalue[i] == INFINITY || _ptvalue == NAN){
            _ptvalue[i] = 0.0;
        }
     //   if (size>1) {
   //         _ptvalue[i] *= fabs(cos(thet*pi/180.0));
     //   }
    }
}