__kernel void CylinderKernel(__global const float *qx, global const float *qy, __global float *_ptvalue, const float sub,
const float rr, const float h, const float scale, const float radius_weight, const float length_weight,
const float theta_weight, const float phi_weight, const float cyl_theta,
const float cyl_phi, const int count, const int size)
{
	// qq is the q-value for the calculation (1/A)
	// rr is the radius of the cylinder (A)
	// h is the -LENGTH of the cylinder = L (A)
    int i = get_global_id(0);

    if(i < count)
    {
        float qq = sqrt(qx[i]*qx[i]+qy[i]*qy[i]);

        float pi = 4.0*atan(1.0);
        float theta = cyl_theta*pi/180.0;
        float phi = cyl_phi*pi/180.0;

        float cyl_x = cos(theta)*cos(phi);
        float cyl_y = sin(theta);
        float cos_val = cyl_x*(qx[i]/qq) + cyl_y*(qy[i]/qq);

        float alpha = acos(cos_val);
        if(alpha == 0.0){
            alpha = 1.0e-26;
        }
        float besarg = qq*rr*sin(alpha);
        float siarg = qq*h/2*cos(alpha);

	    float xx=0.0; float y=0.0; float bj=0.0; float ans1=0.0; float ans2=0.0; float z=0.0; float answer=0.0;
	    float contrast=0.0; float form=0.0; float be=0.0; float si=0.0;

        float ax = fabs(besarg);

        if(ax < 8.0) {
            y=besarg*besarg;
            ans1=besarg*(72362614232.0+y*(-7895059235.0+y*(242396853.1+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
            ans2=144725228442.0+y*(2300535178.0+y*(18583304.74+y*(99447.43394+y*(376.9991397+y*1.0))));
            bj=ans1/ans2;
        }
        else{
            z=8.0/ax;
            y=z*z;
            xx=ax - 2.356194491;
            ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))));
            ans2=0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
            bj=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);

            if(besarg < 0.0){bj*=-1;}
        }

        float d1 = qq*rr*sin(alpha);

        if (besarg == 0.0) {be = sin(alpha);}
        else {be = bj*bj*4.0*sin(alpha)/(d1*d1);}
        if(siarg == 0.0) {si = 1.0;}
        else{si = sin(siarg)*sin(siarg)/(siarg*siarg);}

        form = be*si/sin(alpha);
        answer = sub*sub*form*acos(-1.0)*rr*rr*h*1.0e8*scale;

        _ptvalue[i] = radius_weight*length_weight*theta_weight*phi_weight*answer*pow(rr,2)*h;
        if (size>1) {
            _ptvalue[i] *= fabs(cos(cyl_theta*pi/180.0));
        }
}
}





























