__kernel void CylinderKernel(__global const real *qx, global const real *qy, __global real *_ptvalue, const real sub,
const real rr, const real h, const real scale, const real radius_weight, const real length_weight,
const real theta_weight, const real phi_weight, const real cyl_theta,
const real cyl_phi, const int count, const int size)
{
	// qq is the q-value for the calculation (1/A)
	// rr is the radius of the cylinder (A)
	// h is the -LENGTH of the cylinder = L (A)
    int i = get_global_id(0);

    if(i < count)
    {
        real qq = sqrt(qx[i]*qx[i]+qy[i]*qy[i]);

        real pi = 4.0*atan(1.0);
        real theta = cyl_theta*pi/180.0;
        real phi = cyl_phi*pi/180.0;

        real cyl_x = cos(theta)*cos(phi);
        real cyl_y = sin(theta);
        real cos_val = cyl_x*(qx[i]/qq) + cyl_y*(qy[i]/qq);

        real alpha = acos(cos_val);
        if(alpha == 0.0){
            alpha = 1.0e-26;
        }
        real besarg = qq*rr*sin(alpha);
        real siarg = qq*h/2*cos(alpha);
        real be=0.0; real si=0.0;

        real bj = NR_BessJ1(besarg);

        real d1 = qq*rr*sin(alpha);

        if (besarg == 0.0){
            be = sin(alpha);
        }
        else{
            be = bj*bj*4.0*sin(alpha)/(d1*d1);
        }
        if(siarg == 0.0){
            si = 1.0;
        }
        else{
            si = sin(siarg)*sin(siarg)/(siarg*siarg);
        }

        real form = be*si/sin(alpha);
        real answer = sub*sub*form*acos(-1.0)*rr*rr*h*1.0e8*scale;

        _ptvalue[i] = radius_weight*length_weight*theta_weight*phi_weight*answer*pow(rr,2)*h;
        if (size>1) {
            _ptvalue[i] *= fabs(cos(cyl_theta*pi/180.0));
        }
    }
}





























