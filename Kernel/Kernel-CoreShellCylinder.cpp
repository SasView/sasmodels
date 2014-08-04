__kernel void CoreShellCylinderKernel(__global const real *qx, __global const real *qy, __global real *_ptvalue,
const real axis_theta, const real axis_phi, const real thickness, const real length, const real radius,
const real scale, const real radius_weight, const real length_weight, const real thickness_weight,
const real theta_weight, const real phi_weight, const real core_sld, const real shell_sld, const real solvent_sld,
const int size, const int total)
{
    int i = get_global_id(0);
    if(i < total)
    {
        real q = sqrt(qx[i]*qx[i]+qy[i]*qy[i]);
        real pi = 4.0*atan(1.0);
        real theta = axis_theta*pi/180.0;
        real alpha = acos(cos(theta)*cos(axis_phi*pi/180.0)*qx[i]/q + sin(theta)*qy[i]/q);

        if (alpha == 0.0){
        alpha = 1.0e-26;
        }

        real si1=0; real si2=0; real be1=0; real be2=0;

        real vol2 = pi*(radius+thickness)*(radius+thickness)*(length+2.0*thickness);

        real besarg1 = q*radius*sin(alpha);
        real besarg2 = q*(radius+thickness)*sin(alpha);
        real sinarg1 = q*length/2.0*cos(alpha);
        real sinarg2 = q*(length/2.0+thickness)*cos(alpha);


        if (besarg1 == 0.0){be1 = 0.5;}
        else{be1 = NR_BessJ1(besarg1)/besarg1;}

        if (besarg2 == 0.0){be2 = 0.5;}
        else{be2 = NR_BessJ1(besarg2)/besarg2;}

        if (sinarg1 == 0.0){si1 = 1.0;}
        else{si1 = sin(sinarg1)/sinarg1;}

        if (sinarg2 == 0.0){si2 = 1.0;}
        else{si2 = sin(sinarg2)/sinarg2;}

        real tt = 2.0*vol2*(shell_sld-solvent_sld)*si2*be2+2.0*(pi*radius*radius*(length))*(core_sld-shell_sld)*si1*be1;

        real answer = (tt*tt)*sin(alpha)/fabs(sin(alpha));
        answer *= answer/(pi*(radius+thickness)*(radius+thickness)*(length+2.0*thickness))*1.0e8*scale;

        _ptvalue[i] += radius_weight*length_weight*thickness_weight*theta_weight*phi_weight*answer*pow(radius+thickness,2)*(length+2.0*thickness);
     //   if (size>1) {
       // _ptvalue[i] *= fabs(cos(axis_theta*pi/180.0));
        //}

    }
}



