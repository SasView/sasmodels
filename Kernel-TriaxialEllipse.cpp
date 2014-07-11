__kernel void TriaxialEllipseKernel(__global const real *qx, __global const real *qy, __global real *_ptvalue,
const real sub, const real scale, const real axisA, const real axisB, const real axisC, const real axis_phi,
const real axis_theta, const real axis_psi, const real axisA_weight, const real axisB_weight, const real axisC_weight,
const real psi_weight, const real phi_weight, const real theta_weight, const int total, const int size)
{
    int i = get_global_id(0);
    if(i < total)
    {
        real q = sqrt(qx[i]*qx[i]+qy[i]*qy[i]);
        real q_x = qx[i]/q;
        real q_y = qy[i]/q;
        real pi = 4.0*atan(1.0);
        
        //convert angle degree to radian
        real theta = axis_theta*pi/180.0;
        real phi = axis_phi*pi/180.0;
        real psi = axis_psi*pi/180.0;
        real cos_val = cos(theta)*cos(phi)*q_x + sin(theta)*q_y;
        
        // calculate the axis of the ellipse wrt q-coord.
        real cos_nu = (-cos(phi)*sin(psi)*sin(theta)+sin(phi)*cos(psi))*q_x + sin(psi)*cos(theta)*q_y;
        real cos_mu = (-sin(theta)*cos(psi)*cos(phi)-sin(psi)*sin(phi))*q_x + cos(theta)*cos(psi)*q_y;

        real answer=0;
        real t = q*sqrt(axisA*axisA*cos_nu*cos_nu+axisB*axisB*cos_mu*cos_mu+axisC*axisC*cos_val*cos_val);

        if (t==0.0){
            answer = 1.0;
        }
        else{
            answer  = 3.0*(sin(t)-t*cos(t))/(t*t*t);
        }
        answer*=answer*sub*sub*4.0*pi/3.0*axisA*axisB*axisC*1.0e8*scale;

        _ptvalue[i] = axisA_weight*axisB_weight*axisC_weight*theta_weight*phi_weight*psi_weight*answer*axisA*axisB*axisC;
        if (size>1)
        {
            _ptvalue[i] *= fabs(cos(theta*pi/180.0));
        }
    }
}




