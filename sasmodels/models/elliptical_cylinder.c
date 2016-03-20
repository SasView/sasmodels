double form_volume(double r_minor, double r_ratio, double length);
double Iq(double q, double r_minor, double r_ratio, double length,
          double sld, double solvent_sld);
double Iqxy(double qx, double qy, double r_minor, double r_ratio, double length,
            double sld, double solvent_sld, double theta, double phi, double psi);


double _elliptical_cylinder_kernel(double q, double r_minor, double r_ratio, double theta);

double _elliptical_cylinder_kernel(double q, double r_minor, double r_ratio, double theta)
{
    // This is the function LAMBDA1^2 in Feigin's notation
    // q is the q-value for the calculation (1/A)
    // r_minor is the transformed radius"a" in Feigin's notation
    // r_ratio is the ratio (major radius)/(minor radius) of the Ellipsoid [=] ---
    // theta is the dummy variable of the integration

    double retval,arg;

    arg = q*r_minor*sqrt((1.0+r_ratio*r_ratio)/2+(1.0-r_ratio*r_ratio)*cos(theta)/2);
    if (arg == 0.0){
        retval = 1.0;
    }else{
        //retval = 2.0*NR_BessJ1(arg)/arg;
        retval = sas_J1c(arg);
    }
    return retval*retval ;
}


double form_volume(double r_minor, double r_ratio, double length)
{
    return M_PI * r_minor * r_minor * r_ratio * length;
}

double Iq(double q, double r_minor, double r_ratio, double length,
          double sld, double solvent_sld) {

    const int nordi=76; //order of integration
    const int nordj=20;
    double va,vb;       //upper and lower integration limits
    double summ,zi,yyy,answer;         //running tally of integration
    double summj,vaj,vbj,zij,arg,si;            //for the inner integration

    // orientational average limits
    va = 0.0;
    vb = 1.0;
    // inner integral limits
    vaj=0.0;
    vbj=M_PI;

    //initialize integral
    summ = 0.0;

    const double delrho = sld - solvent_sld;

    for(int i=0;i<nordi;i++) {
        //setup inner integral over the ellipsoidal cross-section
        summj=0;
        zi = ( Gauss76Z[i]*(vb-va) + va + vb )/2.0;     //the "x" dummy
        arg = r_minor*sqrt(1.0-zi*zi);
        for(int j=0;j<nordj;j++) {
            //20 gauss points for the inner integral
            zij = ( Gauss20Z[j]*(vbj-vaj) + vaj + vbj )/2.0;        //the "y" dummy
            yyy = Gauss20Wt[j] * _elliptical_cylinder_kernel(q, arg, r_ratio, zij);
            summj += yyy;
        }
        //now calculate the value of the inner integral
        answer = (vbj-vaj)/2.0*summj;
        //divide integral by Pi
        answer /= M_PI;

        //now calculate outer integral
        arg = q*length*zi/2.0;
        if (arg == 0.0){
            si = 1.0;
        }else{
            si = sin(arg) * sin(arg) / arg / arg;
        }
        yyy = Gauss76Wt[i] * answer * si;
        summ += yyy;
    }

    answer = (vb-va)/2.0*summ;
    // Multiply by contrast^2
    answer *= delrho*delrho;

    const double vol = form_volume(r_minor, r_ratio, length);
    return answer*vol*vol*1e-4;
}


double Iqxy(double qx, double qy, double r_minor, double r_ratio, double length,
            double sld, double solvent_sld, double theta, double phi, double psi) {
    const double _theta = theta * M_PI / 180.0;
    const double _phi = phi * M_PI / 180.0;
    const double _psi = psi * M_PI / 180.0;
    const double q = sqrt(qx*qx+qy*qy);
    const double q_x = qx/q;
    const double q_y = qy/q;

    //Cylinder orientation
    double cyl_x = cos(_theta) * cos(_phi);
    double cyl_y = sin(_theta);

    //cyl_z = -cos(_theta) * sin(_phi);

    // q vector
    //q_z = 0;

    // Note: cos(alpha) = 0 and 1 will get an
    // undefined value from CylKernel
    //alpha = acos( cos_val );

    //ellipse orientation:
    // the elliptical corss section was transformed and projected
    // into the detector plane already through sin(alpha)and furthermore psi remains as same
    // on the detector plane.
    // So, all we need is to calculate the angle (nu) of the minor axis of the ellipse wrt
    // the wave vector q.

    //x- y- component on the detector plane.
    const double ella_x =  -cos(_phi)*sin(_psi) * sin(_theta)+sin(_phi)*cos(_psi);
    const double ella_y =  sin(_psi)*cos(_theta);
    const double ellb_x =  -sin(_theta)*cos(_psi)*cos(_phi)-sin(_psi)*sin(_phi);
    const double ellb_y =  cos(_theta)*cos(_psi);

    // Compute the angle btw vector q and the
    // axis of the cylinder
    double cos_val = cyl_x*q_x + cyl_y*q_y;// + cyl_z*q_z;

    // calculate the axis of the ellipse wrt q-coord.
    double cos_nu = ella_x*q_x + ella_y*q_y;
    double cos_mu = ellb_x*q_x + ellb_y*q_y;

    // The following test should always pass
    if (fabs(cos_val)>1.0) {
      //printf("cyl_ana_2D: Unexpected error: cos(alpha)>1\n");
      cos_val = 1.0;
    }
    if (fabs(cos_nu)>1.0) {
      //printf("cyl_ana_2D: Unexpected error: cos(nu)>1\n");
      cos_nu = 1.0;
    }
    if (fabs(cos_mu)>1.0) {
      //printf("cyl_ana_2D: Unexpected error: cos(nu)>1\n");
      cos_mu = 1.0;
    }

    const double r_major = r_ratio * r_minor;
    const double qr = q*sqrt( r_major*r_major*cos_nu*cos_nu + r_minor*r_minor*cos_mu*cos_mu );
    const double qL = q*length*cos_val/2.0;

    double Be;
    if (qr==0){
      Be = 0.5;
    }else{
      //Be = NR_BessJ1(qr)/qr;
      Be = 0.5*sas_J1c(qr);
    }

    double Si;
    if (qL==0){
      Si = 1.0;
    }else{
      Si = sin(qL)/qL;
    }

    const double k = 2.0 * Be * Si;
    const double vol = form_volume(r_minor, r_ratio, length);
    return (sld - solvent_sld) * (sld - solvent_sld) * k * k *vol*vol*1.0e-4;
}
