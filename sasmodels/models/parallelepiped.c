double form_volume(double length_a, double length_b, double length_c);
double Iq(double q, double sld, double solvent_sld, double length_a, double length_b, double length_c);
double Iqxy(double qx, double qy, double sld, double solvent_sld,
    double length_a, double length_b, double length_c, double theta, double phi, double psi);

// From Igor library
double _pkernel(double a, double b,double c, double ala, double alb, double alc);
double _pkernel(double a, double b,double c, double ala, double alb, double alc){
    double argA,argB,argC,tmp1,tmp2,tmp3;
    //handle arg=0 separately, as sin(t)/t -> 1 as t->0
    argA = 0.5*a*ala;
    argB = 0.5*b*alb;
    argC = 0.5*c*alc;
    if(argA==0.0) {
        tmp1 = 1.0;
    } else {
        tmp1 = sin(argA)*sin(argA)/argA/argA;
    }
    if (argB==0.0) {
        tmp2 = 1.0;
    } else {
        tmp2 = sin(argB)*sin(argB)/argB/argB;
    }
    if (argC==0.0) {
        tmp3 = 1.0;
    } else {
        tmp3 = sin(argC)*sin(argC)/argC/argC;
    }
    return (tmp1*tmp2*tmp3);

}


double form_volume(double length_a, double length_b, double length_c)
{
    return length_a * length_b * length_c;
}


double Iq(double q,
    double sld,
    double solvent_sld,
    double length_a,
    double length_b,
    double length_c)
{
    double tmp1, tmp2;
    
    double mu = q * length_b;
    
    // Scale sides by B
    double a_scaled = length_a / length_b;
    double c_scaled = length_c / length_b;
        
    //Order of integration
    int nordi=76;			        
    int nordj=76;

    // outer integral (with nordi gauss points), integration limits = 0, 1
    double summ = 0; //initialize integral

    for( int i=0; i<nordi; i++) {
		
        // inner integral (with nordj gauss points), integration limits = 0, 1
	
        double summj = 0.0;
	    double sigma = 0.5 * ( Gauss76Z[i] + 1.0 );		
		
	    for(int j=0; j<nordj; j++) {

            double uu = 0.5 * ( Gauss76Z[j] + 1.0 );
            double mudum = mu * sqrt(1.0-sigma*sigma);
	        double arg1 = 0.5 * mudum * cos(0.5*M_PI*uu);
	        double arg2 = 0.5 * mudum * a_scaled * sin(0.5*M_PI*uu);
            if(arg1==0.0) {
	        tmp1 = 1.0;
            } else {
                tmp1 = sin(arg1)*sin(arg1)/arg1/arg1;
            }
            if (arg2==0.0) {
                tmp2 = 1.0;
            } else {
                tmp2 = sin(arg2)*sin(arg2)/arg2/arg2;
            }

            summj += Gauss76Wt[j] * tmp1 * tmp2;
        }
		
        // value of the inner integral
        double answer = 0.5 * summj;

        double arg = 0.5 * mu * c_scaled * sigma;
        if ( arg == 0.0 ) {
            answer *= 1.0;
        } else {
            answer *= sin(arg)*sin(arg)/arg/arg;
        }
		
	    // sum of outer integral
        summ += Gauss76Wt[i] * answer;
        
    }	
   
    const double vd = (sld-solvent_sld) * form_volume(length_a, length_b, length_c);
    
    // convert from [1e-12 A-1] to [cm-1] and 0.5 factor for outer integral
    return 1.0e-4 * 0.5 * vd * vd * summ;
    
}


double Iqxy(double qx, double qy,
    double sld,
    double solvent_sld,
    double length_a,
    double length_b,
    double length_c,
    double theta,
    double phi,
    double psi)
{
    double q = sqrt(qx*qx+qy*qy);
    double qx_scaled = qx/q;
    double qy_scaled = qy/q;

    // Convert angles given in degrees to radians
    theta *= M_PI_180;
    phi   *= M_PI_180;
    psi   *= M_PI_180;
    
    // Parallelepiped c axis orientation
    double cparallel_x = cos(theta) * cos(phi);
    double cparallel_y = sin(theta);
    
    // Compute angle between q and parallelepiped axis
    double cos_val_c = cparallel_x*qx_scaled + cparallel_y*qy_scaled;// + cparallel_z*qz;

    // Parallelepiped a axis orientation
    double parallel_x = -cos(phi)*sin(psi) * sin(theta)+sin(phi)*cos(psi);
    double parallel_y = sin(psi)*cos(theta);
    double cos_val_a = parallel_x*qx_scaled + parallel_y*qy_scaled;

    // Parallelepiped b axis orientation
    double bparallel_x = -sin(theta)*cos(psi)*cos(phi)-sin(psi)*sin(phi);
    double bparallel_y = cos(theta)*cos(psi);
    double cos_val_b = bparallel_x*qx_scaled + bparallel_y*qy_scaled;

    // The following tests should always pass
    if (fabs(cos_val_c)>1.0) {
      //printf("parallel_ana_2D: Unexpected error: cos(alpha)>1\n");
      cos_val_c = 1.0;
    }
    if (fabs(cos_val_a)>1.0) {
      //printf("parallel_ana_2D: Unexpected error: cos(alpha)>1\n");
      cos_val_a = 1.0;
    }
    if (fabs(cos_val_b)>1.0) {
      //printf("parallel_ana_2D: Unexpected error: cos(alpha)>1\n");
      cos_val_b = 1.0;
    }
    
    // Call the IGOR library function to get the kernel
    double form = _pkernel( q*length_a, q*length_b, q*length_c, cos_val_a, cos_val_b, cos_val_c);
  
    // Multiply by contrast^2
    const double vd = (sld - solvent_sld) * form_volume(length_a, length_b, length_c);
    return 1.0e-4 * vd * vd * form;
}
