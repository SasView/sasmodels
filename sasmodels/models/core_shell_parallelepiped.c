double form_volume(double a_side, double b_side, double c_side, 
                   double arim_thickness, double brim_thickness, double crim_thickness);
double Iq(double q, double core_sld, double arim_sld, double brim_sld, double crim_sld,
          double solvent_sld, double a_side, double b_side, double c_side,
          double arim_thickness, double brim_thickness, double crim_thickness);
double Iqxy(double qx, double qy, double core_sld, double arim_sld, double brim_sld,
            double crim_sld, double solvent_sld, double a_side, double b_side,
            double c_side, double arim_thickness, double brim_thickness,
            double crim_thickness, double theta, double phi, double psi);

double form_volume(double a_side, double b_side, double c_side, 
                   double arim_thickness, double brim_thickness, double crim_thickness)
{
    //return a_side * b_side * c_side;
    return a_side * b_side * c_side + 
           2.0 * arim_thickness * b_side * c_side + 
           2.0 * brim_thickness * a_side * c_side +
           2.0 * crim_thickness * a_side * b_side;
}

double Iq(double q,
    double core_sld,
    double arim_sld,
    double brim_sld,
    double crim_sld,
    double solvent_sld,
    double a_side,
    double b_side,
    double c_side,
    double arim_thickness,
    double brim_thickness,
    double crim_thickness)
{
    // Code converted from functions CSPPKernel and CSParallelepiped in libCylinder.c_scaled
    // Did not understand the code completely, it should be rechecked (Miguel Gonzalez)
    
    double t1, t2, t3, t4, tmp, answer;   
    double mu = q * b_side;
    
    //calculate volume before rescaling (in original code, but not used)
    //double vol = form_volume(a_side, b_side, c_side, arim_thickness, brim_thickness, crim_thickness);		
    //double vol = a_side * b_side * c_side + 
    //       2.0 * arim_thickness * b_side * c_side + 
    //       2.0 * brim_thickness * a_side * c_side +
    //       2.0 * crim_thickness * a_side * b_side;
    
    // Scale sides by B
    double a_scaled = a_side / b_side;
    double c_scaled = c_side / b_side;
    double arim_scaled = arim_thickness / b_side;
    double brim_scaled = brim_thickness / b_side;
       
    // DelRho values (note that drC is not used later)       
	double dr0 = core_sld-solvent_sld;
	double drA = arim_sld-solvent_sld;
	double drB = brim_sld-solvent_sld;
	//double drC = crim_sld-solvent_sld;

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

	        double Vin = a_side * b_side * c_side;
	        //double Vot = (a_side * b_side * c_side +
            //            2.0 * arim_thickness * b_side * c_side +
            //            2.0 * a_side * brim_thickness * c_side +
            //            2.0 * a_side * b_side * crim_thickness);
	        double V1 = (2.0 * arim_thickness * b_side * c_side);    // incorrect V1 (aa*bb*cc+2*ta*bb*cc)
	        double V2 = (2.0 * a_side * brim_thickness * c_side);    // incorrect V2(aa*bb*cc+2*aa*tb*cc)
	
            // ta and tb correspond to the definitions in CSPPKernel, but they don't make sense to me (MG)
            // the a_scaled in the definition of tb was present in CSPPKernel in libCylinder.c,
            // while in cspkernel in csparallelepiped.cpp (used for the 2D), all the definitions
            // for ta, tb, tc use also A + 2*rim_thickness (but not scaled by B!!!)            
            double ta = (a_scaled+2.0*arim_thickness)/b_side; 
            double tb = (a_scaled+2.0*brim_thickness)/b_side;
    
	        double arg1 = (0.5*mudum*a_scaled) * sin(0.5*M_PI*uu);
	        double arg2 = (0.5*mudum) * cos(0.5*M_PI*uu);
	        double arg3=  (0.5*mudum*ta) * sin(0.5*M_PI*uu);
	        double arg4=  (0.5*mudum*tb) * cos(0.5*M_PI*uu);

	        if(arg1==0.0){
		        t1 = 1.0;
	        } else {
		        t1 = (sin(arg1)/arg1);                //defn for CSPP model sin(arg1)/arg1    test:  (sin(arg1)/arg1)*(sin(arg1)/arg1)   
	        }
	        if(arg2==0.0){
		        t2 = 1.0;
	        } else {
		        t2 = (sin(arg2)/arg2);           //defn for CSPP model sin(arg2)/arg2   test: (sin(arg2)/arg2)*(sin(arg2)/arg2)    
	        }	
	        if(arg3==0.0){
		        t3 = 1.0;
	        } else {
		        t3 = sin(arg3)/arg3;
	        }
	        if(arg4==0.0){
		        t4 = 1.0;
	        } else {
		        t4 = sin(arg4)/arg4;
	        }
            
            // Expression in libCylinder.c (neither drC nor Vot are used)
	        tmp =( dr0*t1*t2*Vin + drA*(t3-t1)*t2*V1+ drB*t1*(t4-t2)*V2 )*( dr0*t1*t2*Vin + drA*(t3-t1)*t2*V1+ drB*t1*(t4-t2)*V2 );   //  correct FF : square of sum of phase factors            
            
            // To note also that in csparallelepiped.cpp, there is a function called
            // cspkernel, which seems to make more sense and has the following comment:
            //   This expression is different from NIST/IGOR package (I strongly believe the IGOR is wrong!!!). 10/15/2010.
            //   tmp =( dr0*tmp1*tmp2*tmp3*Vin + drA*(tmpt1-tmp1)*tmp2*tmp3*V1+ drB*tmp1*(tmpt2-tmp2)*tmp3*V2 + drC*tmp1*tmp2*(tmpt3-tmp3)*V3)*
            //   ( dr0*tmp1*tmp2*tmp3*Vin + drA*(tmpt1-tmp1)*tmp2*tmp3*V1+ drB*tmp1*(tmpt2-tmp2)*tmp3*V2 + drC*tmp1*tmp2*(tmpt3-tmp3)*V3);   //  correct FF : square of sum of phase factors
            // This is the function called by csparallelepiped_analytical_2D_scaled,
            // while CSParallelepipedModel calls CSParallelepiped in libCylinder.c        
            
            summj += Gauss76Wt[j] * tmp;

        }
		
        // value of the inner integral
        answer = 0.5 * summj;

		// finish the outer integral 
		double arg = 0.5 * mu* c_scaled *sigma;
		if ( arg == 0.0 ) {
			answer *= 1.0;
		} else {
	        answer *= sin(arg)*sin(arg)/arg/arg;
		}
		
		// now sum up the outer integral
		summ += Gauss76Wt[i] * answer;

    }
    
	answer = 0.5 * summ;

	//convert from [1e-12 A-1] to [cm-1]
	answer *= 1.0e-4;
	
	return answer;
}

double Iqxy(double qx, double qy,
    double core_sld,
    double arim_sld,
    double brim_sld,
    double crim_sld,
    double solvent_sld,
    double a_side,
    double b_side,
    double c_side,
    double arim_thickness,
    double brim_thickness,
    double crim_thickness,
    double theta,
    double phi,
    double psi)
{
    double tmp1, tmp2, tmp3, tmpt1, tmpt2, tmpt3;   

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
    
    // cspkernel in csparallelepiped recoded here 
    double dr0 = core_sld-solvent_sld;
	double drA = arim_sld-solvent_sld;
	double drB = brim_sld-solvent_sld;
	double drC = crim_sld-solvent_sld;
	double Vin = a_side * b_side * c_side;
    // As for the 1D case, Vot is not used
	//double Vot = (a_side * b_side * c_side +
    //              2.0 * arim_thickness * b_side * c_side +
    //              2.0 * a_side * brim_thickness * c_side +
    //              2.0 * a_side * b_side * crim_thickness);
	double V1 = (2.0 * arim_thickness * b_side * c_side);    // incorrect V1 (aa*bb*cc+2*ta*bb*cc)
	double V2 = (2.0 * a_side * brim_thickness * c_side);    // incorrect V2(aa*bb*cc+2*aa*tb*cc)
    double V3 = (2.0 * a_side * b_side * crim_thickness);
    // The definitions of ta, tb, tc are not the same as in the 1D case because there is no
    // the scaling by B. The use of a_side for the 3 of them seems clearly a mistake to me,
    // but for the moment I let it like this until understanding better the code.
	double ta = a_side + 2.0*arim_thickness;
    double tb = a_side + 2.0*brim_thickness;
    double tc = a_side + 2.0*crim_thickness;
    //handle arg=0 separately, as sin(t)/t -> 1 as t->0
    double argA = 0.5*q*a_side*cos_val_a;
    double argB = 0.5*q*b_side*cos_val_b;
    double argC = 0.5*q*c_side*cos_val_c;
    double argtA = 0.5*q*ta*cos_val_a;
    double argtB = 0.5*q*tb*cos_val_b;
    double argtC = 0.5*q*tc*cos_val_c;
    
    if(argA==0.0) {
        tmp1 = 1.0;
    } else {
        tmp1 = sin(argA)/argA;
    }
    if (argB==0.0) {
        tmp2 = 1.0;
    } else {
        tmp2 = sin(argB)/argB;
    }
    if (argC==0.0) {
        tmp3 = 1.0;
    } else {
        tmp3 = sin(argC)/argC;
    }
    if(argtA==0.0) {
        tmpt1 = 1.0;
    } else {
        tmpt1 = sin(argtA)/argtA;
    }
    if (argtB==0.0) {
        tmpt2 = 1.0;
    } else {
        tmpt2 = sin(argtB)/argtB;
    }
    if (argtC==0.0) {
        tmpt3 = 1.0;
    } else {
        tmpt3 = sin(argtC)*sin(argtC)/argtC/argtC;
    }

    // f uses Vin, V1, V2, and V3 and it seems to have more sense than the value computed
    // in the 1D code, but should be checked!
    double f = ( dr0*tmp1*tmp2*tmp3*Vin + drA*(tmpt1-tmp1)*tmp2*tmp3*V1 + 
               drB*tmp1*(tmpt2-tmp2)*tmp3*V2 + drC*tmp1*tmp2*(tmpt3-tmp3)*V3);
   
    return 1.0e-4 * f * f;
}
