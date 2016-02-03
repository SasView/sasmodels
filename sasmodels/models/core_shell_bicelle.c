double form_volume(double radius, double rim_thickness, double face_thickness, double length);
double Iq(double q,
          double radius,
          double rim_thickness,
          double face_thickness,
          double length,
          double core_sld,
          double face_sld,
          double rim_sld,
          double solvent_sld);


double Iqxy(double qx, double qy,
          double radius,
          double rim_thickness,
          double face_thickness,
          double length,
          double core_sld,
          double face_sld,
          double rim_sld,
          double solvent_sld,
          double theta,
          double phi);


double form_volume(double radius, double rim_thickness, double face_thickness, double length)
{
    return M_PI*(radius+rim_thickness)*(radius+rim_thickness)*(length+2*face_thickness);
}

static double
bicelle_kernel(double qq,
              double rad,
              double radthick,
              double facthick,
              double length,
              double rhoc,
              double rhoh,
              double rhor,
              double rhosolv,
              double dum)
{
	double dr1,dr2,dr3;
	double besarg1,besarg2;
	double vol1,vol2,vol3;
	double sinarg1,sinarg2;
	double t1,t2,t3;
	double retval,si1,si2,be1,be2;

	dr1 = rhoc-rhoh;
	dr2 = rhor-rhosolv;
	dr3=  rhoh-rhor;
	vol1 = M_PI*rad*rad*(2.0*length);
	vol2 = M_PI*(rad+radthick)*(rad+radthick)*(2.0*length+2.0*facthick);
	vol3= M_PI*(rad)*(rad)*(2.0*length+2.0*facthick);
	besarg1 = qq*rad*sin(dum);
	besarg2 = qq*(rad+radthick)*sin(dum);
	sinarg1 = qq*length*cos(dum);
	sinarg2 = qq*(length+facthick)*cos(dum);

	if(besarg1 == 0) {
		be1 = 0.5;
	} else {
		be1 = J1(besarg1)/besarg1;
	}
	if(besarg2 == 0) {
		be2 = 0.5;
	} else {
		be2 = J1(besarg2)/besarg2;
	}
	if(sinarg1 == 0) {
		si1 = 1.0;
	} else {
		si1 = sin(sinarg1)/sinarg1;
	}
	if(sinarg2 == 0) {
		si2 = 1.0;
	} else {
		si2 = sin(sinarg2)/sinarg2;
	}
	t1 = 2.0*vol1*dr1*si1*be1;
	t2 = 2.0*vol2*dr2*si2*be2;
	t3 = 2.0*vol3*dr3*si2*be1;

	retval = ((t1+t2+t3)*(t1+t2+t3))*sin(dum);
	return(retval);

}

static double
bicelle_integration(double qq,
                   double rad,
                   double radthick,
                   double facthick,
                   double length,
                   double rhoc,
                   double rhoh,
                   double rhor,
                   double rhosolv)
{


	double answer,halfheight;
	double lolim,uplim,summ,yyy,zi;
	int nord,i;

	// set up the integration end points
	nord = 76;
	lolim = 0.0;
	uplim = M_PI/2;
	halfheight = length/2.0;

	summ = 0.0;
	i=0;
	for(i=0;i<nord;i++) {
		zi = ( Gauss76Z[i]*(uplim-lolim) + uplim + lolim )/2.0;
		yyy = Gauss76Wt[i] * bicelle_kernel(qq, rad, radthick, facthick,
		                     halfheight, rhoc, rhoh, rhor,rhosolv, zi);
		summ += yyy;
	}

	// calculate value of integral to return
	answer = (uplim-lolim)/2.0*summ;
	return(answer);
}

static double
bicelle_kernel_2d(double q, double q_x, double q_y,
          double radius,
          double rim_thickness,
          double face_thickness,
          double length,
          double core_sld,
          double face_sld,
          double rim_sld,
          double solvent_sld,
          double theta,
          double phi)
{
    double cyl_x, cyl_y;
    double alpha, cos_val;
    double answer;

    //convert angle degree to radian
    theta *= M_PI/180.0;
    phi *= M_PI/180.0;

    // Cylinder orientation
    cyl_x = cos(theta) * cos(phi);
    cyl_y = sin(theta);

    // Compute the angle btw vector q and the axis of the cylinder
    cos_val = cyl_x*q_x + cyl_y*q_y;
    alpha = acos( cos_val );

    // Get the kernel
    answer = bicelle_kernel(q, radius, rim_thickness, face_thickness,
                           length/2.0, core_sld, face_sld, rim_sld,
                           solvent_sld, alpha) / fabs(sin(alpha));

    answer *= 1.0e-4;

    return answer;
}

double Iq(double q,
          double radius,
          double rim_thickness,
          double face_thickness,
          double length,
          double core_sld,
          double face_sld,
          double rim_sld,
          double solvent_sld)
{
    double intensity = bicelle_integration(q, radius, rim_thickness, face_thickness,
                       length, core_sld, face_sld, rim_sld, solvent_sld);
    return intensity*1.0e-4;
}


double Iqxy(double qx, double qy,
          double radius,
          double rim_thickness,
          double face_thickness,
          double length,
          double core_sld,
          double face_sld,
          double rim_sld,
          double solvent_sld,
          double theta,
          double phi)
{
    double q;
    q = sqrt(qx*qx+qy*qy);
    double intensity = bicelle_kernel_2d(q, qx/q, qy/q,
                      radius,
                      rim_thickness,
                      face_thickness,
                      length,
                      core_sld,
                      face_sld,
                      rim_sld,
                      solvent_sld,
                      theta,
                      phi);

    return intensity;
}
