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
    double si1,si2,be1,be2;

    const double dr1 = rhoc-rhoh;
    const double dr2 = rhor-rhosolv;
    const double dr3 = rhoh-rhor;
    const double vol1 = M_PI*rad*rad*(2.0*length);
    const double vol2 = M_PI*(rad+radthick)*(rad+radthick)*2.0*(length+facthick);
    const double vol3 = M_PI*rad*rad*2.0*(length+facthick);
    double sn,cn;
    SINCOS(dum, sn, cn);
    double besarg1 = qq*rad*sn;
    double besarg2 = qq*(rad+radthick)*sn;
    double sinarg1 = qq*length*cn;
    double sinarg2 = qq*(length+facthick)*cn;

    be1 = sas_J1c(besarg1);
    be2 = sas_J1c(besarg2);
    si1 = sinc(sinarg1);
    si2 = sinc(sinarg2);

    const double t = vol1*dr1*si1*be1 +
                     vol2*dr2*si2*be2 +
                     vol3*dr3*si2*be1;

    const double retval = t*t*sn;

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
    // set up the integration end points
    const double uplim = M_PI/4;
    const double halfheight = length/2.0;

    double summ = 0.0;
    for(int i=0;i<N_POINTS_76;i++) {
        double zi = (Gauss76Z[i] + 1.0)*uplim;
        double yyy = Gauss76Wt[i] * bicelle_kernel(qq, rad, radthick, facthick,
                             halfheight, rhoc, rhoh, rhor,rhosolv, zi);
        summ += yyy;
    }

    // calculate value of integral to return
    double answer = uplim*summ;
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
    //convert angle degree to radian
    theta *= M_PI_180;
    phi *= M_PI_180;

    // Cylinder orientation
    const double cyl_x = cos(theta) * cos(phi);
    const double cyl_y = sin(theta);

    // Compute the angle btw vector q and the axis of the cylinder
    const double cos_val = cyl_x*q_x + cyl_y*q_y;
    const double alpha = acos( cos_val );

    // Get the kernel
    double answer = bicelle_kernel(q, radius, rim_thickness, face_thickness,
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
