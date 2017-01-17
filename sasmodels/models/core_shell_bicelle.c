double form_volume(double radius, double thick_rim, double thick_face, double length);
double Iq(double q,
          double radius,
          double thick_rim,
          double thick_face,
          double length,
          double core_sld,
          double face_sld,
          double rim_sld,
          double solvent_sld);


double Iqxy(double qx, double qy,
          double radius,
          double thick_rim,
          double thick_face,
          double length,
          double core_sld,
          double face_sld,
          double rim_sld,
          double solvent_sld,
          double theta,
          double phi);


double form_volume(double radius, double thick_rim, double thick_face, double length)
{
    return M_PI*(radius+thick_rim)*(radius+thick_rim)*(length+2.0*thick_face);
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
              double sin_alpha,
              double cos_alpha)
{
    double si1,si2,be1,be2;

    const double dr1 = rhoc-rhoh;
    const double dr2 = rhor-rhosolv;
    const double dr3 = rhoh-rhor;
    const double vol1 = M_PI*rad*rad*(2.0*length);
    const double vol2 = M_PI*(rad+radthick)*(rad+radthick)*2.0*(length+facthick);
    const double vol3 = M_PI*rad*rad*2.0*(length+facthick);
    double besarg1 = qq*rad*sin_alpha;
    double besarg2 = qq*(rad+radthick)*sin_alpha;
    double sinarg1 = qq*length*cos_alpha;
    double sinarg2 = qq*(length+facthick)*cos_alpha;

    be1 = sas_2J1x_x(besarg1);
    be2 = sas_2J1x_x(besarg2);
    si1 = sas_sinx_x(sinarg1);
    si2 = sas_sinx_x(sinarg2);

    const double t = vol1*dr1*si1*be1 +
                     vol2*dr2*si2*be2 +
                     vol3*dr3*si2*be1;

    const double retval = t*t*sin_alpha;

    return retval;

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
    const double uplim = M_PI_4;
    const double halfheight = 0.5*length;

    double summ = 0.0;
    for(int i=0;i<N_POINTS_76;i++) {
        double alpha = (Gauss76Z[i] + 1.0)*uplim;
        double sin_alpha, cos_alpha; // slots to hold sincos function output
        SINCOS(alpha, sin_alpha, cos_alpha);
        double yyy = Gauss76Wt[i] * bicelle_kernel(qq, rad, radthick, facthick,
                             halfheight, rhoc, rhoh, rhor, rhosolv,
                             sin_alpha, cos_alpha);
        summ += yyy;
    }

    // calculate value of integral to return
    double answer = uplim*summ;
    return answer;
}

static double
bicelle_kernel_2d(double qx, double qy,
          double radius,
          double thick_rim,
          double thick_face,
          double length,
          double core_sld,
          double face_sld,
          double rim_sld,
          double solvent_sld,
          double theta,
          double phi)
{
    double q, sin_alpha, cos_alpha;
    ORIENT_SYMMETRIC(qx, qy, theta, phi, q, sin_alpha, cos_alpha);

    double answer = bicelle_kernel(q, radius, thick_rim, thick_face,
                           0.5*length, core_sld, face_sld, rim_sld,
                           solvent_sld, sin_alpha, cos_alpha) / fabs(sin_alpha);

    answer *= 1.0e-4;

    return answer;
}

double Iq(double q,
          double radius,
          double thick_rim,
          double thick_face,
          double length,
          double core_sld,
          double face_sld,
          double rim_sld,
          double solvent_sld)
{
    double intensity = bicelle_integration(q, radius, thick_rim, thick_face,
                       length, core_sld, face_sld, rim_sld, solvent_sld);
    return intensity*1.0e-4;
}


double Iqxy(double qx, double qy,
          double radius,
          double thick_rim,
          double thick_face,
          double length,
          double core_sld,
          double face_sld,
          double rim_sld,
          double solvent_sld,
          double theta,
          double phi)
{
    double intensity = bicelle_kernel_2d(qx, qy,
                      radius,
                      thick_rim,
                      thick_face,
                      length,
                      core_sld,
                      face_sld,
                      rim_sld,
                      solvent_sld,
                      theta,
                      phi);

    return intensity;
}
