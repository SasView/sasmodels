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
bicelle_kernel(double q,
              double rad,
              double radthick,
              double facthick,
              double halflength,
              double rhoc,
              double rhoh,
              double rhor,
              double rhosolv,
              double sin_alpha,
              double cos_alpha)
{
    const double dr1 = rhoc-rhoh;
    const double dr2 = rhor-rhosolv;
    const double dr3 = rhoh-rhor;
    const double vol1 = M_PI*square(rad)*2.0*(halflength);
    const double vol2 = M_PI*square(rad+radthick)*2.0*(halflength+facthick);
    const double vol3 = M_PI*square(rad)*2.0*(halflength+facthick);

    const double be1 = sas_2J1x_x(q*(rad)*sin_alpha);
    const double be2 = sas_2J1x_x(q*(rad+radthick)*sin_alpha);
    const double si1 = sas_sinx_x(q*(halflength)*cos_alpha);
    const double si2 = sas_sinx_x(q*(halflength+facthick)*cos_alpha);

    const double t = vol1*dr1*si1*be1 +
                     vol2*dr2*si2*be2 +
                     vol3*dr3*si2*be1;

    const double retval = t*t;

    return retval;

}

static double
bicelle_integration(double q,
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
    const double halflength = 0.5*length;

    double summ = 0.0;
    for(int i=0;i<N_POINTS_76;i++) {
        double alpha = (Gauss76Z[i] + 1.0)*uplim;
        double sin_alpha, cos_alpha; // slots to hold sincos function output
        SINCOS(alpha, sin_alpha, cos_alpha);
        double yyy = Gauss76Wt[i] * bicelle_kernel(q, rad, radthick, facthick,
                             halflength, rhoc, rhoh, rhor, rhosolv,
                             sin_alpha, cos_alpha);
        summ += yyy*sin_alpha;
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
                           solvent_sld, sin_alpha, cos_alpha);
    return 1.0e-4*answer;
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
