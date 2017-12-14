static double
form_volume(double radius, double thick_rim, double thick_face, double length)
{
    return M_PI*square(radius+thick_rim)*(length+2.0*thick_face);
}

static double
bicelle_kernel(double qab,
    double qc,
    double radius,
    double thick_radius,
    double thick_face,
    double halflength,
    double sld_core,
    double sld_face,
    double sld_rim,
    double sld_solvent)
{
    const double dr1 = sld_core-sld_face;
    const double dr2 = sld_rim-sld_solvent;
    const double dr3 = sld_face-sld_rim;
    const double vol1 = M_PI*square(radius)*2.0*(halflength);
    const double vol2 = M_PI*square(radius+thick_radius)*2.0*(halflength+thick_face);
    const double vol3 = M_PI*square(radius)*2.0*(halflength+thick_face);

    const double be1 = sas_2J1x_x((radius)*qab);
    const double be2 = sas_2J1x_x((radius+thick_radius)*qab);
    const double si1 = sas_sinx_x((halflength)*qc);
    const double si2 = sas_sinx_x((halflength+thick_face)*qc);

    const double t = vol1*dr1*si1*be1 +
                     vol2*dr2*si2*be2 +
                     vol3*dr3*si2*be1;

    return t;
}

static double
Iq(double q,
    double radius,
    double thick_radius,
    double thick_face,
    double length,
    double sld_core,
    double sld_face,
    double sld_rim,
    double sld_solvent)
{
    // set up the integration end points
    const double uplim = M_PI_4;
    const double halflength = 0.5*length;

    double total = 0.0;
    for(int i=0;i<GAUSS_N;i++) {
        double theta = (GAUSS_Z[i] + 1.0)*uplim;
        double sin_theta, cos_theta; // slots to hold sincos function output
        SINCOS(theta, sin_theta, cos_theta);
        double fq = bicelle_kernel(q*sin_theta, q*cos_theta, radius, thick_radius, thick_face,
                                   halflength, sld_core, sld_face, sld_rim, sld_solvent);
        total += GAUSS_W[i]*fq*fq*sin_theta;
    }

    // calculate value of integral to return
    double answer = total*uplim;
    return 1.0e-4*answer;
}

static double
Iqac(double qab, double qc,
    double radius,
    double thick_rim,
    double thick_face,
    double length,
    double core_sld,
    double face_sld,
    double rim_sld,
    double solvent_sld)
{
    double fq = bicelle_kernel(qab, qc, radius, thick_rim, thick_face,
                           0.5*length, core_sld, face_sld, rim_sld,
                           solvent_sld);
    return 1.0e-4*fq*fq;
}