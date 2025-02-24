

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
radius_from_excluded_volume(double radius, double thick_rim, double thick_face, double length)
{
    const double radius_tot = radius + thick_rim;
    const double length_tot = length + 2.0*thick_face;
    return 0.5*cbrt(0.75*radius_tot*(2.0*radius_tot*length_tot + (radius_tot + length_tot)*(M_PI*radius_tot + length_tot)));
}

static double
radius_from_volume(double radius, double thick_rim, double thick_face, double length)
{
    const double volume_bicelle = form_volume(radius,thick_rim,thick_face,length);
    return cbrt(volume_bicelle/M_4PI_3);
}

static double
radius_from_diagonal(double radius, double thick_rim, double thick_face, double length)
{
    const double radius_tot = radius + thick_rim;
    const double length_tot = length + 2.0*thick_face;
    return sqrt(radius_tot*radius_tot + 0.25*length_tot*length_tot);
}

static double
radius_effective(int mode, double radius, double thick_rim, double thick_face, double length)
{
    switch (mode) {
    default:
    case 1: // equivalent cylinder excluded volume
        return radius_from_excluded_volume(radius, thick_rim, thick_face, length);
    case 2: // equivalent sphere
        return radius_from_volume(radius, thick_rim, thick_face, length);
    case 3: // outer rim radius
        return radius + thick_rim;
    case 4: // half outer thickness
        return 0.5*length + thick_face;
    case 5: // half diagonal
        return radius_from_diagonal(radius,thick_rim,thick_face,length);
    }
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

static void
integrand_core_shell_bicelle(double x, double q, double radius, double thick_radius, double thick_face, double halflength, double sld_core, double sld_face, double sld_rim, double sld_solvent, double* res1, double* res2, int n, int i){
    double sin_theta, cos_theta;
    get_sin_x(n, i, &sin_theta);
    get_cos_x(n, i, &cos_theta);

    const double form = bicelle_kernel(q*sin_theta, q*cos_theta, radius, thick_radius, thick_face, halflength, sld_core, sld_face, sld_rim, sld_solvent);
    *res1 = form * sin_theta;
    *res2 = form * form * sin_theta;
}

void integrate_core_shell_bicelle(
    double a,
    double b,
    double q,
    double radius,
    double thick_radius,
    double thick_face,
    double halflength,
    double sld_core,
    double sld_face,
    double sld_rim,
    double sld_solvent,
    double* res1,
    double* res2
    ){


    const double A = q*halflength;
    const double B = q*radius;
    const double C = q*(halflength+thick_face);
    const double D = q*(radius+thick_radius);

    int expo = (int)(eval_poly(log2m(max(limits[0][0],min(limits[0][1], A))), log2m(max(limits[1][0],min(limits[1][1], B))), log2m(max(limits[2][0],min(limits[2][1], C))), log2m(max(limits[3][0],min(limits[3][1], D)))) + 1);
    int n = (int)(pow(2, max(1, min(15, expo))));

    double *xg, *wg;
    get_gauss_points(n, &xg, &wg);

    // Perform the integration
    *res1 = 0;
    *res2 = 0;

    for (int i = 0; i < n; i++){
        double t1, t2;
        integrand_core_shell_bicelle(a + (b - a) * 0.5 * (xg[i] + 1), q, radius, thick_radius, thick_face, halflength, sld_core, sld_face, sld_rim, sld_solvent, &t1, &t2, n, i);
        *res1 += t1 * wg[i];
        *res2 += t2 * wg[i];
    }

    *res1 *= (b - a) * 0.5;
    *res2 *= (b - a) * 0.5;


}


static void
Fq(double q,
    double *F1,
    double *F2,
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
    const double uplim = M_PI_2;
    const double halflength = 0.5*length;
    
    double total_F1, total_F2;
    integrate_core_shell_bicelle(0, uplim, q, radius, thick_radius, thick_face, halflength, sld_core, sld_face, sld_rim, sld_solvent, &total_F1, &total_F2);


    *F1 = 1.0e-2*total_F1;
    *F2 = 1.0e-4*total_F2;


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
