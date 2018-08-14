// NOTE that "length" here is the full height of the core!
static double
form_volume(double r_minor,
    double x_core,
    double thick_rim,
    double thick_face,
    double length)
{
    return M_PI*(r_minor+thick_rim)*(r_minor*x_core+thick_rim)*(length+2.0*thick_face);
}

static void
Fq(double q,
    double *F1,
    double *F2,
    double r_minor,
    double x_core,
    double thick_rim,
    double thick_face,
    double length,
    double sld_core,
    double sld_face,
    double sld_rim,
    double sld_solvent)
{
     // core_shell_bicelle_elliptical, RKH Dec 2016, based on elliptical_cylinder and core_shell_bicelle
     // tested against limiting cases of cylinder, elliptical_cylinder, stacked_discs, and core_shell_bicelle
    const double halfheight = 0.5*length;
    const double r_major = r_minor * x_core;
    const double r2A = 0.5*(square(r_major) + square(r_minor));
    const double r2B = 0.5*(square(r_major) - square(r_minor));
    const double vol1 = M_PI*r_minor*r_major*(2.0*halfheight);
    const double vol2 = M_PI*(r_minor+thick_rim)*(r_major+thick_rim)*2.0*(halfheight+thick_face);
    const double vol3 = M_PI*r_minor*r_major*2.0*(halfheight+thick_face);
    const double dr1 = vol1*(sld_core-sld_face);
    const double dr2 = vol2*(sld_rim-sld_solvent);
    const double dr3 = vol3*(sld_face-sld_rim);

    //initialize integral
    double outer_total_F1 = 0.0;
    double outer_total_F2 = 0.0;
    for(int i=0;i<GAUSS_N;i++) {
        //setup inner integral over the ellipsoidal cross-section
        //const double cos_theta = ( GAUSS_Z[i]*(vb-va) + va + vb )/2.0;
        const double cos_theta = ( GAUSS_Z[i] + 1.0 )/2.0;
        const double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
        const double qab = q*sin_theta;
        const double qc = q*cos_theta;
        const double si1 = sas_sinx_x(halfheight*qc);
        const double si2 = sas_sinx_x((halfheight+thick_face)*qc);
        double inner_total_F1 = 0;
        double inner_total_F2 = 0;
        for(int j=0;j<GAUSS_N;j++) {
            //76 gauss points for the inner integral (WAS 20 points,so this may make unecessarily slow, but playing safe)
            //const double beta = ( GAUSS_Z[j]*(vbj-vaj) + vaj + vbj )/2.0;
            const double beta = ( GAUSS_Z[j] +1.0)*M_PI_2;
            const double rr = sqrt(r2A - r2B*cos(beta));
            const double be1 = sas_2J1x_x(rr*qab);
            const double be2 = sas_2J1x_x((rr+thick_rim)*qab);
            const double f = dr1*si1*be1 + dr2*si2*be2 + dr3*si2*be1;

            inner_total_F1 += GAUSS_W[j] * f;
            inner_total_F2 += GAUSS_W[j] * f * f;
        }
        //now calculate outer integral
        outer_total_F1 += GAUSS_W[i] * inner_total_F1;
        outer_total_F2 += GAUSS_W[i] * inner_total_F2;
    }
    // now complete change of integration variables (1-0)/(1-(-1))= 0.5
    outer_total_F1 *= 0.25;
    outer_total_F2 *= 0.25;

    //convert from [1e-12 A-1] to [cm-1]
    *F1 = 1e-2*outer_total_F1;
    *F2 = 1e-4*outer_total_F2;
}

static double
Iqabc(double qa, double qb, double qc,
    double r_minor,
    double x_core,
    double thick_rim,
    double thick_face,
    double length,
    double sld_core,
    double sld_face,
    double sld_rim,
    double sld_solvent)
{
    const double dr1 = sld_core-sld_face;
    const double dr2 = sld_rim-sld_solvent;
    const double dr3 = sld_face-sld_rim;
    const double r_major = r_minor*x_core;
    const double halfheight = 0.5*length;
    const double vol1 = M_PI*r_minor*r_major*length;
    const double vol2 = M_PI*(r_minor+thick_rim)*(r_major+thick_rim)*2.0*(halfheight+thick_face);
    const double vol3 = M_PI*r_minor*r_major*2.0*(halfheight+thick_face);

    // Compute effective radius in rotated coordinates
    const double qr_hat = sqrt(square(r_major*qb) + square(r_minor*qa));
    const double qrshell_hat = sqrt(square((r_major+thick_rim)*qb)
                                   + square((r_minor+thick_rim)*qa));
    const double be1 = sas_2J1x_x( qr_hat );
    const double be2 = sas_2J1x_x( qrshell_hat );
    const double si1 = sas_sinx_x( halfheight*qc );
    const double si2 = sas_sinx_x( (halfheight + thick_face)*qc );
    const double fq = vol1*dr1*si1*be1 + vol2*dr2*si2*be2 +  vol3*dr3*si2*be1;
    return 1.0e-4 * fq*fq;
}
