// NOTE that "length" here is the full height of the core!
static double
form_volume(double r_minor,
        double x_core,
        double thick_rim,
        double thick_face,
        double length,
        double radius_protein,
        double x_protein,
        double thick_protein)
{
    return M_PI*(  (r_minor + thick_rim)*(r_minor*x_core + thick_rim)* length +
                 square(r_minor)*x_core*2.0*thick_face +
                 square(radius_protein)*x_protein*2.0*thick_protein );
}

static double
radius_from_excluded_volume(double r_minor, double x_core, double thick_rim, double thick_face, double length,
 double radius_protein, double x_protein, double thick_protein)
{
    const double r_equiv     = sqrt((r_minor + thick_rim)*(r_minor*x_core + thick_rim));
    const double length_tot  = length + 2.0*thick_face +2.0*thick_protein;
    return 0.5*cbrt(0.75*r_equiv*(2.0*r_equiv*length_tot + (r_equiv + length_tot)*(M_PI*r_equiv + length_tot)));
}

static double
radius_from_volume(double r_minor, double x_core, double thick_rim, double thick_face, double length,
 double radius_protein, double x_protein, double thick_protein)
{
    const double volume_bicelle = form_volume(r_minor, x_core, thick_rim,thick_face,length, radius_protein,
    x_protein, thick_protein);
    return cbrt(volume_bicelle/M_4PI_3);
}

static double
radius_from_diagonal(double r_minor, double x_core, double thick_rim, double thick_face, double length,
 double radius_protein, double x_protein, double thick_protein)
 // RKH SORT THIS - do we want to allow for diagonal of protein ???
{
    const double radius_max = (x_core < 1.0 ? r_minor : x_core*r_minor);
    const double radius_max_tot = radius_max + thick_rim;
    const double length_tot = length + 2.0*thick_face;
    return sqrt(radius_max_tot*radius_max_tot + 0.25*length_tot*length_tot);
}

static double
radius_effective(int mode, double r_minor, double x_core, double thick_rim, double thick_face, double length,
 double radius_protein, double x_protein, double thick_protein)
{
    switch (mode) {
    default:
    case 1: // equivalent cylinder excluded volume
        return radius_from_excluded_volume(r_minor, x_core, thick_rim, thick_face,
                                length, radius_protein, x_protein, thick_protein);
    case 2: // equivalent sphere
        return radius_from_volume(r_minor, x_core, thick_rim, thick_face,
                                length, radius_protein, x_protein, thick_protein);
    case 3: // outer rim average radius
        return 0.5*r_minor*(1.0 + x_core) + thick_rim;
    case 4: // outer rim min radius
        return (x_core < 1.0 ? x_core*r_minor+thick_rim : r_minor+thick_rim);
    case 5: // outer max radius
        return (x_core > 1.0 ? x_core*r_minor+thick_rim : r_minor+thick_rim);
    case 6: // half outer thickness
        return 0.5*length + thick_face;
    case 7: // half diagonal (this ignores the missing "corners", so may give unexpected answer if thick_face
            //  or thick_rim is extremely large, ALSO ignores the protein)
        return radius_from_diagonal(r_minor,x_core,thick_rim,thick_face,
                                length, radius_protein, x_protein, thick_protein);
    }
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
        double rhoc,
        double rhoh,
        double rhor,
        double rhosolv,
        double sigma,
        double r_prot_minor,
        double x_prot,
        double thick_prot,
        double rhoprot)
{
      // core_shell_bicelle_elliptical_protein, RKH 24th April 2019, based on ..
     // core_shell_bicelle_elliptical_belt_rough, RKH 5th Oct 2017, based on core_shell_bicelle_elliptical
    const double halfheight = 0.5*length;
    const double r_major = r_minor * x_core;
    const double r2A = 0.5*(square(r_major) + square(r_minor));
    const double r2B = 0.5*(square(r_major) - square(r_minor));
    const double vol1 = M_PI*r_minor*r_major*(2.0*halfheight);
    const double vol2 = M_PI*(r_minor+thick_rim)*(r_major+thick_rim)*2.0*halfheight;
    const double vol3 = M_PI*r_minor*r_major*2.0*(halfheight+thick_face);
    // dr1,2,3 are now for Vcore, Vcore+rim, Vcore+face,
    const double dr1 = vol1*(-rhor - rhoh + rhoc + rhosolv);
    const double dr2 = vol2*(rhor-rhosolv);
    const double dr3 = vol3*(rhoh-rhosolv);
    // three more terms for the protein, which can stick out of the face (on both sides simultaneously)
    const double r_prot_major = r_prot_minor * x_prot;
    const double r2A_prot = 0.5*(square(r_prot_major) + square(r_prot_minor));
    const double r2B_prot = 0.5*(square(r_prot_major) - square(r_prot_minor));
    const double area_prot = M_PI*r_prot_minor*r_prot_major;
    const double dr4 = area_prot*2.0*(halfheight + thick_face + thick_prot)*(rhoprot - rhosolv);
    // note this next term really does need (-rhoh +rhosolv)
    const double dr5 = area_prot*2.0*(halfheight + thick_face)*(rhosolv -rhoh);
    const double dr6 = area_prot*2.0*halfheight*(rhoh-rhoc);

    //initialize integral
    double outer_total_F1 = 0.0;
    double outer_total_F2 = 0.0;
    for(int i=0;i<GAUSS_N;i++) {
        //setup inner integral over the ellipsoidal cross-section
        // since we generate these lots of times, why not store them somewhere?
        //const double cos_theta = ( GAUSS_Z[i]*(vb-va) + va + vb )/2.0;
        const double cos_theta = ( GAUSS_Z[i] + 1.0 )/2.0;
        const double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
        const double qab = q*sin_theta;
        const double qc = q*cos_theta;
        const double si1 = sas_sinx_x(halfheight*qc);
        const double si2 = sas_sinx_x((halfheight + thick_face)*qc);
        const double si3 = sas_sinx_x((halfheight + thick_face + thick_prot)*qc);
        double inner_total_F1 = 0;
        double inner_total_F2 = 0;
        for(int j=0;j<GAUSS_N;j++) {
            //76 gauss points for the inner integral (WAS 20 points,so this may make unnecessarily slow, but playing safe)
            //const double beta = ( GAUSS_Z[j]*(vbj-vaj) + vaj + vbj )/2.0;
            const double beta = ( GAUSS_Z[j] +1.0)*M_PI_2;
            const double rr = sqrt(r2A - r2B*cos(beta));
            const double be1 = sas_2J1x_x(rr*qab);
            const double be2 = sas_2J1x_x((rr+thick_rim)*qab);
            const double rrp = sqrt(r2A_prot - r2B_prot*cos(beta));
            const double be4 = sas_2J1x_x(rrp*qab);
            const double f = dr1*si1*be1 + dr2*si1*be2 + dr3*si2*be1 +
                             be4*( dr4*si3 + dr5*si2 +dr6*si1);

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
    *F1 = 1e-2*outer_total_F1*exp(-0.25*square(q*sigma));
    *F2 = 1e-4*outer_total_F2*exp(-0.5*square(q*sigma));
}

static double
Iqabc(double qa, double qb, double qc,
          double r_minor,
          double x_core,
          double thick_rim,
          double thick_face,
          double length,
          double rhoc,
          double rhoh,
          double rhor,
          double rhosolv,
          double sigma,
          double r_prot_minor,
          double x_prot,
          double thick_prot,
          double rhoprot)
{
    // Vol1,2,3 and dr1,2,3 are now for Vcore, Vcore+rim, Vcore+face,
    const double dr1 = -rhor - rhoh + rhoc + rhosolv;
    const double dr2 = rhor-rhosolv;
    const double dr3 = rhoh-rhosolv;
    const double r_major = r_minor*x_core;
    const double halfheight = 0.5*length;
    const double vol1 = M_PI*r_minor*r_major*length;
    const double vol2 = M_PI*(r_minor+thick_rim)*(r_major+thick_rim)*2.0*halfheight;
    const double vol3 = M_PI*r_minor*r_major*2.0*(halfheight+thick_face);
    // three more terms for the protein
    const double r_prot_major = r_prot_minor * x_prot;
    const double area_prot = M_PI*r_prot_minor*r_prot_major;
    const double dr4 = area_prot*2.0*(halfheight + thick_face + thick_prot)*(rhoprot - rhosolv);
    // note this next term really does need (-rhoh +rhosolv)
    const double dr5 = area_prot*2.0*(halfheight + thick_face)*(rhosolv -rhoh);
    const double dr6 = area_prot*2.0*halfheight*(rhoh-rhoc);

    // Compute effective radius in rotated coordinates
    const double qr_hat = sqrt(square(r_major*qb) + square(r_minor*qa));
    // does this need to be changed for the "missing corners" where there there is no "belt" ?
    const double qrshell_hat = sqrt(square((r_major+thick_rim)*qb)
                                   + square((r_minor+thick_rim)*qa));
    const double qr_prot_hat = sqrt(square(r_prot_major*qb) + square(r_prot_minor*qa));
    const double be1 = sas_2J1x_x( qr_hat );
    const double be2 = sas_2J1x_x( qrshell_hat );
    const double be4 = sas_2J1x_x( qr_prot_hat );
    const double si1 = sas_sinx_x( halfheight*qc );
    const double si2 = sas_sinx_x( (halfheight + thick_face)*qc );
    const double si3 = sas_sinx_x( (halfheight + thick_face + thick_prot)*qc );
    const double fq = vol1*dr1*si1*be1 + vol2*dr2*si1*be2 +  vol3*dr3*si2*be1 +
                           be4*( dr4*si3 + dr5*si2 +dr6*si1) ;
    const double atten = exp(-0.5*(qa*qa + qb*qb + qc*qc)*(sigma*sigma));
    return 1.0e-4 * fq*fq * atten;
}
