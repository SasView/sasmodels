static double
form_volume(double length_a, double length_b, double length_c)
{
    return length_a * length_b * length_c;
}

static double
radius_from_excluded_volume(double length_a, double length_b, double length_c)
{
    double r_equiv, length;
    double lengths[3] = {length_a, length_b, length_c};
    double lengthmax = fmax(lengths[0],fmax(lengths[1],lengths[2]));
    double length_1 = lengthmax;
    double length_2 = lengthmax;
    double length_3 = lengthmax;

    for(int ilen=0; ilen<3; ilen++) {
        if (lengths[ilen] < length_1) {
            length_2 = length_1;
            length_1 = lengths[ilen];
            } else {
                if (lengths[ilen] < length_2) {
                        length_2 = lengths[ilen];
                }
            }
    }
    if(length_2-length_1 > length_3-length_2) {
        r_equiv = sqrt(length_2*length_3/M_PI);
        length  = length_1;
    } else  {
        r_equiv = sqrt(length_1*length_2/M_PI);
        length  = length_3;
    }

    return 0.5*cbrt(0.75*r_equiv*(2.0*r_equiv*length + (r_equiv + length)*(M_PI*r_equiv + length)));
}

static double
radius_effective(int mode, double length_a, double length_b, double length_c)
{
    switch (mode) {
    default:
    case 1: // equivalent cylinder excluded volume
        return radius_from_excluded_volume(length_a,length_b,length_c);
    case 2: // equivalent volume sphere
        return cbrt(length_a*length_b*length_c/M_4PI_3);
    case 3: // half length_a
        return 0.5 * length_a;
    case 4: // half length_b
        return 0.5 * length_b;
    case 5: // half length_c
        return 0.5 * length_c;
    case 6: // equivalent circular cross-section
        return sqrt(length_a*length_b/M_PI);
    case 7: // half ab diagonal
        return 0.5*sqrt(length_a*length_a + length_b*length_b);
    case 8: // half diagonal
        return 0.5*sqrt(length_a*length_a + length_b*length_b + length_c*length_c);
    }
}

static void
Fq(double q,
    double *F1,
    double *F2,
    double sld,
    double solvent_sld,
    double length_a,
    double length_b,
    double length_c)
{
    const double mu = 0.5 * q * length_b;

    // Scale sides by B
    const double a_scaled = length_a / length_b;
    const double c_scaled = length_c / length_b;

    // outer integral (with gauss points), integration limits = 0, 1
    double outer_total_F1 = 0.0; //initialize integral
    double outer_total_F2 = 0.0; //initialize integral
    for( int i=0; i<GAUSS_N; i++) {
        const double sigma = 0.5 * ( GAUSS_Z[i] + 1.0 );
        const double mu_proj = mu * sqrt(1.0-sigma*sigma);

        // inner integral (with gauss points), integration limits = 0, 1
        // corresponding to angles from 0 to pi/2.
        double inner_total_F1 = 0.0;
        double inner_total_F2 = 0.0;
        for(int j=0; j<GAUSS_N; j++) {
            const double uu = 0.5 * ( GAUSS_Z[j] + 1.0 );
            double sin_uu, cos_uu;
            SINCOS(M_PI_2*uu, sin_uu, cos_uu);
            const double si1 = sas_sinx_x(mu_proj * sin_uu * a_scaled);
            const double si2 = sas_sinx_x(mu_proj * cos_uu);
            const double fq = si1 * si2;
            inner_total_F1 += GAUSS_W[j] * fq;
            inner_total_F2 += GAUSS_W[j] * fq * fq;
        }
        // now complete change of inner integration variable (1-0)/(1-(-1))= 0.5
        inner_total_F1 *= 0.5;
        inner_total_F2 *= 0.5;

        const double si = sas_sinx_x(mu * c_scaled * sigma);
        outer_total_F1 += GAUSS_W[i] * inner_total_F1 * si;
        outer_total_F2 += GAUSS_W[i] * inner_total_F2 * si * si;
    }
    // now complete change of outer integration variable (1-0)/(1-(-1))= 0.5
    outer_total_F1 *= 0.5;
    outer_total_F2 *= 0.5;

    // Multiply by contrast^2 and convert from [1e-12 A-1] to [cm-1]
    const double V = form_volume(length_a, length_b, length_c);
    const double contrast = (sld-solvent_sld);
    const double s = contrast * V;
    *F1 = 1.0e-2 * s * outer_total_F1;
    *F2 = 1.0e-4 * s * s * outer_total_F2;
}

static double
Iqabc(double qa, double qb, double qc,
    double sld,
    double solvent_sld,
    double mag_sld,
    double length_a,
    double length_b,
    double length_c,
    double Length, double Singlets, double Doubles, double Trimers, double Quadramers, double Pentamers)
{
    const double siA = sas_sinx_x(0.5*length_a*qa);
    const double siB = sas_sinx_x(0.5*length_b*qb);
    const double siC = sas_sinx_x(0.5*length_c*qc);
    const double V = form_volume(length_a, length_b, length_c);
    const double drho = (sld - solvent_sld);
    const double mrho = mag_sld
    const double Amp = V * drho * siA * siB * siC;
    const double MAmp = V * mrho * siA * siB * siC;

    // adding several particles in the chain:
    double VolumeFraction = 1.0;
    double SingletFraction = Singlets;
    double DimerFraction = Doubles;
    double TrimerFraction = Trimers;
    double QuadramerFraction = Quadramers;
    double PentamerFraction = Pentamers;

    // leave out polydispersity and angular orientation

    // Just start parameters:
    double SingletIntensity = 0;
    double MSingletIntensity = 0;
    double DimerIntensity = 0;
    double MDimerIntensity = 0;
    double TrimerIntensity = 0;
    double MTrimerIntensity = 0;
    double QuadramerIntensity = 0;
    double MQuadramerIntensity = 0;
    double PentamerIntensity = 0;
    double MPentamerIntensity = 0;

    // we take out: (i) polydispersity of chains (called "anglewt/sigma" before), (ii) orientation of chain (called "Viewing angle")
    // Now: chain only oriented along Qx
    // What does it mean for the magnetic orientation of the particles?: at the moment not oriented / random

    double Q_X = q*cos(0);
    double Q_Y = q*sin(0);


    double ChainProjX = cos(0);
    double ChainProjY = cos(phi);
    double ChainProjZ = sin(phi);


    double real_phase = 1.0;
    double img_phase = 0.0;
    double mreal_phase = 1.0;
    double mimg_phase = 0.0;

    for(int k=1; k<5; k++){
        real_phase += cos(k*Length*(Q_X*ChainProjX + Q_Y*ChainProjY));
        img_phase += sin(k*Length*(Q_X*ChainProjX + Q_Y*ChainProjY));

    if(k==1){
    DimerIntensity  += pow(Amp*img_phase,2))/(2.0*Vol);
    MDimerIntensity  += pow(MAmp*mimg_phase,2))/(2.0*Vol);
    }
    if(k==2){
    TrimerIntensity  += pow(Amp*img_phase,2))/(3.0*Vol);
    MTrimerIntensity  += pow(MAmp*mimg_phase,2))/(3.0*Vol);
    }
    if(k==3){
    QuadramerIntensity  += pow(Amp*img_phase,2))/(4.0*Vol);
    MQuadramerIntensity  += pow(MAmp*mimg_phase,2))/(4.0*Vol);
    }
    if(k==4){
    PentamerIntensity  += pow(Amp*img_phase,2))/(5.0*Vol);
    MPentamerIntensity  += pow(MAmp*mimg_phase,2))/(5.0*Vol);
    }
    }
    //end k loop for dimers

    double FractionScale = SingletFraction + DimerFraction + TrimerFraction + QuadramerFraction + PentamerFraction;
    if(FractionScale == 0){FractionScale = 1.0;}

    double SIntensity = SingletFraction*SingletIntensity + DimerFraction*DimerIntensity + TrimerFraction*TrimerIntensity + QuadramerFraction*QuadramerIntensity + PentamerFraction*PentamerIntensity;
    double MIntensity = 0.0;
    if(MVar <= 1){
        MIntensity = MSingletIntensity*(SingletFraction + DimerFraction + TrimerFraction + QuadramerFraction + PentamerFraction);
    }
    else{
        MIntensity = SingletFraction*MSingletIntensity + DimerFraction*MDimerIntensity + TrimerFraction*MTrimerIntensity + QuadramerFraction*MQuadramerIntensity + PentamerFraction*MPentamerIntensity;
    }

    double Intensity = (SIntensity+MIntensity)*(1E4)/FractionScale;

    return Intensity;

}

