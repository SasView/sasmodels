double form_volume(double radius);

double Iq(double q,
          double volfraction,
          double radius,
          double fractal_dim,
          double cor_length,
          double sld_block,
          double sld_solvent);

double Iqxy(double qx, double qy,
          double volfraction,
          double radius,
          double fractal_dim,
          double cor_length,
          double sld_block,
          double sld_solvent);


double Iq(double q,
          double volfraction,
          double radius,
          double fractal_dim,
          double cor_length,
          double sld_block,
          double sld_solvent)
{
    double qr,r0,Df,corr,phi,sldp,sldm;
    double pq,sq,inten;
    
     // Actively check the argument - needed for mass fractal - is it needie
     //here?
    if (fractal_dim <= 0.0){
       return 0.0;
    }
   
    phi = volfraction;        // volume fraction of building block spheres...
    r0 = radius;     //  radius of building block
    Df = fractal_dim;     //  fractal dimension
    corr = cor_length;       //  correlation length of fractal-like aggregates
    sldp = sld_block;       // SLD of building block
    sldm = sld_solvent;       // SLD of matrix or solution
 
     qr=q*r0;
    
    //calculate P(q) for the spherical subunits
    pq = phi*M_4PI_3*r0*r0*r0*(sldp-sldm)*(sldp-sldm)*sph_j1c(qr)*sph_j1c(qr);
    
    //calculate S(q)
    sq = Df*sas_gamma(Df-1.0)*sin((Df-1.0)*atan(q*corr));
    sq /= pow(qr,Df) * pow((1.0 + 1.0/(q*corr)/(q*corr)),((Df-1.0)/2.0));
    sq += 1.0;
    
    //combine, scale to units cm-1 sr-1 (assuming data on absolute scale)
    //and return
    inten = pq*sq;
    // convert I(1/A) to (1/cm)
    inten *= 1.0e8;
    //convert rho^2 in 10^-6 1/A to 1/A
    inten *= 1.0e-12;    
    
    
    return(inten);
}


// Iqxy is never called since no orientation or magnetic parameters.
double Iqxy(double qx, double qy,
          double volfraction,
          double radius,
          double fractal_dim,
          double cor_length,
          double sld_block,
          double sld_solvent)
{
    double q = sqrt(qx*qx + qy*qy);
    return Iq(q,
        volfraction, radius,
        fractal_dim, cor_length,
        sld_block, sld_solvent);

}

