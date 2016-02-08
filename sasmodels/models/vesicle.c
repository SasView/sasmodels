double form_volume(double radius, double thickness);

double Iq(double q, 
          double sld, double solvent_sld,
          double radius, double thickness);

double Iqxy(double qx, double qy,
          double sld, double solvent_sld,
          double radius, double thickness);

double form_volume(double radius, double thickness)
{
    //note that for the vesicle model, the volume is ONLY the shell volume
    double volume;
    volume =4.*M_PI*(radius+thickness)*(radius+thickness)*(radius+thickness)/3;
    volume -=4.*M_PI*radius*radius*radius/3.;
    return volume;
}

double Iq(double q,
    double sld,
    double solvent_sld,
    double radius,
    double thickness)

/*
   scattering from a unilamellar vesicle.
   same functional form as the core-shell sphere, but more intuitive for
   a vesicle
*/

/*
   note that the sph_j1c we are using has been optimized for precision over
   SasView's original implementation. HOWEVER at q==0 that implementation
   set bes=1.0 rather than 0.0 (correct value) on the grounds I believe that 
   bes=0.00 causes Iq to have a divide by 0 error (mostly encountered when
   doing a theory curve in 2D?  We should verify this and if necessary fix
     -PDB Feb 7, 2016 
*/
{
    double bes,vol,contrast,f,f2;

    // core first, then add in shell
    contrast = solvent_sld-sld;
    bes = sph_j1c(q*radius);
    vol = 4.0*M_PI/3.0*radius*radius*radius;
    f = vol*bes*contrast;
 
    //now the shell
    contrast = sld-solvent_sld;
    bes = sph_j1c(q*(radius+thickness));
    vol = 4.0*M_PI/3.0*(radius+thickness)*(radius+thickness)*(radius+thickness);
    f += vol*bes*contrast;

    //rescale to [cm-1]. No volume normalization as this is done by the caller
    f2 = f*f*1.0e-4;
    
    return(f2);
}


double Iqxy(double qx, double qy,
          double sld, double solvent_sld,
          double radius, double thickness)
          
{
    double q = sqrt(qx*qx + qy*qy);
    return Iq(q,
        sld, solvent_sld,
        radius,thickness);

}
