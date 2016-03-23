double form_volume(double radius, double thickness);

double Iq(double q, 
          double sld, double sld_solvent, double volfraction,
          double radius, double thickness);

double Iqxy(double qx, double qy,
          double sld, double sld_solvent, double volfraction,
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
    double sld_solvent,
    double volfraction,
    double radius,
    double thickness)

/*
   scattering from a unilamellar vesicle.
   same functional form as the core-shell sphere, but more intuitive for
   a vesicle
*/

{
    double vol,contrast,f,f2;

    // core first, then add in shell
    contrast = sld_solvent-sld;
    vol = 4.0*M_PI/3.0*radius*radius*radius;
    f = vol*sph_j1c(q*radius)*contrast;
 
    //now the shell. No volume normalization as this is done by the caller
    contrast = sld-sld_solvent;
    vol = 4.0*M_PI/3.0*(radius+thickness)*(radius+thickness)*(radius+thickness);
    f += vol*sph_j1c(q*(radius+thickness))*contrast;

    //rescale to [cm-1]. 
    f2 = volfraction*f*f*1.0e-4;
    
    return(f2);
}


double Iqxy(double qx, double qy,
          double sld, double sld_solvent, double volfraction,
          double radius, double thickness)
          
{
    double q = sqrt(qx*qx + qy*qy);
    return Iq(q,
        sld, sld_solvent, volfraction,
        radius,thickness);

}
