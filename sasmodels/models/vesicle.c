double form_volume(double radius, double thickness);

double Iq(double q, 
          double sld, double sld_solvent, double volfraction,
          double radius, double thickness);

double form_volume(double radius, double thickness)
{
    //note that for the vesicle model, the volume is ONLY the shell volume
    return M_4PI_3*(cube(radius+thickness) - cube(radius));
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
    vol = M_4PI_3*cube(radius);
    f = vol * sas_3j1x_x(q*radius) * contrast;
 
    //now the shell. No volume normalization as this is done by the caller
    contrast = sld-sld_solvent;
    vol = M_4PI_3*cube(radius+thickness);
    f += vol * sas_3j1x_x(q*(radius+thickness)) * contrast;

    //rescale to [cm-1]. 
    f2 = volfraction * f*f*1.0e-4;
    
    return f2;
}
