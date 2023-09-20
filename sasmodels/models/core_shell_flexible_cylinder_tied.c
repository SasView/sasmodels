static double
// compiler needs ALL params here to by of type "volume" in the .py definition, even though we may not really want to be 
// able to makle them all polydisperse.
form_volume(double length, 
    double kuhn_length,
    double radius,
    double vol_dry_shell_over_core,
    double f_solvent_in_shell)
{
// I think that because we are normalising to CORE volume here, that we can get away with returning 1.0 for the volume,
// compare core_shell_ellipsoid_tied where I solved for the shell thickness in form_volume
   return 1.0;
}

static double
Iq(double q,
   double length,
   double kuhn_length,
   double radius,
   double vol_dry_shell_over_core,
   double sld_core,
   double sld_dry_shell,
   double sld_solvent,
   double f_solvent_in_shell)
{
    const double sld_shell = f_solvent_in_shell*sld_solvent + (1-f_solvent_in_shell)*sld_dry_shell;
    const double contrast1 = sld_core - sld_shell;
    const double contrast2 = sld_shell - sld_solvent;
    const double thickness = radius*( sqrt( vol_dry_shell_over_core/(1-f_solvent_in_shell) + 1.0 ) -1);
    //printf(" thickness = %g\n",thickness);
    const double radius2 = radius + thickness;
    const double volume1 = M_PI*radius*radius*length;
    const double volume2 = M_PI*radius2*radius2*length;
    const double flex = Sk_WR(q, length, kuhn_length);
    // note normalised to core volume here,  if contrast2 =0 should give same result as flexible_cylinder
    return 1.0e-4 *square(volume1*contrast1*sas_2J1x_x(q*radius) + volume2*contrast2*sas_2J1x_x(q*radius2)) * flex / volume1;
}
