static double
// following fuzzy_sphere model fuzziness is a "volume" parameter i.e. can be included in polydispersity ( dubious ???)
// so it needs to appear here. Note we assume fuzziness << radius 
form_volume(double radius, double thickness, double fuzziness)
{
    return M_4PI_3 * cube(radius + thickness);
}

static double
radius_effective(int mode, double radius, double thickness, double fuzziness)
{
    switch (mode) {
    default:
    case 1: // outer radius
        return radius + thickness;
    case 2: // core radius
        return radius;
    }
}

static void
Fq(double q, double *F1, double *F2, double radius,
   double thickness, double f_solv_core, double f_solv_shell,
   double core_sld, double shell_sld, double solvent_sld, double fuzziness
   ) {
    const double new_core_sld = f_solv_core*solvent_sld + (1.0- f_solv_core)*core_sld;
    const double new_shell_sld = f_solv_shell*solvent_sld + (1.0- f_solv_shell)*shell_sld;
    const double qf = exp(-0.5*square(q*fuzziness));
    double form = qf*core_shell_fq(q,
                              radius,
                              thickness,
                              new_core_sld,
                              new_shell_sld,
                              solvent_sld);
    *F1 = 1.0e-2*form;
    *F2 = 1.0e-4*form*form;
}
