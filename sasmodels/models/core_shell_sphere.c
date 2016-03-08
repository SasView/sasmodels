double form_volume(double radius, double thickness);
double Iq(double q, double radius, double thickness, double core_sld, double shell_sld, double solvent_sld);
double Iqxy(double qx, double qy, double radius, double thickness, double core_sld, double shell_sld, double solvent_sld);


double Iq(double q, double radius, double thickness, double core_sld, double shell_sld, double solvent_sld) {


    double intensity = core_shell_kernel(q,
                              radius,
                              thickness,
                              core_sld,
                              shell_sld,
                              solvent_sld);
    return intensity;
}

double Iqxy(double qx, double qy, double radius, double thickness, double core_sld, double shell_sld, double solvent_sld) {
    const double q = sqrt(qx*qx+qy*qy);
    return Iq(q, radius, thickness, core_sld, shell_sld, solvent_sld);
}

double form_volume(double radius, double thickness)
{
    return 4.0 * M_PI / 3.0 * pow((radius + thickness), 3);
}
