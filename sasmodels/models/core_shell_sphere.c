double form_volume(double radius, double thickness);
double Iq(double q, double radius, double thickness, double core_sld, double shell_sld, double solvent_sld);
double Iqxy(double qx, double qy, double radius, double thickness, double core_sld, double shell_sld, double solvent_sld);


double Iq(double q, double radius, double thickness, double core_sld, double shell_sld, double solvent_sld) {
    // Core first, then add in shell
    const double core_qr = q * radius;
    const double core_contrast = core_sld - shell_sld;
    const double core_bes = sph_j1c(core_qr);
    const double core_volume = 4.0 * M_PI / 3.0 * radius * radius * radius;
    double f = core_volume * core_bes * core_contrast;

    // Now the shell
    const double shell_qr = q * (radius + thickness);
    const double shell_contrast = shell_sld - solvent_sld;
    const double shell_bes = sph_j1c(shell_qr);
    const double shell_volume = 4.0 * M_PI / 3.0 * pow((radius + thickness), 3);
    f += shell_volume * shell_bes * shell_contrast;
    return f * f * 1.0e-4;
}

double Iqxy(double qx, double qy, double radius, double thickness, double core_sld, double shell_sld, double solvent_sld) {
    const double q = sqrt(qx*qx+qy*qy);
    return Iq(q, radius, thickness, core_sld, shell_sld, solvent_sld);
}

double form_volume(double radius, double thickness)
{
    return 4.0 * M_PI / 3.0 * pow((radius + thickness), 3);
}
