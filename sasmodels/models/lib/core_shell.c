/*******************************************************************

core_shell_kernel

Form factor used in core_shell and fractal_core_shell

********************************************************************/
static
double core_shell_kernel(double q,
                         double radius,
                         double thickness,
                         double core_sld,
                         double shell_sld,
                         double solvent_sld) {
    // Core first, then add in shell
    const double core_qr = q * radius;
    const double core_contrast = core_sld - shell_sld;
    const double core_bes = sas_3j1x_x(core_qr);
    const double core_volume = M_4PI_3 * cube(radius);
    double f = core_volume * core_bes * core_contrast;

    // Now the shell
    const double shell_qr = q * (radius + thickness);
    const double shell_contrast = shell_sld - solvent_sld;
    const double shell_bes = sas_3j1x_x(shell_qr);
    const double shell_volume = M_4PI_3 * cube(radius + thickness);
    f += shell_volume * shell_bes * shell_contrast;
    return f * f * 1.0e-4;
}